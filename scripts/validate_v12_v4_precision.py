#!/usr/bin/env python3
"""Small-cost V4 calibration against block-matrix event-count expectations.

Defaults are a smoke experiment, not the full V4 acceptance protocol. Reports
actual fixed-count and adaptive-stop coverage with binomial uncertainty, and
fixed-tree versus independently resampled tree-population standard errors.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears._ffi import PreparedSampler
from rustbears.uncertainty.checkpoint import model_manifest,task_fingerprint,canonical_manifest
from rustbears.uncertainty.diagnostics import run_bsm_adaptive
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
from rustbears.uncertainty.simulate import sampler_buffers
from scripts.validate_v12_gate1 import ReferenceTree


def exact_event_counts(model,params):
    """Frechet block exponential on each branch, independent dense pruning."""
    reference=ReferenceTree(model,params,root_weights=model.root_sampling_probs)
    likelihood=float(reference.root_weights @ reference.downpass(model.tips)[reference.root])
    n=len(reference.states)
    answer={}
    for kind,increment in (("d",1),("e",-1)):
        reward=reference.q.copy()
        for i,source in enumerate(reference.states):
            for j,target in enumerate(reference.states):
                if i==j or len(target)-len(source)!=increment:
                    reward[i,j]=0
        expected=0.
        for child,transition in reference.transition.items():
            t=model.tree.lengths[child]
            block=np.block([[reference.q,reward],[np.zeros_like(reward),reference.q]])
            moment=expm(block*t)[:n,n:]
            reference.transition[child]=moment
            expected+=float(reference.root_weights @ reference.downpass(model.tips)[reference.root])/likelihood
            reference.transition[child]=transition
        answer[kind]=expected
    return answer


def make_fit(scale=1):
    tree=b.read_newick("((T1:1,T2:1):1,(T3:1.3,T4:1.3):.7);")
    tree.lengths*=scale
    tree.ages*=scale
    model=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=.3),e=dict(type="fixed",init=.2))
    return b.BioGeoBEARS(tree,{"T1":"100","T2":"010","T3":"101","T4":"001"},model,
                        likelihood_conditioning="tip_survival").fit()


def counts(sampler,n,seed):
    raw=sampler.draw(1,n,seed,1,bsm=True)
    table=raw[1].reshape((-1,10))
    return {kind:np.bincount(table[table[:,7]==code,0].astype(int)-1,minlength=n)
            for kind,code in (("d",0),("e",1))}


def report(rows,truth):
    means=np.array([r[0] for r in rows]);ses=np.array([r[1] for r in rows])
    cover=(means-1.96*ses<=truth)&(means+1.96*ses>=truth)
    hit,n=int(cover.sum()),len(cover)
    ci=[0. if hit==0 else float(beta.ppf(.025,hit,n-hit+1)),
        1. if hit==n else float(beta.ppf(.975,hit+1,n-hit))]
    return dict(truth=truth,replicates=n,mean_estimate=float(means.mean()),bias=float(means.mean()-truth),
                empirical_sd=float(means.std(ddof=1)),mean_estimated_se=float(ses.mean()),
                covered=hit,coverage=hit/n,coverage_cp95=ci)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True)
    parser.add_argument("--replicates",type=int,default=20)
    parser.add_argument("--maps",type=int,default=200)
    parser.add_argument("--maps-per-tree",type=int,default=20)
    parser.add_argument("--adaptive-max",type=int,default=1000)
    parser.add_argument("--seed",type=int,default=20261115)
    args=parser.parse_args()
    if min(args.replicates,args.maps,args.maps_per_tree,args.adaptive_max)<2:
        parser.error("all experiment counts must be at least two")
    started=time.monotonic()
    fits=[make_fit(scale) for scale in (.6,.8,1.,1.2,1.4)]
    truths=[exact_event_counts(f.model,f.params) for f in fits]
    tasks=[TaskIdentity("v4_precision","DEC",str(tree),outer_rep=r,role=phase)
           for r in range(1,args.replicates+1) for tree in range(5)
           for phase in ("fixed","fixed_trees","population","adaptive")]
    seeds=derive_task_seeds(args.seed,tasks)
    def seed(r,tree,phase):
        return seeds[TaskIdentity("v4_precision","DEC",str(tree),outer_rep=r,role=phase)]
    fixed={x:[] for x in "de"};adaptive={x:[] for x in "de"}
    fixed_trees={x:[] for x in "de"};population={x:[] for x in "de"}
    stops=[]
    samplers=[PreparedSampler(sampler_buffers(f.model,f.params),f.model.root_sampling_probs,0) for f in fits]
    manifest=dict(protocol="V4-small-calibration-v1",settings=vars(args),models=[model_manifest(f) for f in fits])
    try:
        for r in range(1,args.replicates+1):
            data=counts(samplers[2],args.maps,seed(r,2,"fixed"))
            for kind in "de":
                fixed[kind].append((data[kind].mean(),data[kind].std(ddof=1)/np.sqrt(args.maps)))
            sampled=[counts(s,args.maps_per_tree,seed(r,i,"fixed_trees")) for i,s in enumerate(samplers)]
            for kind in "de":
                means=np.array([d[kind].mean() for d in sampled])
                variances=np.array([d[kind].var(ddof=1) for d in sampled])
                fixed_trees[kind].append((means.mean(),np.sqrt(variances.sum()/(25*args.maps_per_tree))))
            draw=np.random.default_rng(seed(r,0,"population")).integers(0,5,size=5)
            sampled=[counts(samplers[int(tree)],args.maps_per_tree,seed(r,i,"population")) for i,tree in enumerate(draw)]
            for kind in "de":
                means=np.array([d[kind].mean() for d in sampled])
                population[kind].append((means.mean(),means.std(ddof=1)/np.sqrt(5)))
            mapped=run_bsm_adaptive(fits[2],min_maps=min(100,args.adaptive_max),batch_size=100,
                    max_maps=args.adaptive_max,target_rmcse=.1,min_mean_count=.5,seed=seed(r,2,"adaptive"))
            stops.append(dict(replicate=r,maps=mapped.num_maps,reason=mapped.stop_reason))
            for kind in "de":
                row=mapped.diagnostics.event_summary_df.set_index("event_type").loc[kind]
                adaptive[kind].append((row["mean"],row.se_target))
    finally:
        for sampler in samplers:
            sampler.close()
    result=dict(schema="rustbears.v4.precision-calibration.v1",status="experiments_completed",
                task_fingerprint=task_fingerprint(manifest),settings=vars(args),
                manifest=canonical_manifest(manifest),
                fixed_count_400_by_500_completed=(args.replicates>=400 and args.maps>=500),
                full_v4_acceptance_protocol_completed=False,
                reference_boundary="Independent DEC Q and dense pruning; block matrix exponential count moments; established cladogenetic local weight kernel reused.",
                interpretation="Empirical calibration with binomial uncertainty; these intervals are nominal, not anytime-valid. This small experiment does not replace the complete V4 protocol.",
                exact_means_by_tree=truths,
                fixed_count={x:report(fixed[x],truths[2][x]) for x in "de"},
                adaptive_stop={x:report(adaptive[x],truths[2][x]) for x in "de"},
                fixed_tree_sample={x:report(fixed_trees[x],float(np.mean([t[x] for t in truths]))) for x in "de"},
                independent_tree_population={x:report(population[x],float(np.mean([t[x] for t in truths]))) for x in "de"},
                stopping_records=stops,seconds=time.monotonic()-started,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    destination=Path(args.output);destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":result["status"],"seconds":result["seconds"],"output":str(destination),
                      "fixed_coverage":{x:result["fixed_count"][x]["coverage"] for x in "de"},
                      "adaptive_coverage":{x:result["adaptive_stop"][x]["coverage"] for x in "de"}}))


if __name__=="__main__":
    main()
