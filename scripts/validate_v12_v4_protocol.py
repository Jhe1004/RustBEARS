#!/usr/bin/env python3
"""Frozen full V4 experiment protocol; default writes a plan and runs no grid.

Use --mode smoke for a deliberately reduced execution. --mode formal must be
scheduled only after the intended Gate2/AMD release. Results are descriptive
calibration evidence, never rerun-until-pass acceptance decisions.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears._ffi import PreparedSampler
from rustbears.uncertainty.checkpoint import (CheckpointStore,model_manifest,
    canonical_manifest,task_fingerprint,frame_arrays)
from rustbears.uncertainty.diagnostics import _Design,rank_status,run_bsm_adaptive
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
from rustbears.uncertainty.simulate import sampler_buffers
from scripts.validate_v12_gate1 import ReferenceTree
from scripts.validate_v12_v4_precision import make_fit,report


FORMAL=dict(fixed=dict(replicates=400,maps=500),
    fixed_trees=dict(replicates=500,trees=5,maps_per_tree=20),
    independent=dict(replicates=500,trees=10,maps_per_tree=20),
    mcmc=dict(replicates=500,trees=60,maps_per_tree=20,persistence=.8),
    routes=dict(replicates=500,map_grid=[50,200,1000,5000],tie_margin=.25,alpha=.05),
    adaptive=dict(replicates=500,min_maps=200,batch_size=200,max_maps=5000,target_rmcse=.05,min_mean_count=1.0))


def exact_moments(model,params):
    """First/second factorial moments via marked CTMC and tree product rules."""
    ref=ReferenceTree(model,params,root_weights=model.root_sampling_probs)
    states=ref.states;n=len(states)
    dmat=model._layer_arrays(0,params)[0]
    labels=["d","e","y","s","v"]+[f"route:{a}>{z}" for a in range(model.na) for z in range(model.na) if a!=z]
    result={}
    for label in labels:
        reward=np.zeros_like(ref.q);marked=np.zeros_like(ref.clado)
        for i,source in enumerate(states):
            for j,target in enumerate(states):
                if label=="d" and len(target)==len(source)+1 or label=="e" and len(target)==len(source)-1:
                    reward[i,j]=ref.q[i,j]
                if label.startswith("route:"):
                    a,z=map(int,label.split(":")[1].split(">"))
                    if a in source and set(target)-set(source)=={z} and len(target)==len(source)+1:
                        reward[i,j]=dmat[a,z]
        if label in ("y","s","v"):
            for i,parent in enumerate(states):
                for j,left in enumerate(states):
                    for k,right in enumerate(states):
                        kind="y" if left==parent and right==parent else "j" if (set(left)|set(right))-set(parent) else "s" if left==parent or right==parent else "v"
                        if kind==label:
                            marked[i,j,k]=ref.clado[i,j,k]
        branch={}
        block=np.block([[ref.q,reward,np.zeros_like(reward)],
                        [np.zeros_like(reward),ref.q,reward],
                        [np.zeros_like(reward),np.zeros_like(reward),ref.q]])
        for child,p in ref.transition.items():
            value=expm(block*model.tree.lengths[child])
            branch[child]=(p,value[:n,n:2*n],2*value[:n,2*n:])
        moments=np.zeros((ref.nnode,3,n));moments[:ref.ntips,0]=model.tips
        convolve=lambda c,x,y:np.einsum("ijk,j,k->i",c,x,y)
        for parent,(left,right) in zip(ref.parents,ref.steps):
            propagated=[]
            for child in (int(left),int(right)):
                p,h1,h2=branch[child];v0,v1,v2=moments[child]
                propagated.append((p@v0,p@v1+h1@v0,p@v2+2*h1@v1+h2@v0))
            l,r=propagated;c=ref.clado
            moments[parent,0]=convolve(c,l[0],r[0])
            moments[parent,1]=convolve(c,l[1],r[0])+convolve(c,l[0],r[1])+convolve(marked,l[0],r[0])
            moments[parent,2]=(convolve(c,l[2],r[0])+2*convolve(c,l[1],r[1])+convolve(c,l[0],r[2])+
                2*convolve(marked,l[1],r[0])+2*convolve(marked,l[0],r[1]))
        mass,first,second=moments[ref.root]@ref.root_weights
        mean=float(first/mass);variance=max(0.,float(second/mass+mean-mean*mean))
        result[label]=dict(mean=mean,variance=variance)
        if label in ("d","e"):
            per_branch={}
            for child,(p,h1,_h2) in branch.items():
                ref.transition[child]=h1
                per_branch[str(child+1)]=float(ref.root_weights @ ref.downpass(model.tips)[ref.root])/mass
                ref.transition[child]=p
            result[label]["branch_means"]=per_branch
    return result


def raw_counts(raw,n):
    ana=raw[1].reshape((-1,10));clado=raw[0].reshape((-1,12))
    result={}
    for name,code,table,column in (("d",0,ana,7),("e",1,ana,7),("y",4,clado,9),("s",5,clado,9),("v",6,clado,9)):
        result[name]=np.bincount(table[table[:,column]==code,0].astype(int)-1,minlength=n)
    for a in range(3):
        for z in range(3):
            if a!=z:
                rows=ana[(ana[:,7]==0)&(ana[:,8]==a)&(ana[:,9]==z)]
                result[f"route:{a}>{z}"]=np.bincount(rows[:,0].astype(int)-1,minlength=n)
    return result


def cp(success,n):
    return [0. if success==0 else float(beta.ppf(.025,success,n-success+1)),
            1. if success==n else float(beta.ppf(.975,success+1,n-success))]


def symmetric_fit():
    tree="((T1:1,T2:1):1,(T3:1,T4:1):1);"
    template=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=.3),e=dict(type="fixed",init=.2))
    return b.BioGeoBEARS(tree,{t:"11" for t in ("T1","T2","T3","T4")},template,likelihood_conditioning="tip_survival").fit()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",required=True)
    parser.add_argument("--mode",choices=("plan","smoke","formal"),default="plan")
    parser.add_argument("--sections",default=",".join(FORMAL))
    parser.add_argument("--seed",type=int,default=20261115)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--budget-json",type=Path)
    args=parser.parse_args();destination=Path(args.output_dir);destination.mkdir(parents=True,exist_ok=True)
    sections=args.sections.split(",")
    if len(set(sections))!=len(sections) or set(sections)-FORMAL.keys():
        parser.error("sections must be distinct known protocol sections")
    settings=json.loads(json.dumps(FORMAL))
    budget_record=None
    if args.mode=="formal":
        if args.budget_json is None:parser.error("formal execution requires a locked --budget-json")
        from scripts.v12_formal_protocol import load_locked_budget,select_budget_job
        budget,budget_record=load_locked_budget(args.budget_json,source_script=__file__,execute=True)
        job=select_budget_job(budget,analysis="v4")
        if job["seed"]!=args.seed:raise ValueError("V4 seed differs from the locked budget")
        requested=job["settings"]
        if set(requested)!=set(FORMAL):raise ValueError("V4 budget must retain all protocol sections")
        for section,original in FORMAL.items():
            if set(requested[section])!=set(original):raise ValueError("V4 section fields differ from the protocol")
            for key,value in original.items():
                actual=requested[section][key]
                if key=="replicates":
                    if type(actual) is not int or actual<2:raise ValueError("V4 requires at least two planned replicates")
                elif actual!=value:raise ValueError("V4 budget may reduce repetitions, not alter the statistical target or map protocol")
        settings=requested
    elif args.budget_json is not None:
        parser.error("formal budgets cannot relabel plan/smoke executions")
    if args.mode=="smoke":
        for value in settings.values():value["replicates"]=4
        settings["fixed"]["maps"]=50
        settings["routes"]["map_grid"]=[50,200]
        settings["adaptive"].update(min_maps=20,batch_size=20,max_maps=100,target_rmcse=.2)
    fit_list=[make_fit(scale) for scale in (.6,.8,1,1.2,1.4)];sym=symmetric_fit()
    moments=[exact_moments(f.model,f.params) for f in fit_list]
    source_files=[Path(__file__),ROOT/"rustbears/uncertainty/diagnostics.py",ROOT/"rustbears/uncertainty/rng.py",ROOT/"scripts/validate_v12_gate1.py"]
    hashes={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    manifest=dict(protocol="V4-full-v1",mode=args.mode,seed=args.seed,formal_settings=FORMAL,full_proposal_settings=FORMAL,
                  prospectively_amended=bool(args.mode=="formal" and settings!=FORMAL),locked_budget=budget_record,
                  execution_settings=settings,sections=sections,source_sha256=hashes,
                  models=[model_manifest(f) for f in fit_list],symmetric_model=model_manifest(sym),exact_moments=moments,
                  early_stop_definition="stop_reason=target_precision_met while any predeclared true-active event/route (exact mean>=min_mean_count) has sqrt(exact_variance/Mstop)/exact_mean > target_rmcse; budget-limit stops are reported separately.",
                  statistical_interpretation="Normal intervals are pointwise approximations; report empirical coverage and binomial uncertainty without changing seeds or asserting anytime validity.")
    manifest_path=destination/"manifest.json"
    payload=dict(task_fingerprint=task_fingerprint(manifest),manifest=canonical_manifest(manifest))
    if manifest_path.exists() and json.loads(manifest_path.read_text())!=payload:
        raise ValueError("output directory contains a different immutable protocol manifest; use a new directory")
    manifest_path.write_text(json.dumps(payload,indent=2)+"\n")
    if args.mode=="plan":
        print(json.dumps(dict(status="planned_not_executed",manifest=str(destination/"manifest.json"))))
        return
    store=CheckpointStore(destination/"records.jsonl",manifest,resume=args.resume)
    jobs=[];tasks=[]
    for section in sections:
        for rep in range(1,settings[section]["replicates"]+1):
            jobs.append((len(jobs)+1,section,rep))
            nt=settings[section].get("trees",1)
            tasks.extend(TaskIdentity("v4_full",section,str(tree),outer_rep=rep,role="draw") for tree in range(nt))
            tasks.append(TaskIdentity("v4_full",section,"selection",outer_rep=rep,role="tree_chain"))
    seeds=derive_task_seeds(args.seed,tasks)
    def getseed(section,rep,tree=0):return seeds[TaskIdentity("v4_full",section,str(tree),outer_rep=rep,role="draw")]
    samplers=[PreparedSampler(sampler_buffers(f.model,f.params),f.model.root_sampling_probs,0) for f in fit_list+[sym]]
    started=time.monotonic()
    try:
        for rid,section,rep in jobs:
            if rid in store.completed_ids:continue
            cfg=settings[section];record=dict(id=rid,section=section,replicate=rep);arrays={}
            if section=="fixed":
                raw=samplers[2].draw(1,cfg["maps"],getseed(section,rep),1,bsm=True)
                values=raw_counts(raw,cfg["maps"])
                record["estimates"]={kind:dict(mean=float(values[kind].mean()),se=float(values[kind].std(ddof=1)/np.sqrt(cfg["maps"]))) for kind in "de"}
                arrays.update(values)
                record["branch_estimates"]={};ana=raw[1].reshape((-1,10))
                for kind,code in (("d",0),("e",1)):
                    for node in moments[2][kind]["branch_means"]:
                        selected=ana[(ana[:,7]==code)&(ana[:,1]==int(node)-1)]
                        values=np.bincount(selected[:,0].astype(int)-1,minlength=cfg["maps"])
                        key=f"{kind}_node{node}"
                        record["branch_estimates"][key]=dict(mean=float(values.mean()),se=float(values.std(ddof=1)/np.sqrt(cfg["maps"])))
                        arrays[key]=values
            elif section in ("fixed_trees","independent","mcmc"):
                k,m=cfg["trees"],cfg["maps_per_tree"]
                rng=np.random.default_rng(seeds[TaskIdentity("v4_full",section,"selection",outer_rep=rep,role="tree_chain")])
                index=np.arange(5) if section=="fixed_trees" else rng.integers(0,5,size=k)
                if section=="mcmc":
                    for i in range(1,k):
                        if rng.random()<cfg["persistence"]:index[i]=index[i-1]
                record["tree_pool_indices"]=index.tolist()
                values=[raw_counts(samplers[int(t)].draw(1,m,getseed(section,rep,i),1,bsm=True),m) for i,t in enumerate(index)]
                metadata=pd.DataFrame(dict(map=np.arange(1,k*m+1),tree_id=np.repeat(np.arange(1,k+1),m),model_id=1,stratum_weight=1/k))
                mapped=SimpleNamespace(num_maps=k*m,map_metadata=metadata,tree_order=tuple(range(1,k+1)))
                design=_Design(mapped,"fixed_sample" if section=="fixed_trees" else "posterior_population",
                               None if section=="fixed_trees" else "independent" if section=="independent" else "mcmc_chain")
                record["estimates"]={}
                for kind in "de":
                    x=np.concatenate([v[kind] for v in values]);summary=design.summarize(x,.05)
                    record["estimates"][kind]=dict(mean=summary["mean"],se=summary["se_target"],
                        effective_n_trees=summary["effective_n_trees"],var_between_trunc=summary["var_between_trunc"])
                    arrays[kind]=x
            elif section=="routes":
                maximum=max(cfg["map_grid"]);values=raw_counts(samplers[-1].draw(1,maximum,getseed(section,rep),1,bsm=True),maximum)
                difference=values["route:0>1"]-values["route:1>0"]
                record["map_grid"]={}
                for n in cfg["map_grid"]:
                    mean=float(difference[:n].mean());se=float(difference[:n].std(ddof=1)/np.sqrt(n))
                    record["map_grid"][str(n)]=dict(mean=mean,se=se,status=rank_status(mean-1.96*se,mean+1.96*se,cfg["tie_margin"]))
                arrays["paired_route_difference"]=difference
            else:
                mapped=run_bsm_adaptive(fit_list[2],seed=getseed(section,rep),**{key:cfg[key] for key in
                    ("min_maps","batch_size","max_maps","target_rmcse","min_mean_count")})
                table=mapped.diagnostics.event_summary_df.set_index("event_type")
                record.update(stop_reason=mapped.stop_reason,maps=mapped.num_maps,
                    estimates={kind:dict(mean=float(table.loc[kind,"mean"]),se=float(table.loc[kind,"se_target"])) for kind in "de"})
                active={name:value for name,value in moments[2].items() if value["mean"]>=cfg["min_mean_count"]}
                unmet=[name for name,value in active.items() if np.sqrt(value["variance"]/mapped.num_maps)/value["mean"]>cfg["target_rmcse"]]
                record["true_precision_unmet_at_stop"]=unmet
                record["early_stopped_quantities"]=unmet if mapped.stop_reason=="target_precision_met" else []
                for name,frame in (("clado",mapped.clado_events_df),("ana",mapped.ana_events_df)):
                    metadata,values=frame_arrays(frame,name);record[name]=metadata;arrays.update(values)
            store.append(record,arrays)
    finally:
        for sampler in samplers:sampler.close()
    records=store.records;summaries={}
    for section in sections:
        selected=[r for r in records if r["section"]==section]
        if section=="routes":
            summaries[section]={}
            for n in settings[section]["map_grid"]:
                statuses=Counter(r["map_grid"][str(n)]["status"] for r in selected)
                summaries[section][str(n)]={status:dict(count=statuses[status],fraction=statuses[status]/len(selected),cp95=cp(statuses[status],len(selected))) for status in ("resolved","near_tie","unresolved")}
        else:
            summaries[section]={kind:report([(r["estimates"][kind]["mean"],r["estimates"][kind]["se"]) for r in selected],
                moments[2][kind]["mean"] if section in ("fixed","adaptive") else float(np.mean([v[kind]["mean"] for v in moments]))) for kind in "de"}
            if section=="fixed":
                summaries[section]["per_branch"]={}
                for kind in "de":
                    for node,truth in moments[2][kind]["branch_means"].items():
                        key=f"{kind}_node{node}"
                        summaries[section]["per_branch"][key]=report([(r["branch_estimates"][key]["mean"],r["branch_estimates"][key]["se"]) for r in selected],truth)
            if section in ("fixed_trees","independent","mcmc"):
                for kind in "de":
                    summaries[section][kind]["between_variance_quantiles"]=np.quantile([r["estimates"][kind]["var_between_trunc"] for r in selected],[.025,.5,.975]).tolist()
                    if section=="mcmc":summaries[section][kind]["ess_quantiles"]=np.quantile([r["estimates"][kind]["effective_n_trees"] for r in selected],[.025,.5,.975]).tolist()
            if section=="adaptive":
                n=len(selected);early=sum(bool(r["early_stopped_quantities"]) for r in selected);caps=sum(r["stop_reason"]=="max_maps_reached" for r in selected)
                summaries[section].update(stopping_map_quantiles=np.quantile([r["maps"] for r in selected],[0,.25,.5,.75,1]).tolist(),
                    early_stopping=dict(count=early,fraction=early/n,cp95=cp(early,n)),max_maps_reached=dict(count=caps,fraction=caps/n,cp95=cp(caps,n)))
    stable=all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest for path,digest in hashes.items())
    result=dict(status="completed" if stable else "completed_source_changed",formal=args.mode=="formal" and set(sections)==set(FORMAL),
                settings=settings,sections=sections,task_fingerprint=store.task_fingerprint,source_unchanged=stable,
                seconds=time.monotonic()-started,summaries=summaries,
                interpretation="Each rate has binomial uncertainty; neither sequence stopping nor selected routes receive an anytime/simultaneous guarantee.")
    (destination/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(status=result["status"],formal=result["formal"],records=len(records),seconds=result["seconds"])))


if __name__=="__main__":main()
