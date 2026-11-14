"""Predeclared V1/V2b/V2c jobs; manifest-only unless --execute is supplied.

Smoke mode exercises paths with tiny budgets and cannot estimate coverage,
size or power. Formal jobs are separate immutable manifests, never a resume
or relabelling of smoke output. No formal job is started by importing this file.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import hypergeom

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO))
import rustbears as b
from rustbears.uncertainty.checkpoint import CheckpointStore,model_manifest,task_fingerprint
from rustbears.uncertainty.guards import clone_model
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds

FAMILIES=("DEC","DEC+J","DIVALIKE","DIVALIKE+J","BAYAREALIKE","BAYAREALIKE+J")


def balanced_tree(n,height=10.):
    def subtree(labels,age,parent=None):
        branch="" if parent is None else f":{parent-age:.17g}"
        if len(labels)==1: return labels[0]+branch
        mid=len(labels)//2
        left,right=labels[:mid],labels[mid:]
        return "("+subtree(left,age*.65 if len(left)>1 else 0,age)+","+subtree(right,age*.65 if len(right)>1 else 0,age)+")"+branch
    return subtree([f"T{i+1}" for i in range(n)],height)+";"


def finite_json(value):
    if isinstance(value,np.generic): value=value.item()
    if isinstance(value,float) and not np.isfinite(value): return value.hex()
    if isinstance(value,dict): return {str(k):finite_json(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [finite_json(v) for v in value]
    return value


def template(family):
    obj=b.model_template(family).with_parameters(d=dict(min=0,max=2,init=.015),e=dict(min=0,max=2,init=.008))
    if "+J" in family: obj=obj.with_parameters(j=dict(min=0,init=.03))
    return obj


def score_nodes(probabilities,truth):
    onehot=np.eye(probabilities.shape[1])[truth]
    assigned=probabilities[np.arange(len(truth)),truth]
    with np.errstate(divide="ignore"): logs=np.log(assigned)
    calibration=[]
    flatp=probabilities.ravel();flaty=onehot.ravel()
    for low,high in zip(np.linspace(0,1,11)[:-1],np.linspace(0,1,11)[1:]):
        mask=(flatp>=low)&((flatp<high) if high<1 else (flatp<=high))
        calibration.append(dict(lower=float(low),upper=float(high),n=int(mask.sum()),
            sum_prediction=float(flatp[mask].sum()),sum_observed=float(flaty[mask].sum())))
    return dict(brier=float(np.mean(np.sum((probabilities-onehot)**2,axis=1))),
                log_score=float(np.mean(logs)),n_zero_true_probabilities=int((assigned==0).sum()),calibration=calibration)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis",choices=("v1","v2b","v2c"),required=True)
    parser.add_argument("--mode",choices=("smoke","formal"),required=True)
    parser.add_argument("--scenario",default="interior")
    parser.add_argument("--n-tips",type=int)
    parser.add_argument("--areas",type=int)
    parser.add_argument("--tree",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--seed",type=int,default=20261115)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--finite-b-extrapolation",action="store_true")
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    smoke=args.mode=="smoke"
    tips=args.n_tips or (8 if smoke else 40 if args.analysis=="v1" else 30)
    areas=args.areas or (2 if smoke else 4)
    if tips<6 or areas<2: parser.error("at least 6 tips and 2 areas are required")
    if args.analysis=="v1" and not smoke and (tips not in (40,150) or areas not in (4,5)):
        parser.error("formal V1 requires n-tips 40/150 and areas 4/5")
    scenarios={"interior":dict(d=.015,e=.008,j=.05),"e_zero":dict(d=.015,e=0,j=.05),
               "j_zero":dict(d=.015,e=.008,j=0),"double_zero":dict(d=.015,e=0,j=0),
               "null_e0":dict(d=.015,e=0,j=0),"null_e02":dict(d=.015,e=.02,j=0),
               "alt_j001":dict(d=.015,e=.02,j=.01),"alt_j003":dict(d=.015,e=.02,j=.03),
               "alt_j008":dict(d=.015,e=.02,j=.08),"dominant":dict(d=.015,e=.02,j=.08),
               "ambiguous":dict(d=.015,e=.02,j=.01)}
    allowed={"v1":("interior","e_zero","j_zero","double_zero"),
             "v2b":("null_e0","null_e02","alt_j001","alt_j003","alt_j008"),
             "v2c":("dominant","ambiguous")}
    if args.scenario not in allowed[args.analysis]: parser.error("scenario must be one of "+str(allowed[args.analysis]))
    tree=args.tree.read_text() if args.tree else balanced_tree(tips)
    parsed=b.read_newick(tree)
    if parsed.tip_count!=tips: parser.error("input tree tip count disagrees with --n-tips")
    geog=b.TipRanges(tuple(f"A{i+1}" for i in range(areas)),{name:(i%areas,) for i,name in enumerate(parsed.tip_labels)})
    base=b.BioGeoBEARS(tree,geog,template("DEC+J"),likelihood_conditioning="tip_survival")
    outer=2 if smoke else (300 if args.analysis=="v2c" else 200)
    inner=3 if smoke else (399 if args.finite_b_extrapolation else 199)
    plugin=3 if smoke else 199;alt_count=2 if smoke else 199
    truth=scenarios[args.scenario]
    config=dict(analysis=args.analysis,mode=args.mode,scope="workflow smoke only" if smoke else "predeclared formal job",
        scenario=args.scenario,scenario_is_a_design_intent_not_a_guaranteed_model_weight_regime=True,
        n_tips=tips,areas=areas,tree=tree,outer_replicates=outer,inner_bootstrap=inner,plugin_bootstrap=plugin,
        surrogate_alt_replicates=alt_count,truth=truth,seed=args.seed,alpha=.05,
        fit_options=dict(maxiter=50 if smoke else 500),profile_grid=4 if smoke else 16,
        precision="strict_f64",root_prior="uniform_nonnull",generating_contract="tip_survival",
        tested_contracts=["tip_survival","legacy_unconditioned","C_explicit_refit_to_A"] if args.analysis=="v2b" else ["tip_survival"],
        candidate_models=list(FAMILIES) if args.analysis=="v2c" else None,
        score_target="branch_top_at_internal_node_before_cladogenesis" if args.analysis=="v2c" else None,
        finite_b_extrapolation=bool(args.finite_b_extrapolation),
        finite_b_method="hypergeometric conditional rejection probabilities at B=99,199,399; paired outer-replicate linear extrapolation in 1/B",
        finite_b_reference="https://doi.org/10.1080/01621459.2000.10474226",
        failures="retain every planned outer replicate; no silent removal; report missing-result bounds",
        logical_nested_model_problems=outer*(inner+1)*2*3 if args.analysis=="v2b" else None,
        minimum_nested_fit_calls=outer*(inner+1)*3*3 if args.analysis=="v2b" else None,
        budget_excludes_surrogate_retries_profile_and_extra_guard_refits=True,
        resolved_input=model_manifest(base),
        source_sha256={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted((REPO/"rustbears").rglob("*.py"))},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.mkdir(parents=True,exist_ok=True)
    manifest_path=args.output/"preregistered_manifest.json"
    # This review copy contains a canonical fingerprint; CheckpointStore below
    # retains the lossless full manifest including numerical array bytes.
    from rustbears.uncertainty.checkpoint import canonical_manifest
    review=dict(task_fingerprint=task_fingerprint(config),manifest=canonical_manifest(config))
    if manifest_path.exists() and json.loads(manifest_path.read_text())!=review:
        raise ValueError("output contains a different preregistration; use a separate directory")
    manifest_path.write_text(json.dumps(review,indent=2))
    if not args.execute:
        print(f"Wrote {args.mode} preregistration; no simulation started: {manifest_path}")
        return
    store=CheckpointStore(args.output/"outer_records.jsonl",config,resume=args.resume)
    roles=["simulate","profile","A","C","C_to_A"]+[
        "fit_"+family+"_"+contract for family in FAMILIES for contract in ("tip_survival","legacy_unconditioned")]
    identities={(rep,role):TaskIdentity(args.analysis+"_"+args.scenario,
                    role.split("_")[1] if role.startswith("fit_") else "DEC+J","given_tree",rep,role=role)
                for rep in range(1,outer+1) for role in roles}
    seeds=derive_task_seeds(args.seed,identities.values())
    with b.PreparedSimulation(base,truth) as sampler:
        for rep in range(1,outer+1):
            if rep in store.completed_ids: continue
            batch=sampler.draw(1,start_rep=rep,seed=seeds[identities[(rep,"simulate")]])
            data=batch.to_tip_ranges(0);arrays={"true_node_states":batch.node_states,"tip_states":batch.tip_states}
            row=dict(id=rep,rep_id=rep,status="complete",errors=[],outputs={},
                     rng={role:dict(identity=identities[(rep,role)].to_dict(),seed=seeds[identities[(rep,role)]]) for role in roles})
            def fit_family(family,contract="tip_survival"):
                model=clone_model(base,geog=data,model_object=template(family),likelihood_conditioning=contract)
                return model.fit(**config["fit_options"],seed=np.random.default_rng(seeds[identities[(rep,"fit_"+family+"_"+contract)]]))
            try:
                if args.analysis=="v1":
                    fit=fit_family("DEC+J")
                    for method in ("asymptotic_chi2","bootstrap_plugin"):
                        try:
                            result=b.profile_ci(fit,params=["d","e","j"],method=method,n_boot_plugin=plugin,
                                include_onesided_zero_ref=True,n_grid=config["profile_grid"],num_threads=args.workers,
                                seed=seeds[identities[(rep,"profile")]],fit_options=config["fit_options"],
                                checkpoint_file=args.output/f"profile-{rep}-{method}.jsonl",resume=args.resume)
                            for name,interval in result.intervals.items():
                                low=0 if interval.lower is None else interval.lower
                                high=fit.model.model.params[name].max if interval.upper is None else interval.upper
                                se=interval.se_hessian_local
                                row["outputs"][method+"_"+name]=dict(lower=interval.lower,upper=interval.upper,
                                    status_lower=interval.status_lower,status_upper=interval.status_upper,
                                    covered=bool(low<=truth[name]<=high),mle=interval.mle,
                                    one_sided_upper_reference=interval.onesided_upper_from_zero,
                                    one_sided_reference_covers=None if interval.onesided_upper_from_zero is None else bool(truth[name]<=interval.onesided_upper_from_zero),
                                    wald_covers=None if se is None else bool(abs(interval.mle-truth[name])<=1.959963984540054*se))
                        except (ValueError,RuntimeError,FloatingPointError) as error:
                            row["errors"].append(method+": "+str(error))
                elif args.analysis=="v2b":
                    null_c,alt_c=fit_family("DEC","legacy_unconditioned"),fit_family("DEC+J","legacy_unconditioned")
                    null_a,alt_a=fit_family("DEC"),fit_family("DEC+J")
                    routes=(("A",null_a,alt_a,{}),("C",null_c,alt_c,dict(simulation_conditioning="tip_survival")),
                            ("C_to_A",null_c,alt_c,dict(target_conditioning="tip_survival",refit_observed=True)))
                    for route,nf,af,options in routes:
                        result=b.bootstrap_lrt(nf,af,n_boot=inner,calc_surrogate_power=True,n_alt_sim=alt_count,
                            num_threads=args.workers,seed=seeds[identities[(rep,route)]],fit_options=config["fit_options"],
                            checkpoint_file=args.output/f"lrt-{rep}-{route}.jsonl",resume=args.resume,**options)
                        summary=dict(p_conservative=result.pval_mc_conservative,reject=result.pval_mc_conservative<=.05,
                            p_failure_bounds=result.pval_failure_bounds,p_mc_ci95=result.pval_mc_ci95,n_failed=result.n_failed,
                            d_obs=result.d_obs,asymptotic_p_nominal=result.pval_asymptotic,
                            surrogate_power_fixed_threshold=result.surrogate_power_fixed_threshold,
                            recomputed_observed=result.recomputed_observed_under_new_contract)
                        if args.finite_b_extrapolation and not smoke:
                            r=result.n_exceedances;failed=result.n_failed
                            summary["subsample_rejection_probability_bounds"]={str(B):[
                                float(hypergeom.cdf(int(round(.05*(B+1)))-1,inner,r+failed,B)),
                                float(hypergeom.cdf(int(round(.05*(B+1)))-1,inner,r,B))] for B in (99,199,399)}
                        row["outputs"][route]=summary
                else:
                    fits=[fit_family(family) for family in FAMILIES]
                    average=b.model_average(fits,criterion="AICc")
                    best=fits[int(np.argmin([fit.aicc for fit in fits]))]
                    internal=np.arange(base.tree.tip_count,len(base.tree.labels))
                    true=batch.node_states[0,internal]
                    for name,prob in (("model_average",average.ancestral_probs_branch_top),
                                      ("best_model",best.evaluation.ancestral_probs_branch_top)):
                        arrays[name]=prob
                        row["outputs"][name]=score_nodes(prob[internal],true)
                    row["outputs"]["weights"]=average.weights_df.to_dict(orient="records")
            except (ValueError,RuntimeError,FloatingPointError) as error:
                row["errors"].append(str(error))
            if row["errors"]: row["status"]="incomplete"
            store.append(row,arrays)
            print(args.mode,args.analysis,args.scenario,rep,row["status"],flush=True)
    records=store.records
    summary=dict(mode=args.mode,analysis=args.analysis,scenario=args.scenario,
                 n_planned=outer,n_records=len(records),n_incomplete=sum(r["status"]!="complete" for r in records),
                 interpretation="SMOKE: workflow execution only, no inferential conclusion" if smoke else
                 "Finite-sample conditional simulation results; all incomplete outer replicates retained")
    if args.analysis in ("v1","v2b"):
        metrics={}
        targets=([method+"_"+name for method in ("asymptotic_chi2","bootstrap_plugin") for name in ("d","e","j")]
                 if args.analysis=="v1" else ["A","C","C_to_A"])
        field="covered" if args.analysis=="v1" else "reject"
        for target in targets:
            observed=[r["outputs"][target][field] for r in records if target in r["outputs"]]
            hit=sum(observed);missing=outer-len(observed)
            metrics[target]=dict(n_observed=len(observed),n_missing=missing,n_true=hit,
                                 missing_outcome_bounds=[hit/outer,(hit+missing)/outer])
        summary["metrics"]=metrics
        if args.analysis=="v2b" and args.finite_b_extrapolation and not smoke:
            extrapolation={};design=np.column_stack([np.ones(3),1/np.array([99,199,399])])
            intercept_weights=np.linalg.pinv(design)[0]
            for route in ("A","C","C_to_A"):
                paired=[]
                for r in records:
                    values=r["outputs"].get(route,{}).get("subsample_rejection_probability_bounds")
                    if values is not None and all(values[str(B)][0]==values[str(B)][1] for B in (99,199,399)):
                        paired.append([values[str(B)][0] for B in (99,199,399)])
                estimates=np.asarray(paired)@intercept_weights if paired else np.array([])
                extrapolation[route]=dict(n_complete=len(estimates),n_missing=outer-len(estimates),
                    status="complete" if len(estimates)==outer else "incomplete_not_a_formal_estimate",
                    intercept_unclipped=None if not len(estimates) else float(estimates.mean()),
                    paired_outer_mcse=None if len(estimates)<2 else float(estimates.std(ddof=1)/np.sqrt(len(estimates))))
            summary["finite_B_extrapolation"]=extrapolation
    (args.output/"summary.json").write_text(json.dumps(finite_json(summary),indent=2,allow_nan=False))


if __name__=="__main__": main()
