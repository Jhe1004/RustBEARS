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
from scripts.v12_formal_protocol import load_locked_budget,select_budget_job
from scripts.v12_workflow_timing import WorkflowTimer

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


def read_locked_budget(path, analysis, scenario, tips, areas, *, execute=False, finite_b_extrapolation=False):
    """Read a prospective budget; formal execution requires its hashed release evidence."""
    budget,metadata=load_locked_budget(path,source_script=__file__,execute=execute)
    job=select_budget_job(budget,analysis=analysis,scenario=scenario,n_tips=tips,areas=areas,
                          finite_b_extrapolation=finite_b_extrapolation)
    for key in ("outer_replicates","inner_bootstrap","plugin_bootstrap","surrogate_alt_replicates"):
        value=job.get(key)
        if isinstance(value,bool) or not isinstance(value,int) or value<(0 if key=="surrogate_alt_replicates" else 1):
            raise ValueError("invalid locked budget field "+key)
    if job["outer_replicates"]<2:raise ValueError("at least two outer replicates are required; small N must retain its uncertainty")
    if isinstance(job.get("seed"),bool) or not isinstance(job.get("seed"),int) or job["seed"]<0:
        raise ValueError("locked budget requires a nonnegative integer seed")
    return job,metadata


def proportion_summary(values,planned,*,unavailable=0):
    from scipy.stats import beta
    values=[bool(v) for v in values];hit=sum(values);missing=planned-len(values)
    if missing<0:raise ValueError("more observations than planned replicates")
    return dict(n_observed=len(values),n_unavailable=int(unavailable),n_missing=missing,
        n_true=hit,observed_fraction=hit/len(values) if values else None,
        all_planned_fraction_counting_unavailable_as_noncoverage=hit/planned,
        missing_outcome_bounds=[hit/planned,(hit+missing)/planned],
        mc_cp95_envelope=[0. if hit==0 else float(beta.ppf(.025,hit,planned-hit+1)),
            1. if hit+missing==planned else float(beta.ppf(.975,hit+missing+1,planned-hit-missing))])


def distribution(values):
    values=[float(v) for v in values if v is not None]
    finite=[v for v in values if np.isfinite(v)]
    return dict(n=len(values),n_nonfinite=len(values)-len(finite),
        quantiles=None if not finite else dict(zip(("min","q025","median","q975","max"),
            np.quantile(finite,[0,.025,.5,.975,1]).tolist())))


def interval_record(interval,truth,upper_box):
    se=interval.se_hessian_local
    se=None if se is None or not np.isfinite(se) else float(se)
    wald_lower=None if se is None else interval.mle-1.959963984540054*se
    wald_upper=None if se is None else interval.mle+1.959963984540054*se
    low=0 if interval.lower is None else interval.lower
    high=upper_box if interval.upper is None else interval.upper
    return dict(lower=interval.lower,upper=interval.upper,status_lower=interval.status_lower,
        status_upper=interval.status_upper,covered=bool(low<=truth<=high),mle=interval.mle,
        interval_kind=interval.interval_kind,critical_value=interval.critical_value,
        endpoint_deviance_residual_max=interval.endpoint_deviance_residual_max,
        se_hessian_local=se,wald_lower=wald_lower,wald_upper=wald_upper,
        wald_available=se is not None,wald_unavailability_reason=None if se is not None else "boundary_or_nonpositive_local_hessian",
        wald_covers=None if se is None else bool(wald_lower<=truth<=wald_upper),
        one_sided_upper_reference=interval.onesided_upper_from_zero,
        one_sided_reference_covers=None if interval.onesided_upper_from_zero is None else bool(truth<=interval.onesided_upper_from_zero))


def summarize_v1(records,planned):
    metrics={}
    for parameter in ("d","e","j"):
        definitions=[(method,method+"_"+parameter,"covered","lower","upper")
                     for method in ("asymptotic_chi2","bootstrap_plugin")]
        definitions += [("onesided_zero_reference","asymptotic_chi2_"+parameter,"one_sided_reference_covers",None,"one_sided_upper_reference"),
                        ("wald_local","asymptotic_chi2_"+parameter,"wald_covers","wald_lower","wald_upper")]
        for method,key,field,low,high in definitions:
            rows=[r["outputs"][key] for r in records if key in r["outputs"]]
            unavailable=sum(row.get(field) is None for row in rows)
            observed=[row[field] for row in rows if row.get(field) is not None]
            summary=proportion_summary(observed,planned,unavailable=unavailable)
            summary.update(method=method,parameter=parameter,source_profile=key,
                lower_endpoint=distribution([0. if low is None and row.get(field) is not None else row.get(low) for row in rows]),
                upper_endpoint=distribution([row.get(high) for row in rows]),
                interpretation="Undefined intervals remain unavailable; no artificial Wald/one-sided interval is imputed.")
            metrics[method+"_"+parameter]=summary
    return metrics


def summarize_v2c(records,planned):
    """Use independent datasets as the replication unit, not correlated nodes."""
    paired=[r for r in records if all(name in r["outputs"] for name in ("model_average","best_model"))]
    out=dict(n_planned=planned,n_paired=len(paired),n_missing=planned-len(paired),
             replication_unit="simulated_dataset",log_score_orientation="mean log probability; larger is better; zero true probability gives -infinity",methods={})
    for name in ("model_average","best_model"):
        rows=[r["outputs"][name] for r in paired]
        brier=np.array([r["brier"] for r in rows]);logs=np.array([r["log_score"] for r in rows])
        n=len(rows);calibration=[]
        for index in range(10):
            bins=[r["calibration"][index] for r in rows]
            counts=np.array([r["n"] for r in bins],float)
            pred=np.array([r["sum_prediction"] for r in bins],float)
            obs=np.array([r["sum_observed"] for r in bins],float);total=float(counts.sum())
            frequency=float(obs.sum()/total) if total else None
            residual=obs-frequency*counts if total else np.zeros(n)
            mcse=float(np.sqrt(n/(n-1)*np.sum(residual**2))/total) if n>1 and total else None
            calibration.append(dict(lower=index/10,upper=(index+1)/10,n=int(total),
                mean_prediction=float(pred.sum()/total) if total else None,observed_frequency=frequency,
                mcse_dataset_cluster_ratio=mcse,n_dataset_clusters=n))
        total_brier=float(brier.sum())
        out["methods"][name]=dict(n_datasets=n,mean_brier=float(brier.mean()) if n else None,
            mcse_brier=float(brier.std(ddof=1)/np.sqrt(n)) if n>1 else None,
            brier_missing_result_bounds=[total_brier/planned,(total_brier+2*(planned-n))/planned],
            mean_log_probability=float(logs.mean()) if n else None,
            mcse_log_probability=float(logs.std(ddof=1)/np.sqrt(n)) if n>1 and np.isfinite(logs).all() else None,
            n_nonfinite_log_datasets=int((~np.isfinite(logs)).sum()),
            n_zero_true_probabilities=sum(r["n_zero_true_probabilities"] for r in rows),calibration=calibration)
    differences=np.array([r["outputs"]["model_average"]["brier"]-r["outputs"]["best_model"]["brier"] for r in paired])
    out["paired_brier_difference_average_minus_best"]=dict(mean=float(differences.mean()) if len(differences) else None,
        mcse=float(differences.std(ddof=1)/np.sqrt(len(differences))) if len(differences)>1 else None,
        n_datasets=len(differences),interpretation="negative favors model averaging; complete-case paired dataset contrast")
    weights=[r["outputs"].get("weights",[]) for r in paired]
    maxima=[max(float(row["weight"]) for row in table) for table in weights if table]
    out["realized_maximum_model_weight"]=distribution(maxima)
    out["scenario_interpretation"]="Prespecified higher/lower generating jump rate; realized weights describe the actual competition regime."
    return out


def summarize_v2b(records,planned):
    output={}
    for route in ("A","C","C_to_A"):
        rows=[r["outputs"][route] for r in records if route in r["outputs"]]
        available=[row for row in rows if row.get("surrogate_power_fixed_threshold") is not None]
        surrogate=np.array([row["surrogate_power_fixed_threshold"] for row in available])
        difference=np.array([row["surrogate_power_fixed_threshold"]-float(row["reject"]) for row in available])
        known_reject=sum(row["p_failure_bounds"][1]<=.05 for row in rows)
        possible_reject=sum(row["p_failure_bounds"][0]<=.05 for row in rows)+(planned-len(rows))
        output[route]=dict(bootstrap_rejection=proportion_summary([r["reject"] for r in rows],planned),
            nominal_chi2_rejection=proportion_summary([r["asymptotic_p_nominal"]<=.05 for r in rows],planned),
            rejection_bounds_including_inner_and_outer_failures=[known_reject/planned,possible_reject/planned],
            n_outer_with_inner_failures=sum(row["n_failed"]>0 for row in rows),
            surrogate=dict(n_outer=len(available),mean=float(surrogate.mean()) if len(surrogate) else None,
                mcse_across_outer_datasets=float(surrogate.std(ddof=1)/np.sqrt(len(surrogate))) if len(surrogate)>1 else None,
                missing_result_bounds=[sum(row["surrogate_power_failure_bounds"][0] for row in available)/planned,
                    (sum(row["surrogate_power_failure_bounds"][1] for row in available)+planned-len(available))/planned],
                paired_gap_mean=float(difference.mean()) if len(difference) else None,
                paired_gap_mcse=float(difference.std(ddof=1)/np.sqrt(len(difference))) if len(difference)>1 else None,
                interpretation="Descriptive gap: per-dataset fitted-alternative fixed-threshold surrogate minus repeated-sampling nested rejection. Different generating-parameter targets prevent interpreting this as isolated fixed-threshold bias."))
    return output


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
    parser.add_argument("--budget-json",type=Path,help="prospectively locked budget; required for formal execution")
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    smoke=args.mode=="smoke"
    tips=args.n_tips or (8 if smoke else 40 if args.analysis=="v1" else 30)
    areas=args.areas or (2 if smoke else 4)
    if tips<6 or areas<2: parser.error("at least 6 tips and 2 areas are required")
    if args.analysis=="v1" and not smoke and args.budget_json is None and (tips not in (40,150) or areas not in (4,5)):
        parser.error("formal V1 requires n-tips 40/150 and areas 4/5")
    scenarios={"interior":dict(d=.015,e=.008,j=.05),"e_zero":dict(d=.015,e=0,j=.05),
               "j_zero":dict(d=.015,e=.008,j=0),"double_zero":dict(d=.015,e=0,j=0),
               "null_e0":dict(d=.015,e=0,j=0),"null_e02":dict(d=.015,e=.02,j=0),
               "alt_j001":dict(d=.015,e=.02,j=.01),"alt_j003":dict(d=.015,e=.02,j=.03),
               "alt_j008":dict(d=.015,e=.02,j=.08),"higher_jump":dict(d=.015,e=.02,j=.08),
               "lower_jump":dict(d=.015,e=.02,j=.01)}
    allowed={"v1":("interior","e_zero","j_zero","double_zero"),
             "v2b":("null_e0","null_e02","alt_j001","alt_j003","alt_j008"),
             "v2c":("higher_jump","lower_jump")}
    if args.scenario not in allowed[args.analysis]: parser.error("scenario must be one of "+str(allowed[args.analysis]))
    tree=args.tree.read_text() if args.tree else balanced_tree(tips)
    parsed=b.read_newick(tree)
    if parsed.tip_count!=tips: parser.error("input tree tip count disagrees with --n-tips")
    geog=b.TipRanges(tuple(f"A{i+1}" for i in range(areas)),{name:(i%areas,) for i,name in enumerate(parsed.tip_labels)})
    base=b.BioGeoBEARS(tree,geog,template("DEC+J"),likelihood_conditioning="tip_survival")
    outer=2 if smoke else (300 if args.analysis=="v2c" else 200)
    inner=3 if smoke else (399 if args.finite_b_extrapolation else 199)
    plugin=3 if smoke else 199;alt_count=2 if smoke else 199
    budget_record=None
    if smoke and args.budget_json is not None:parser.error("locked formal budgets cannot relabel smoke output")
    if not smoke and args.execute and args.budget_json is None:
        parser.error("formal execution requires --budget-json with actual Gate 2/3 evidence")
    if args.budget_json is not None:
        job,budget_record=read_locked_budget(args.budget_json,args.analysis,args.scenario,tips,areas,execute=args.execute,finite_b_extrapolation=args.finite_b_extrapolation)
        if job["seed"]!=args.seed:parser.error("--seed differs from the prospectively locked budget")
        outer,inner,plugin,alt_count=(job[k] for k in ("outer_replicates","inner_bootstrap","plugin_bootstrap","surrogate_alt_replicates"))
        if args.finite_b_extrapolation and inner<399:parser.error("finite-B extrapolation requires at least 399 inner draws")
    proposal_counts=(300 if args.analysis=="v2c" else 200,399 if args.finite_b_extrapolation else 199,199,199)
    same_counts=(outer,inner,plugin,alt_count)==proposal_counts
    same_dimensions=(tips in (40,150) and areas in (4,5)) if args.analysis=="v1" else (tips==30 and areas==4)
    protocol_scope="workflow_smoke" if smoke else "full_proposal_cell" if same_counts and same_dimensions else "prospectively_amended_formal_cell"
    truth=scenarios[args.scenario]
    config=dict(analysis=args.analysis,mode=args.mode,scope="workflow smoke only" if smoke else "predeclared formal job",
        protocol_scope=protocol_scope,full_proposal_cell_dimensions=bool(not smoke and same_counts and same_dimensions),
        scenario=args.scenario,scenario_is_a_design_intent_not_a_guaranteed_model_weight_regime=True,locked_budget=budget_record,
        has_locked_formal_budget=bool(not smoke and budget_record),
        proposal_outer_replicates=300 if args.analysis=="v2c" else 200,
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
        actual_no_retry_nested_fit_calls_including_surrogate=outer*(6+9*(1+inner+alt_count)) if args.analysis=="v2b" else None,
        resolved_input=model_manifest(base),
        source_sha256={p.relative_to(REPO).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
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
    if not smoke:
        from datetime import datetime, timezone
        started_path=args.output/"execution_started.json"
        started=dict(budget_sha256=budget_record["sha256"],first_started_utc=datetime.now(timezone.utc).isoformat(),gate2_gate3_verified=True)
        if started_path.exists():
            existing=json.loads(started_path.read_text())
            if existing["budget_sha256"]!=budget_record["sha256"]:raise ValueError("execution began under a different locked budget")
            if datetime.fromisoformat(budget_record["locked_utc"])>datetime.fromisoformat(existing["first_started_utc"]):
                raise ValueError("budget was locked after execution began")
        else:
            if (args.output/"outer_records.jsonl").exists():raise ValueError("formal resume requires its original execution-start evidence")
            started_path.write_text(json.dumps(started,indent=2)+"\n")
    store=CheckpointStore(args.output/"outer_records.jsonl",config,resume=args.resume)
    timer=WorkflowTimer(args.output/'workflow_timings.jsonl',mode=args.mode,analysis=args.analysis,
        scenario=args.scenario,n_tips=tips,areas=areas,states=len(base.states),workers=args.workers,
        native_threads=1,resume_requested=args.resume)
    roles=["simulate","profile","A","C","C_to_A"]+[
        "fit_"+family+"_"+contract for family in FAMILIES for contract in ("tip_survival","legacy_unconditioned")]
    rng_scenario={"higher_jump":"dominant","lower_jump":"ambiguous"}.get(args.scenario,args.scenario)
    identities={(rep,role):TaskIdentity(args.analysis+"_"+rng_scenario,
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
                return timer.run('fit',lambda:model.fit(**config["fit_options"],
                    seed=np.random.default_rng(seeds[identities[(rep,"fit_"+family+"_"+contract)]])),
                    replicate=rep,model=family,contract=contract,call_workers=1)
            try:
                if args.analysis=="v1":
                    fit=fit_family("DEC+J")
                    for method in ("asymptotic_chi2","bootstrap_plugin"):
                        try:
                            checkpoint=args.output/f"profile-{rep}-{method}.jsonl"
                            result=timer.run('profile_'+method,lambda:b.profile_ci(fit,params=["d","e","j"],method=method,n_boot_plugin=plugin,
                                include_onesided_zero_ref=True,n_grid=config["profile_grid"],num_threads=args.workers,
                                seed=seeds[identities[(rep,"profile")]],fit_options=config["fit_options"],
                                checkpoint_file=checkpoint,resume=args.resume),checkpoint=checkpoint,
                                replicate=rep,params=['d','e','j'],plugin_bootstrap=plugin)
                            final_fit=result.updated_fit or fit
                            row.setdefault("profile_fits",{})[method]=dict(params=final_fit.params,lnL=final_fit.lnL,
                                rebased_mle=result.rebased_mle,interval_kind=result.interval_kind)
                            for name,interval in result.intervals.items():
                                row["outputs"][method+"_"+name]=interval_record(interval,truth[name],final_fit.model.model.params[name].max)
                            if result.updated_fit is not None:fit=result.updated_fit
                        except (ValueError,RuntimeError,FloatingPointError) as error:
                            row["errors"].append(method+": "+str(error))
                elif args.analysis=="v2b":
                    null_c,alt_c=fit_family("DEC","legacy_unconditioned"),fit_family("DEC+J","legacy_unconditioned")
                    null_a,alt_a=fit_family("DEC"),fit_family("DEC+J")
                    routes=(("A",null_a,alt_a,{}),("C",null_c,alt_c,dict(simulation_conditioning="tip_survival")),
                            ("C_to_A",null_c,alt_c,dict(target_conditioning="tip_survival",refit_observed=True)))
                    for route,nf,af,options in routes:
                        checkpoint=args.output/f"lrt-{rep}-{route}.jsonl"
                        result=timer.run('bootstrap_lrt',lambda:b.bootstrap_lrt(nf,af,n_boot=inner,calc_surrogate_power=alt_count>0,n_alt_sim=alt_count,
                            num_threads=args.workers,seed=seeds[identities[(rep,route)]],fit_options=config["fit_options"],
                            checkpoint_file=checkpoint,resume=args.resume,**options),checkpoint=checkpoint,
                            replicate=rep,route=route,n_boot=inner,n_alt_sim=alt_count)
                        summary=dict(p_conservative=result.pval_mc_conservative,reject=result.pval_mc_conservative<=.05,
                            p_failure_bounds=result.pval_failure_bounds,p_mc_ci95=result.pval_mc_ci95,n_failed=result.n_failed,
                            d_obs=result.d_obs,asymptotic_p_nominal=result.pval_asymptotic,
                            surrogate_power_fixed_threshold=result.surrogate_power_fixed_threshold,
                            surrogate_power_failure_bounds=result.surrogate_power_failure_bounds,
                            surrogate_power_mc_ci95=result.surrogate_power_mc_ci95,
                            n_alt_failed=None if result.surrogate_replicates_df is None else int((result.surrogate_replicates_df.status!="converged").sum()),
                            recomputed_observed=result.recomputed_observed_under_new_contract)
                        if args.finite_b_extrapolation and not smoke:
                            r=result.n_exceedances;failed=result.n_failed
                            summary["subsample_rejection_probability_bounds"]={str(B):[
                                float(hypergeom.cdf(int(round(.05*(B+1)))-1,inner,r+failed,B)),
                                float(hypergeom.cdf(int(round(.05*(B+1)))-1,inner,r,B))] for B in (99,199,399)}
                        row["outputs"][route]=summary
                else:
                    fits=[fit_family(family) for family in FAMILIES]
                    average=timer.run('model_average',lambda:b.model_average(fits,criterion="AICc"),
                        replicate=rep,n_models=len(fits),scope='aggregation only; six model fits have separate fit records')
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
                 protocol_scope=protocol_scope,n_planned=outer,n_records=len(records),n_incomplete=sum(r["status"]!="complete" for r in records),
                 interpretation="SMOKE: workflow execution only, no inferential conclusion" if smoke else
                 "Finite-sample conditional simulation results; all incomplete outer replicates retained")
    summary["locked_budget_sha256"]=None if budget_record is None else budget_record["sha256"]
    if args.analysis=="v1":
        summary["metrics"]=summarize_v1(records,outer)
    elif args.analysis=="v2c":
        summary["prediction"]=summarize_v2c(records,outer)
    else:
        summary["metrics"]={route:proportion_summary([r["outputs"][route]["reject"] for r in records if route in r["outputs"]],outer)
                            for route in ("A","C","C_to_A")}
        summary["comparisons"]=summarize_v2b(records,outer)
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
