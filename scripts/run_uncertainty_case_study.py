#!/usr/bin/env python3
"""Freeze and execute an empirical uncertainty workflow without selecting outcomes.

By default only the input/analysis manifest is written. --execute runs a small
smoke workflow unless --mode formal is explicitly selected. The latter requires
recorded Gate 2/3 evidence and uses the predeclared 399 bootstrap replicates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import rustbears as b
from rustbears.uncertainty.rng import TaskIdentity, derive_task_seeds
from rustbears.uncertainty.checkpoint import model_manifest, canonical_manifest
from scripts.v12_formal_protocol import load_locked_budget, select_budget_job
from scripts.v12_workflow_timing import WorkflowTimer


def dump(path, value):
    def clean(x):
        if isinstance(x, dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x, (tuple,list)):return [clean(v) for v in x]
        if isinstance(x, np.ndarray):return clean(x.tolist())
        if isinstance(x, (float,np.floating)):return float(x) if np.isfinite(x) else None
        if isinstance(x,np.integer):return int(x)
        return x
    path.write_text(json.dumps(clean(value),ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+'\n')


def case_model_template(name):
    """Predeclare physical zeros before fitting either member of a model pair.

    Profile rebasings can then improve the MLE without making the common
    d/e parameter boxes non-nested between a model and its +J alternative.
    """
    template=b.model_template(name)
    changes={key:dict(min=0.) for key in ('d','e','j')
             if template.params[key].type=='free'}
    return template.with_parameters(**changes)


def write_bsm_outputs(maps,directory,prefix=""):
    """Export complete events, design metadata and explicitly fixed-sample diagnostics."""
    if maps is None:return dict(status="no_complete_tree_maps",precision_target_met=False)
    maps.clado_events_df.to_csv(directory/(prefix+'cladogenetic_events.csv'),index=False)
    maps.ana_events_df.to_csv(directory/(prefix+'anagenetic_events.csv'),index=False)
    maps.map_metadata.to_csv(directory/(prefix+'map_metadata.csv'),index=False)
    diag=b.bsm_diagnostics(maps,tree_target='fixed_sample')
    for suffix,frame in [('event_precision',diag.event_summary_df),('route_precision',diag.route_summary_df),
                         ('route_comparisons',diag.route_rank_comparisons_df),('event_through_time',diag.ett_df)]:
        frame.to_csv(directory/(prefix+suffix+'.csv'),index=False)
    result=dict(status='complete',maps=maps.num_maps,tree_target='fixed_sample',
        summary_scope=getattr(maps,'summary_scope','all_planned_strata'),
        n_trees_missing=getattr(maps,'n_trees_missing',0),precision_target_met=diag.precision_target_met,
        interpretation='Conditional finite tree/model collection; event-count missingness has no bounded-probability guarantee.')
    dump(directory/(prefix+'bsm_summary.json'),result)
    return result


def fit_posterior_case(trees,geog,model_names,tree,contract,max_range_size,max_trees,maps_per_tree,directory,workers,seed):
    """Use exactly the physical-zero family definitions used on the reference tree."""
    return b.fit_over_trees(trees,geog,models=[case_model_template(name) for name in model_names],
        summary_tree=tree,likelihood_conditioning=contract,max_range_size=max_range_size,
        max_trees=max_trees,bsm_maps_per_tree=maps_per_tree,
        checkpoint_file=str(directory/'multitree.jsonl'),resume=True,num_threads=workers,seed=seed,
        on_excess_failure='conservative_bounds')


def fit_record(fit):
    return dict(model=fit.model.model.name,lnL=fit.lnL,params=fit.params,
        converged=fit.converged,message=fit.message,
        parameter_bounds={key:[p.min,p.max] for key,p in fit.model.model.params.items()
                          if p.type=='free'})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--tree',type=Path,default=ROOT/'r_package/RustBEARS/inst/extdata/Psychotria_5.2_ultrametric.newick')
    parser.add_argument('--geog',type=Path,default=ROOT/'r_package/RustBEARS/inst/extdata/Psychotria_geog.data')
    parser.add_argument('--posterior-trees',type=Path)
    parser.add_argument('--max-range-size',type=int)
    parser.add_argument('--case-id',default='D1_Psychotria')
    parser.add_argument('--mode',choices=['smoke','formal'],default='smoke')
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--gate-evidence',type=Path)
    parser.add_argument('--budget-json',type=Path)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--seed',type=int,default=20261115)
    args=parser.parse_args()
    tree=b.read_newick(args.tree);geog=b.read_geog(args.geog)
    max_range_size=len(geog.areas) if args.max_range_size is None else args.max_range_size
    if not 1<=max_range_size<=len(geog.areas):
        raise ValueError('max-range-size must be between 1 and the number of areas')
    hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in
        [('tree',args.tree),('geography',args.geog)]+([('posterior_trees',args.posterior_trees)] if args.posterior_trees else [])}
    models=['DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J']
    settings=dict(n_boot=399,bsm_maps=5000,max_trees=200,bsm_maps_per_tree=500)
    budget_record=None
    if args.mode=='smoke':
        if args.budget_json is not None:raise ValueError('Formal budgets cannot relabel smoke results')
        settings=dict(n_boot=3,bsm_maps=40,max_trees=3,bsm_maps_per_tree=40)
    elif args.budget_json is not None:
        budget,budget_record=load_locked_budget(args.budget_json,source_script=__file__,execute=args.execute)
        job=select_budget_job(budget,analysis='d2' if args.posterior_trees else 'd1',case_id=args.case_id)
        if job.get('seed')!=args.seed:raise ValueError('Seed differs from the locked budget')
        for name in settings:
            value=job.get(name)
            if isinstance(value,bool) or not isinstance(value,int) or value<1:raise ValueError('Invalid budget field '+name)
            settings[name]=value
        if args.gate_evidence is not None and hashlib.sha256(args.gate_evidence.read_bytes()).hexdigest()!=budget['gate_evidence']['sha256']:
            raise ValueError('Explicit gate evidence disagrees with locked budget')
    elif args.execute:
        raise ValueError('Formal execution requires --budget-json with matching Gate 2/3 evidence')
    formal_full=settings==dict(n_boot=399,bsm_maps=5000,max_trees=200,bsm_maps_per_tree=500)
    manifest=dict(case_id=args.case_id,software=b.__version__,mode=args.mode,formal=args.mode=='formal',
        input_sha256=hashes,tip_count=tree.tip_count,areas=list(geog.areas),models=models,
        likelihood_contracts=['tip_survival','legacy_unconditioned'],root_prior='uniform_nonnull',
        historical_root_weights='legacy_flat_weights',max_range_size=max_range_size,criterion='AICc',
        parameter_domain='Common d/e and free j lower bounds are physical zero, before all fits',
        **settings,seed=args.seed,workers=args.workers,locked_budget=budget_record,
        protocol_scope='workflow_smoke' if args.mode=='smoke' else 'full_proposal_case' if formal_full else 'prospectively_amended_formal_case',
        workflow_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'rustbears').rglob('*.py'))},
        resolved_reference_model=canonical_manifest(model_manifest(b.BioGeoBEARS(tree,geog,case_model_template('DEC'),
            likelihood_conditioning='tip_survival',max_range_size=max_range_size))),
        model_parameter_bounds={name:{k:[p.min,p.max] for k,p in case_model_template(name).params.items() if p.type=='free'} for name in models},
        selection='Input availability and predefined criteria only; all model outcomes retained',
        original_posterior_required=bool(args.posterior_trees),
        excess_tree_failure_policy='conservative_bounds',
        tree_failure_interpretation='Retain all failed trees and explicit unavailable mapping if none complete',
        interpretation='Conditional analysis; historical and normalized contracts are never pooled in AIC weights')
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/'manifest.json'
    if path.exists() and json.loads(path.read_text())!=manifest:raise SystemExit('Input or workflow changed; use a new output directory.')
    dump(path,manifest)
    if not args.execute:return
    if args.mode=='formal':
        from datetime import datetime, timezone
        start_path=args.out/'execution_started.json'
        if start_path.exists():
            previous=json.loads(start_path.read_text())
            if previous['budget_sha256']!=budget_record['sha256']:raise ValueError('Execution began under a different budget')
        else:
            if (args.out/'results.json').exists():raise ValueError('Formal resume requires original execution-start evidence')
            dump(start_path,dict(budget_sha256=budget_record['sha256'],first_started_utc=datetime.now(timezone.utc).isoformat(),gate2_gate3_verified=True))
    resumed=(args.out/'results.json').exists()
    timer=WorkflowTimer(args.out/'workflow_timings.jsonl',mode=args.mode,analysis='d2' if args.posterior_trees else 'd1',
        case_id=args.case_id,n_tips=tree.tip_count,areas=len(geog.areas),
        states=len(b.BioGeoBEARS(tree,geog,case_model_template('DEC'),max_range_size=max_range_size).states),
        workers=args.workers,native_threads=1,resume_requested=resumed)
    report=dict(status='running',formal=manifest['formal'],contracts={})
    tasks=[TaskIdentity('empirical_fit',str(i+1),args.case_id,role=c)
           for c in manifest['likelihood_contracts'] for i in range(len(models))]
    seeds=derive_task_seeds(args.seed,tasks)
    dump(args.out/'results.json',report)
    for contract in manifest['likelihood_contracts']:
        directory=args.out/contract;directory.mkdir(exist_ok=True)
        fits=[]
        for i,name in enumerate(models):
            seed=seeds[TaskIdentity('empirical_fit',str(i+1),args.case_id,role=contract)]
            fit=timer.run('fit',lambda:b.BioGeoBEARS(tree,geog,case_model_template(name),likelihood_conditioning=contract,
                max_range_size=max_range_size).fit(seed=seed,num_threads=1),model=name,contract=contract,call_workers=1)
            dump(directory/f'initial_fit_{i+1}.json',fit_record(fit))
            if not fit.converged:raise RuntimeError(f'{contract}/{name}: fit did not converge; outcome retained as a failure')
            fits.append(fit)
        dump(directory/'initial_fitted_parameters.json',[fit_record(f) for f in fits])
        profiles=[]
        for i,fit in enumerate(fits):
            try:
                checkpoint=directory/f'profile_{i+1}.jsonl'
                ci=timer.run('profile_asymptotic_chi2',lambda:b.profile_ci(fit,n_grid=8 if args.mode=='smoke' else 16,num_threads=args.workers,
                    checkpoint_file=str(directory/f'profile_{i+1}.jsonl'),
                    resume=checkpoint.exists()),checkpoint=checkpoint,model=fit.model.model.name,contract=contract)
                ci.to_dataframe().to_csv(directory/f'profile_{i+1}.csv',index=False)
                for param,interval in ci.intervals.items():interval.profile_curve_df.to_csv(directory/f'profile_{i+1}_{param}_curve.csv',index=False)
                if ci.updated_fit is not None:fits[i]=ci.updated_fit
                profiles.append(dict(model=fit.model.model.name,status='complete',rebased_mle=ci.rebased_mle))
            except b.ProfileOptimizationError as error:
                profiles.append(dict(model=fit.model.model.name,status='numerically_uncertified',error=str(error)))
            dump(directory/'profile_status.json',profiles)
        dump(directory/'fitted_parameters.json',[fit_record(f) for f in fits])
        average=timer.run('model_average',lambda:b.model_average(fits,criterion='AICc'),contract=contract,
            n_models=len(fits),scope='aggregation only; initial fits and profiles are timed separately')
        average.weights_df.to_csv(directory/'models.csv',index=False)
        np.savez_compressed(directory/'ancestral_arrays.npz',average_top=average.ancestral_probs_branch_top,
            average_bottom=average.ancestral_probs_branch_bottom,
            **{f'model_{i+1}_top':f.evaluation.ancestral_probs_branch_top for i,f in enumerate(fits)})
        bootstrap=[]
        for pair in range(3):
            checkpoint=directory/f'bootstrap_{pair+1}.jsonl'
            result=timer.run('bootstrap_lrt',lambda:b.bootstrap_lrt(fits[pair*2],fits[pair*2+1],n_boot=manifest['n_boot'],
                simulation_conditioning='tip_survival' if contract=='legacy_unconditioned' else 'from_model',
                checkpoint_file=str(checkpoint),resume=True,
                num_threads=args.workers,seed=args.seed+pair),checkpoint=checkpoint,contract=contract,
                null=models[pair*2],alt=models[pair*2+1],n_boot=manifest['n_boot'])
            result.replicates_df.to_csv(directory/f'bootstrap_{pair+1}.csv',index=False)
            bootstrap.append(dict(null=result.null_model,alt=result.alt_model,d_obs=result.d_obs,
                pval=result.pval_mc_conservative,failure_bounds=result.pval_failure_bounds,
                mc_ci95=result.pval_mc_ci95,n_failed=result.n_failed))
        maps=timer.run('run_bsm',lambda:average.run_bsm(num_maps=manifest['bsm_maps'],num_threads=args.workers,seed=args.seed),
            contract=contract,planned_maps=manifest['bsm_maps'],native_threads=args.workers)
        mapping=write_bsm_outputs(maps,directory)
        report['contracts'][contract]=dict(bootstrap=bootstrap,profiles=profiles,
            precision_target_met=mapping['precision_target_met'],single_tree_mapping=mapping)
        dump(args.out/'results.json',report)
        if args.posterior_trees:
            multi=timer.run('fit_over_trees',lambda:fit_posterior_case(args.posterior_trees,geog,models,tree,contract,max_range_size,
                manifest['max_trees'],manifest['bsm_maps_per_tree'],directory,args.workers,args.seed),
                checkpoint=directory/'multitree.jsonl',contract=contract,planned_trees=manifest['max_trees'],
                maps_per_tree=manifest['bsm_maps_per_tree'],scope='all geographic fits, cold audits and joint maps')
            multi.to_dataframe().to_csv(directory/'multitree_probability_bounds.csv',index=False)
            multi.per_tree_df.to_csv(directory/'per_tree_models.csv',index=False)
            multi.param_tree_sensitivity_df.to_csv(directory/'tree_parameter_sensitivity.csv',index=False)
            multi.model_winner_freq.to_csv(directory/'tree_model_winner_frequency.csv')
            joint=write_bsm_outputs(multi.joint_bsm,directory,'joint_')
            report['contracts'][contract]['multitree']=dict(n_planned=multi.n_trees_planned,
                n_complete=multi.n_trees_completed,n_incomplete=multi.n_trees_incomplete,
                fit_method=multi.fit_method,joint_mapping=joint)
            dump(args.out/'results.json',report)
    profiles_complete=all(p['status']=='complete' for c in report['contracts'].values() for p in c['profiles'])
    trees_complete=all(c.get('multitree',{}).get('n_incomplete',0)==0 for c in report['contracts'].values())
    report['status']='complete' if profiles_complete and trees_complete else 'incomplete_profiles' if not profiles_complete else 'complete_with_tree_failures'
    dump(args.out/'results.json',report)


if __name__=='__main__':main()
