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
    if args.mode=='formal':
        if args.gate_evidence is None:raise SystemExit('Formal workflow requires actual Gate 2/3 evidence, not a CI configuration.')
        gate=json.loads(args.gate_evidence.read_text())
        if not gate.get('gate2_passed') or not gate.get('gate3_passed'):raise SystemExit('Gate 2 and Gate 3 must both be passed.')
    manifest=dict(case_id=args.case_id,software=b.__version__,mode=args.mode,formal=args.mode=='formal',
        input_sha256=hashes,tip_count=tree.tip_count,areas=list(geog.areas),models=models,
        likelihood_contracts=['tip_survival','legacy_unconditioned'],root_prior='uniform_nonnull',
        historical_root_weights='legacy_flat_weights',max_range_size=max_range_size,criterion='AICc',
        parameter_domain='Common d/e and free j lower bounds are physical zero, before all fits',
        n_boot=399 if args.mode=='formal' else 3,bsm_maps=5000 if args.mode=='formal' else 40,
        max_trees=200 if args.mode=='formal' else 3,seed=args.seed,workers=args.workers,
        selection='Input availability and predefined criteria only; all model outcomes retained',
        original_posterior_required=bool(args.posterior_trees),
        interpretation='Conditional analysis; historical and normalized contracts are never pooled in AIC weights')
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/'manifest.json'
    if path.exists() and json.loads(path.read_text())!=manifest:raise SystemExit('Input or workflow changed; use a new output directory.')
    dump(path,manifest)
    if not args.execute:return
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
            fit=b.BioGeoBEARS(tree,geog,case_model_template(name),likelihood_conditioning=contract,
                max_range_size=max_range_size).fit(seed=seed,num_threads=1)
            dump(directory/f'initial_fit_{i+1}.json',fit_record(fit))
            if not fit.converged:raise RuntimeError(f'{contract}/{name}: fit did not converge; outcome retained as a failure')
            fits.append(fit)
        dump(directory/'initial_fitted_parameters.json',[fit_record(f) for f in fits])
        profiles=[]
        for i,fit in enumerate(fits):
            try:
                ci=b.profile_ci(fit,n_grid=8 if args.mode=='smoke' else 16,num_threads=args.workers,
                    checkpoint_file=str(directory/f'profile_{i+1}.jsonl'),
                    resume=(directory/f'profile_{i+1}.jsonl').exists())
                ci.to_dataframe().to_csv(directory/f'profile_{i+1}.csv',index=False)
                for param,interval in ci.intervals.items():interval.profile_curve_df.to_csv(directory/f'profile_{i+1}_{param}_curve.csv',index=False)
                if ci.updated_fit is not None:fits[i]=ci.updated_fit
                profiles.append(dict(model=fit.model.model.name,status='complete',rebased_mle=ci.rebased_mle))
            except b.ProfileOptimizationError as error:
                profiles.append(dict(model=fit.model.model.name,status='numerically_uncertified',error=str(error)))
            dump(directory/'profile_status.json',profiles)
        dump(directory/'fitted_parameters.json',[fit_record(f) for f in fits])
        average=b.model_average(fits,criterion='AICc')
        average.weights_df.to_csv(directory/'models.csv',index=False)
        np.savez_compressed(directory/'ancestral_arrays.npz',average_top=average.ancestral_probs_branch_top,
            average_bottom=average.ancestral_probs_branch_bottom,
            **{f'model_{i+1}_top':f.evaluation.ancestral_probs_branch_top for i,f in enumerate(fits)})
        bootstrap=[]
        for pair in range(3):
            result=b.bootstrap_lrt(fits[pair*2],fits[pair*2+1],n_boot=manifest['n_boot'],
                simulation_conditioning='tip_survival' if contract=='legacy_unconditioned' else 'from_model',
                checkpoint_file=str(directory/f'bootstrap_{pair+1}.jsonl'),resume=True,
                num_threads=args.workers,seed=args.seed+pair)
            result.replicates_df.to_csv(directory/f'bootstrap_{pair+1}.csv',index=False)
            bootstrap.append(dict(null=result.null_model,alt=result.alt_model,d_obs=result.d_obs,
                pval=result.pval_mc_conservative,failure_bounds=result.pval_failure_bounds,
                mc_ci95=result.pval_mc_ci95,n_failed=result.n_failed))
        maps=average.run_bsm(num_maps=manifest['bsm_maps'],num_threads=args.workers,seed=args.seed)
        maps.clado_events_df.to_csv(directory/'cladogenetic_events.csv',index=False)
        maps.ana_events_df.to_csv(directory/'anagenetic_events.csv',index=False)
        diag=b.bsm_diagnostics(maps)
        diag.event_summary_df.to_csv(directory/'event_precision.csv',index=False)
        diag.route_rank_comparisons_df.to_csv(directory/'route_comparisons.csv',index=False)
        report['contracts'][contract]=dict(bootstrap=bootstrap,profiles=profiles,precision_target_met=diag.precision_target_met)
        dump(args.out/'results.json',report)
        if args.posterior_trees:
            multi=b.fit_over_trees(args.posterior_trees,geog,models=models,summary_tree=tree,
                likelihood_conditioning=contract,max_range_size=max_range_size,
                max_trees=manifest['max_trees'],bsm_maps_per_tree=40 if args.mode=='smoke' else 500,
                checkpoint_file=str(directory/'multitree.jsonl'),resume=True,num_threads=args.workers,seed=args.seed)
            multi.to_dataframe().to_csv(directory/'multitree_probability_bounds.csv',index=False)
            multi.per_tree_df.to_csv(directory/'per_tree_models.csv',index=False)
    report['status']='complete' if all(p['status']=='complete' for c in report['contracts'].values() for p in c['profiles']) else 'incomplete_profiles'
    dump(args.out/'results.json',report)


if __name__=='__main__':main()
