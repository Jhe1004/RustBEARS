#!/usr/bin/env python3
"""Execute one locked V3 study and retain every planned dataset's outcome.

Parallelism is over datasets; each external chain and geographic fit uses one
thread. Failed inference/diagnostics/fits are recorded and not rerun until they
pass. Posterior inference uses only sequences; true trees and geographic states
are confined to simulation/evaluation. Run with the validation Python that has
ArviZ installed. Default is a manifest-only operation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.v12_formal_protocol import load_locked_budget,select_budget_job
from scripts.validate_v12_v3_external_posterior import sha,_load_dataset


def write_json(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def verify_prepared_study(output,job,budget_sha):
    plan=json.loads((output/'study_plan.json').read_text())
    expected={key:job[key] for key in ('replicates','tips','areas','sites')}
    if plan.get('mode')!='formal' or (plan.get('locked_budget') or {}).get('sha256')!=budget_sha or plan.get('requested_dimensions')!=expected:
        raise ValueError('prepared V3 study is not this locked formal design; pilot promotion is forbidden')
    if [row['dataset'] for row in plan['datasets']]!=list(range(1,job['replicates']+1)):
        raise ValueError('prepared dataset inventory differs from the locked design')
    for index in range(1,job['replicates']+1):
        _folder,simulation=_load_dataset(output/f'dataset_{index:03d}')
        if (simulation.get('mode')!='formal' or (simulation.get('locked_budget') or {}).get('sha256')!=budget_sha
            or simulation['dataset']!=index or any(simulation[key]!=job[key] for key in ('tips','areas','clock_rate'))):
            raise ValueError('a prepared dataset is not bound to this locked formal design')


def verify_completed_dataset(folder,record):
    """Verify the retained input → external inference → import → score chain."""
    _load_dataset(folder)
    score_path=folder/'v3_scores.json'
    if sha(score_path)!=record['score_sha256']:raise ValueError('completed V3 score was modified')
    score=json.loads(score_path.read_text());import_path=folder/'posterior_import.json'
    if score['posterior_import_sha256']!=sha(import_path):raise ValueError('completed V3 posterior import changed')
    imported=json.loads(import_path.read_text());run_path=folder/'external_run.json'
    if imported['source_run_sha256']!=sha(run_path):raise ValueError('completed V3 external-run provenance changed')
    run=json.loads(run_path.read_text())
    if run['status']!='completed' or run['sequence_sha256']!=sha(folder/'alignment.nex'):
        raise ValueError('external inference is not complete for this alignment')
    for name,digest in imported['files'].items():
        if sha(folder/name)!=digest:raise ValueError('imported V3 tree bytes changed')
    for row in run['outputs'].values():
        if sha(row['path'])!=row['sha256']:raise ValueError('external V3 output bytes changed')
    diag=json.loads(Path(run['outputs']['diagnostics_file']['path']).read_text())
    if not diag.get('diagnostic_pass') or diag['min_ess']<200 or diag['max_rhat']>1.01:
        raise ValueError('completed V3 dataset lacks passing chain diagnostics')
    for name,digest in diag.get('raw_files',{}).items():
        if sha(folder/name)!=digest:raise ValueError('raw molecular chain provenance changed')


def summarize(records,output,seed):
    completed=[r for r in records if r['status']=='complete']
    cells={};rng=np.random.default_rng(seed)
    for support in ('all','<0.7','0.7-0.95','>0.95'):
        datasets=[];calibration={name:[dict(lower=i/10,upper=(i+1)/10,n=0,sum_prediction=0.,sum_observed=0.) for i in range(10)] for name in ('mcc','multitree')}
        for record in completed:
            score=json.loads((output/f'dataset_{record["dataset"]:03d}'/'v3_scores.json').read_text())
            selected=[r for r in score['paired_clades'] if support=='all' or r['support_bin']==support]
            if not selected:continue
            row=dict(dataset=record['dataset'],n_clades=len(selected))
            for name in ('mcc','multitree'):
                zero=any(r[name]['zero_probability'] for r in selected)
                row[name]=dict(brier=float(np.mean([r[name]['brier'] for r in selected])),
                    mean_negative_log_probability=None if zero else float(np.mean([r[name]['log_score'] for r in selected])),
                    infinite_log_score=zero)
                for r in selected:
                    for state,p in enumerate(r[name]['probabilities']):
                        bin_index=min(9,int(p*10));cell=calibration[name][bin_index]
                        cell['n']+=1;cell['sum_prediction']+=p;cell['sum_observed']+=int(state==r['truth_state'])
            row['paired_brier_difference']=row['multitree']['brier']-row['mcc']['brier']
            datasets.append(row)
        difference=np.array([r['paired_brier_difference'] for r in datasets]);ci=None
        if len(difference)>=2:
            boot=difference[rng.integers(0,len(difference),size=(5000,len(difference)))].mean(axis=1)
            ci=np.quantile(boot,[.025,.975]).tolist()
        cells[support]=dict(n_datasets_with_scored_clades=len(datasets),n_clades=sum(r['n_clades'] for r in datasets),
            dataset_mean_brier={name:None if not datasets else float(np.mean([r[name]['brier'] for r in datasets])) for name in ('mcc','multitree')},
            paired_dataset_mean_brier_difference=None if not len(difference) else float(difference.mean()),
            paired_dataset_bootstrap_ci95=ci,bootstrap_draws=5000,calibration_bins=calibration,dataset_scores=datasets)
        cells[support]['negative_log_probability']={name:dict(
            mean=None if not datasets or any(r[name]['infinite_log_score'] for r in datasets) else float(np.mean([r[name]['mean_negative_log_probability'] for r in datasets])),
            infinite=any(r[name]['infinite_log_score'] for r in datasets),
            n_infinite_datasets=sum(r[name]['infinite_log_score'] for r in datasets)) for name in ('mcc','multitree')}
        cells[support]['brier_missing_result_bounds']=None
        cells[support]['paired_brier_missing_result_bounds']=None
        if support=='all' and records:
            missing=len(records)-len(datasets);n=len(records)
            cells[support]['brier_missing_result_bounds']={name:[sum(r[name]['brier'] for r in datasets)/n,
                (sum(r[name]['brier'] for r in datasets)+2*missing)/n] for name in ('mcc','multitree')}
            cells[support]['paired_brier_missing_result_bounds']=[(float(difference.sum())-2*missing)/n,(float(difference.sum())+2*missing)/n]
    return dict(n_planned=len(records),n_complete=len(completed),n_missing=len(records)-len(completed),
        missing_datasets=[dict(dataset=r['dataset'],status=r['status'],stage=r.get('stage'),error=r.get('error')) for r in records if r['status']!='complete'],
        support_bins=cells,
        interpretation='Paired Brier differences average shared-clade scores within each dataset, then average datasets. Bootstrap units are whole datasets. Negative difference favors multiple-tree summaries. Calibration counts pool scored node-state pairs and are descriptive, not independent Bernoulli trials.',
        missing_result_boundary='Scores condition on completed datasets and shared true/MCC clades. Failed datasets are retained in the planned denominator; no finite upper bound exists for missing negative log probabilities. Do not infer unconditional accuracy from complete-case scores.',
        missing_bounds_interpretation='Brier bounds use its [0,2] range and the planned dataset denominator; they are deterministic missing-result bounds, not confidence intervals. Only the all-clade target has these bounds: support-bin membership for failed datasets is unknown.',
        convergence_boundary='Passing observed numeric/clade diagnostics is required for import but does not prove exploration of unvisited topology modes.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--budget-json',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--executable',type=Path,required=True);p.add_argument('--workers',type=int,default=1)
    p.add_argument('--execute',action='store_true');p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    if args.workers<1:p.error('workers must be positive')
    budget,metadata=load_locked_budget(args.budget_json,source_script=__file__,execute=args.execute)
    job=select_budget_job(budget,analysis='v3')
    for field in ('replicates','tips','areas','sites','posterior_trees','generations','warmup','sample_every'):
        if type(job[field]) is not int or job[field]<1:raise ValueError('invalid V3 field '+field)
    if job['chains']!=2 or job['failure_policy']!='retain_without_rerun':
        raise ValueError('V3 requires two chains and retained failed datasets without outcome-driven reruns')
    if job['min_ess']!=200 or job['max_rhat']!=1.01:
        raise ValueError('V3 diagnostics are fixed at ESS >=200 and rank R-hat <=1.01')
    for name in ('scripts/validate_v12_v3_external_posterior.py','scripts/run_v12_revbayes_posterior.py','scripts/v3_jc69_coalescent.Rev.template'):
        if budget['source_sha256'].get(name)!=sha(ROOT/name):raise ValueError('V3 source differs from locked budget: '+name)
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    executable=args.executable.resolve()
    if not executable.is_file():raise FileNotFoundError(executable)
    manifest=dict(schema='rustbears.v3_study.v1',budget=metadata,settings=job,workers=args.workers,
        executable=str(executable),executable_sha256=sha(executable),python=sys.executable,
        failure_policy='Retain every planned dataset, including failed convergence; no repeated sampling until success.')
    from rustbears._ffi import _resolve_native_library_path
    import rustbears
    package_root=Path(rustbears.__file__).resolve().parent
    production_hashes={path.relative_to(package_root).as_posix():sha(path) for path in sorted(package_root.rglob('*.py'))}
    native_path=Path(_resolve_native_library_path())
    manifest.update(production_source_sha256=production_hashes,native_sha256=sha(native_path))
    path=output/'execution_manifest.json'
    if path.exists() and json.loads(path.read_text())!=manifest:raise ValueError('V3 output already has a different immutable manifest')
    if (output/'study_plan.json').exists() and not args.resume:raise ValueError('use --resume for an existing prepared study')
    if (output/'study_plan.json').exists():verify_prepared_study(output,job,metadata['sha256'])
    write_json(path,manifest)
    if not args.execute:return
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    cache=output/'runtime_cache';cache.mkdir(exist_ok=True);env.update(MPLCONFIGDIR=str(cache),TMPDIR=str(cache))
    def run(command,folder,name):
        started=time.time()
        with (folder/f'{name}.stdout.log').open('w') as out,(folder/f'{name}.stderr.log').open('w') as err:
            result=subprocess.run(command,cwd=ROOT,env=env,stdout=out,stderr=err,check=False)
        record=dict(argv=command,returncode=result.returncode,wall_seconds=time.time()-started,
                    stdout_sha256=sha(folder/f'{name}.stdout.log'),stderr_sha256=sha(folder/f'{name}.stderr.log'))
        if result.returncode:raise RuntimeError(json.dumps(record))
        return record
    if not (output/'study_plan.json').exists():
        run([sys.executable,str(ROOT/'scripts/validate_v12_v3_external_posterior.py'),'prepare','--output-dir',str(output),
            '--mode','formal','--budget-json',str(args.budget_json.resolve()),
            *[part for key in ('replicates','tips','areas','sites','clock_rate','seed') for part in ('--'+key.replace('_','-'),str(job[key]))]],output,'prepare')
    verify_prepared_study(output,job,metadata['sha256'])
    def one(index):
        folder=output/f'dataset_{index:03d}';ledger=folder/'study_outcome.json'
        if ledger.exists():
            record=json.loads(ledger.read_text())
            if record['budget_sha256']!=metadata['sha256']:raise ValueError('saved outcome uses a different budget')
            if record['status']=='complete':verify_completed_dataset(folder,record)
            return record
        record=dict(dataset=index,budget_sha256=metadata['sha256'],status='running',commands=[])
        try:
            record['stage']='molecular_inference'
            record['commands'].append(run([sys.executable,str(ROOT/'scripts/run_v12_revbayes_posterior.py'),
                '--dataset',str(folder),'--executable',str(executable),
                *[part for key in ('generations','warmup','sample_every','burnin_fraction','seed') for part in ('--'+key.replace('_','-'),str(job[key]))]],folder,'molecular_driver'))
            record['stage']='diagnostics_and_import'
            record['commands'].append(run([sys.executable,str(ROOT/'scripts/validate_v12_v3_external_posterior.py'),
                'collect','--dataset',str(folder),'--burnin-fraction','0','--n-trees',str(job['posterior_trees'])],folder,'import_driver'))
            record['stage']='biogeographic_evaluation'
            record['commands'].append(run([sys.executable,str(ROOT/'scripts/validate_v12_v3_external_posterior.py'),
                'evaluate','--dataset',str(folder),'--threads','1'],folder,'evaluate_driver'))
            record.update(status='complete',score_sha256=sha(folder/'v3_scores.json'))
        except Exception as exc:
            record.update(status='failed_retained',error=f'{type(exc).__name__}: {exc}')
        write_json(ledger,record)
        return record
    rows=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(one,i):i for i in range(1,job['replicates']+1)}
        for future in as_completed(pending):
            row=future.result();rows.append(row);print(json.dumps(dict(dataset=row['dataset'],status=row['status'],stage=row.get('stage'))),flush=True)
            write_json(output/'progress.json',dict(planned=job['replicates'],outcomes=sorted(rows,key=lambda r:r['dataset'])))
    rows.sort(key=lambda r:r['dataset'])
    verify_prepared_study(output,job,metadata['sha256'])
    for row in rows:
        if row['status']=='complete':verify_completed_dataset(output/f'dataset_{row["dataset"]:03d}',row)
    if sha(native_path)!=manifest['native_sha256'] or any(sha(package_root/name)!=digest for name,digest in production_hashes.items()):
        raise ValueError('V3 implementation changed during the study; no aggregate is released')
    if sha(args.budget_json)!=metadata['sha256'] or any(sha(ROOT/name)!=digest for name,digest in budget['source_sha256'].items()):
        raise ValueError('locked V3 budget or workflow sources changed during execution')
    report=summarize(rows,output,job['seed'])
    report.update(status='all_planned_datasets_processed',formal=True,settings=job,budget_sha256=metadata['sha256'],
        all_datasets_scored=all(row['status']=='complete' for row in rows),outcomes=rows)
    write_json(output/'results.json',report)


if __name__=='__main__':
    main()
