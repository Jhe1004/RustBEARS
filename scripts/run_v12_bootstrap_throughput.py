#!/usr/bin/env python3
"""Predeclared B20 thread-throughput pilot; never a formal calibration study."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits


def safe(x):
    if isinstance(x,np.generic):x=x.item()
    if isinstance(x,float) and not np.isfinite(x):return None
    if isinstance(x,np.ndarray):return safe(x.tolist())
    if isinstance(x,dict):return {str(k):safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [safe(v) for v in x]
    return x


def write(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(safe(value),sort_keys=True,indent=2,allow_nan=False)+'\n');tmp.replace(path)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def template(b,name):
    obj=b.model_template(name).with_parameters(d=dict(min=0,max=2,init=.015),e=dict(min=0,max=2,init=.008))
    return obj.with_parameters(j=dict(min=0,init=.03)) if '+J' in name else obj


def setup(args):
    sys.path.insert(0,str(args.package_source))
    import rustbears as b
    from rustbears._ffi import _resolve_native_library_path
    from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
    roles=['fit_DEC_tip_survival','fit_DEC+J_tip_survival','bootstrap']
    tasks={r:TaskIdentity('v2b_alt_j003',r.split('_')[1] if r.startswith('fit_') else 'DEC+J','given_tree',1,role=r) for r in roles}
    derived=derive_task_seeds(20261115,tasks.values())
    seeds={r:derived[t] for r,t in tasks.items()}
    native=_resolve_native_library_path()
    manifest=dict(schema='rustbears.bootstrap_throughput.v1',mode='cost_pilot_not_formal',n_boot=20,
        worker_grid=[1,2,4,8],execution_order=[1,2,4,8],repetitions_per_worker_setting=1,
        inner_native_threads=1,blas_threads=1,maxiter=500,contract='tip_survival',models=['DEC','DEC+J'],
        source_commit=args.source_commit,software=b.__version__,package_source=str(args.package_source),native_path=str(native),native_sha256=sha(native),
        source_sha256={p.relative_to(args.package_source).as_posix():sha(p) for p in sorted((args.package_source/'rustbears').rglob('*.py'))},
        script_sha256=sha(__file__),input_sha256={n:sha(args.input_dir/n) for n in ['tree.newick','geography.data']},
        initial_templates={name:{key:dict(type=p.type,init=p.init,min=p.min,max=p.max) for key,p in template(b,name).params.items()} for name in ['DEC','DEC+J']},
        seed=20261115,rng={r:dict(identity=t.to_dict(),seed=seeds[r]) for r,t in tasks.items()},
        timing_scope='bootstrap_lrt call only, including its observed safeguard and all planned20 simulated pair refits; excludes independently repeated initial observed fits and interpreter startup',
        failures='Retain all planned replicas, optimizer attempts, exceptions and completed outputs. Never replace a seed or drop a slow/failed fit.',
        parity='Same task identity/seed per replica across worker settings; compare full inference records and p/CP values, plus report maximum score/parameter differences.',
        limitations='One fixed30tip16state dataset on a shared loaded host; one run per worker setting; not a formal size/power experiment or universal speedup.',
        versions={n:importlib.metadata.version(n) for n in ['rustbears','numpy','scipy','pandas','threadpoolctl']})
    return b,seeds,manifest


def group(args,b,seeds,manifest):
    import resource
    out=args.output
    if out.exists():raise FileExistsError('Refusing to overwrite a worker run')
    out.mkdir(parents=True)
    write(out/'manifest.json',manifest)
    calls=[];lock=threading.Lock();original=b.BioGeoBEARS.fit;phase='observed_initial'
    def measured(self,*pos,**kw):
        tick=time.perf_counter()
        with lock:identifier=len(calls)+1;calls.append(None)
        row=dict(id=identifier,phase=phase,model=self.model.name,inner_threads=kw.get('num_threads',1),maxiter=kw.get('maxiter',500))
        try:
            result=original(self,*pos,**kw)
            row.update(status='converged' if result.converged else 'nonconverged',lnL=result.lnL,params=result.params,n_evals=result.n_evals)
            return result
        except Exception as error:
            row.update(status='error',error=type(error).__name__+': '+str(error));raise
        finally:
            row['elapsed_sec']=time.perf_counter()-tick
            with lock:
                calls[identifier-1]=row
                with (out/'fit_calls.jsonl').open('a') as f:f.write(json.dumps(safe(row),sort_keys=True,allow_nan=False)+'\n')
    b.BioGeoBEARS.fit=measured
    start_load=os.getloadavg();result=None;error=None;elapsed=None
    try:
        with threadpool_limits(limits=1):
            fits=[b.BioGeoBEARS(args.input_dir/'tree.newick',args.input_dir/'geography.data',template(b,name),likelihood_conditioning='tip_survival').fit(maxiter=500,num_threads=1,seed=seeds[f'fit_{name}_tip_survival']) for name in ['DEC','DEC+J']]
            write(out/'observed_fits.json',[dict(lnL=f.lnL,params=f.params,converged=f.converged) for f in fits])
            if not all(f.converged for f in fits):raise RuntimeError('Observed initial fit did not converge')
            phase='bootstrap_B20';tick=time.perf_counter()
            try:
                result=b.bootstrap_lrt(*fits,n_boot=20,num_threads=args.group,seed=seeds['bootstrap'],fit_options=dict(maxiter=500),checkpoint_file=out/'bootstrap.jsonl')
            finally:elapsed=time.perf_counter()-tick
            result.replicates_df.to_csv(out/'replicates.csv',index=False)
            write(out/'replicates.json',result.replicates_df.to_dict('records'))
            write(out/'inference_summary.json',dict(n_planned=result.n_planned,n_failed=result.n_failed,d_obs=result.d_obs,pval_mc_conservative=result.pval_mc_conservative,pval_failure_bounds=result.pval_failure_bounds,pval_mc_ci95=result.pval_mc_ci95))
            info=threadpool_info()
    except Exception as exception:
        error=type(exception).__name__+': '+str(exception);info=threadpool_info()
    finally:b.BioGeoBEARS.fit=original
    measured_calls=[c for c in calls if c['phase']=='bootstrap_B20']
    report=dict(status='failed' if error else 'completed_with_replica_failures' if result.n_failed else 'complete',error=error,workers=args.group,
        bootstrap_wall_sec=elapsed,bootstrap_fit_calls=len(measured_calls),bootstrap_fit_evaluations=sum(c.get('n_evals',0) for c in measured_calls),
        bootstrap_replicates_per_sec=20/elapsed if elapsed else None,
        all_inner_threads_one=all(c['inner_threads']==1 for c in calls),fit_failures=sum(c['status']!='converged' for c in calls),
        process_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        output_bytes_before_summary=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()),
        output_scope='Includes checkpoint,all fit-call audit,replicate CSV+JSON and manifests; no input/source/runtime copies',
        load_average_start=start_load,load_average_end=os.getloadavg(),threadpool_info=info,
        source_stable=all(sha(args.package_source/name)==value for name,value in manifest['source_sha256'].items()) and sha(manifest['native_path'])==manifest['native_sha256'])
    write(out/'result.json',report);print(json.dumps(safe(report)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True);parser.add_argument('--package-source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--source-commit',required=True)
    parser.add_argument('--execute',action='store_true');parser.add_argument('--group',type=int,choices=[1,2,4,8],help=argparse.SUPPRESS)
    args=parser.parse_args()
    for name in ['input_dir','package_source','output']:setattr(args,name,getattr(args,name).resolve())
    b,seeds,manifest=setup(args)
    if args.group:return group(args,b,seeds,manifest)
    args.output.mkdir(parents=True,exist_ok=True);path=args.output/'manifest.json'
    if path.exists() and json.loads(path.read_text())!=safe(manifest):raise ValueError('Manifest changed; use a new directory')
    write(path,manifest)
    if not args.execute:print('Manifest only; no B20 calculation started');return
    results=[]
    for worker in manifest['execution_order']:
        output=args.output/f'workers{worker}'
        command=[sys.executable,__file__,'--input-dir',str(args.input_dir),'--package-source',str(args.package_source),'--output',str(output),'--source-commit',args.source_commit,'--group',str(worker)]
        with (args.output/f'workers{worker}.log').open('w') as f:process=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=False)
        record=json.loads((output/'result.json').read_text()) if (output/'result.json').exists() else dict(status='process_failed',returncode=process.returncode,workers=worker)
        results.append(record);write(args.output/'progress.json',results);print(json.dumps(safe(record)),flush=True)
    parity=[]
    for worker in [2,4,8]:
        left=args.output/'workers1';right=args.output/f'workers{worker}'
        if not all((d/'replicates.json').exists() for d in [left,right]):parity.append(dict(workers=worker,status='unavailable'));continue
        a=json.loads((left/'replicates.json').read_text());c=json.loads((right/'replicates.json').read_text())
        sa=json.loads((left/'inference_summary.json').read_text());sc=json.loads((right/'inference_summary.json').read_text())
        score_errors=[abs(x[k]-y[k]) for x,y in zip(a,c) for k in ['d_stat','obj_null','obj_alt'] if x[k] is not None and y[k] is not None]
        parameter_errors=[abs(x[k][p]-y[k][p]) for x,y in zip(a,c) for k in ['params_null','params_alt'] for p in (x.get(k) or {}) if p in (y.get(k) or {})]
        parity.append(dict(workers=worker,full_records_equal=a==c,inference_summary_equal=sa==sc,max_score_difference=max(score_errors,default=0),max_parameter_difference=max(parameter_errors,default=0),task_seeds_equal=all(x['rng']==y['rng'] for x,y in zip(a,c))))
    baseline=results[0].get('bootstrap_wall_sec')
    for row in results:
        wall=row.get('bootstrap_wall_sec')
        row['wall_speedup_relative_to_workers1']=baseline/wall if baseline and wall else None
    write(args.output/'results.json',dict(mode='cost_pilot_not_formal',worker_results=results,parity=parity,
        inference_parity_passed=len(parity)==3 and all(p.get('full_records_equal') and p.get('inference_summary_equal') for p in parity),
        scope=manifest['limitations']))


if __name__=='__main__':main()
