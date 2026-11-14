#!/usr/bin/env python3
"""Measure complete fit cost on the same declared inputs on Mac and Linux.

This is a hardware pilot, not a statistical calibration. Each executor group
runs in a fresh subprocess so peak RSS does not inherit earlier groups. A
single process owns all outer threads and limits native and BLAS threads to 1.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears.uncertainty.checkpoint import canonical_manifest,model_manifest,task_fingerprint
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds

CASES={16:(4,4),57:(6,4),128:(7,7)}


def write_json(path,value):
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,np.ndarray):return clean(x.tolist())
        if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
        if isinstance(x,np.integer):return int(x)
        return x
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(clean(value),indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def make_model(states,contract):
    areas,maximum=CASES[states]
    names=tuple(f'A{i+1}' for i in range(areas))
    geog=b.TipRanges(names,dict(a=(0,),b=(1,),c=(0,1),d=(2,)))
    model=b.BioGeoBEARS('((a:1,b:1):1,(c:1,d:1):1);',geog,'DEC+J',
        max_range_size=maximum,likelihood_conditioning=contract)
    if len(model.states)!=states:raise AssertionError('P0 state-space declaration drifted')
    return model


def run_group(args):
    # resource is unavailable on Windows; P0 is explicitly a Mac/Linux task.
    import resource
    model=make_model(args.states,args.contract)
    source=model_manifest(model)
    tasks=[TaskIdentity('p0_refit','DEC+J',str(args.states),outer_rep=i,role='fit')
           for i in range(args.repetitions)]
    seeds=derive_task_seeds(args.seed,tasks)
    def fit_one(index):
        local=make_model(args.states,args.contract)
        start=time.perf_counter()
        fit=local.fit(seed=seeds[tasks[index]],num_threads=1,maxiter=args.maxiter)
        return dict(replicate=index,seed=f'{seeds[tasks[index]]:016x}',
                    elapsed_sec=time.perf_counter()-start,converged=fit.converged,
                    lnL=fit.lnL,params=fit.params,n_evals=fit.n_evals,message=fit.message)
    with threadpool_limits(limits=1):
        start=time.perf_counter()
        if args.workers==1:rows=list(map(fit_one,range(args.repetitions)))
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                rows=list(executor.map(fit_one,range(args.repetitions)))
        elapsed=time.perf_counter()-start
    usage=resource.getrusage(resource.RUSAGE_SELF)
    rss=int(usage.ru_maxrss)*(1 if sys.platform=='darwin' else 1024)
    result=dict(schema='rustbears-p0-v1',scope='four-tip complete fit cost; no large-tree extrapolation',
        states=args.states,areas=CASES[args.states][0],max_range_size=CASES[args.states][1],
        contract=args.contract,workers=args.workers,native_threads=1,blas_threads=1,
        repetitions=args.repetitions,maxiter=args.maxiter,wall_sec=elapsed,
        refits_per_sec=args.repetitions/elapsed,peak_rss_bytes=rss,
        median_fit_sec=float(np.median([r['elapsed_sec'] for r in rows])),
        p95_fit_sec=float(np.quantile([r['elapsed_sec'] for r in rows],.95)),
        n_converged=sum(r['converged'] for r in rows),rows=rows,
        task_fingerprint=task_fingerprint(source),model_manifest=canonical_manifest(source),
        source_stable=task_fingerprint(source)==task_fingerprint(model_manifest(model)),
        software=b.__version__,platform=platform.platform(),cpu_count=os.cpu_count(),
        python=sys.version,load_average=list(os.getloadavg()) if hasattr(os,'getloadavg') else None,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_json(args.output,result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repetitions',type=int,default=20)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--worker-grid',default='1,2,4')
    parser.add_argument('--states',type=int,choices=CASES,default=16)
    parser.add_argument('--contract',choices=['tip_survival','legacy_unconditioned'],default='tip_survival')
    parser.add_argument('--maxiter',type=int,default=500)
    parser.add_argument('--seed',type=int,default=20261115)
    parser.add_argument('--group',action='store_true',help=argparse.SUPPRESS)
    parser.add_argument('--smoke',action='store_true',help='Two fits per group; never a completed P0')
    args=parser.parse_args()
    if min(args.repetitions,args.workers,args.maxiter)<1:raise ValueError('counts must be positive')
    if args.smoke:args.repetitions=2
    if args.group:return run_group(args)
    workers=[int(w) for w in args.worker_grid.split(',')]
    if not workers or min(workers)<1 or max(workers)>256 or len(set(workers))!=len(workers):
        raise ValueError('worker-grid must contain distinct counts from 1 to 256')
    folder=args.output;folder.mkdir(parents=True,exist_ok=True)
    if (folder/'results.json').exists():raise FileExistsError('use a new P0 output directory')
    records=[]
    for states in CASES:
        for contract,worker in [('tip_survival',w) for w in workers]+[('legacy_unconditioned',1)]:
            path=folder/f'{states}_{contract}_{worker}.json'
            command=[sys.executable,str(Path(__file__).resolve()),'--group','--output',str(path),
                '--states',str(states),'--contract',contract,'--workers',str(worker),
                '--repetitions',str(args.repetitions),'--maxiter',str(args.maxiter),'--seed',str(args.seed)]
            with path.with_suffix('.log').open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            record=json.loads(path.read_text());record['command']=command;records.append(record)
            print(json.dumps(dict(states=states,contract=contract,workers=worker,
                converged=record['n_converged'],refits_per_sec=record['refits_per_sec'])),flush=True)
            write_json(folder/'progress.json',dict(status='running',groups=records))
    parity=[]
    for states in CASES:
        selected=[r for r in records if r['states']==states and r['contract']=='tip_survival']
        base=selected[0]['rows']
        differences=[abs(a['lnL']-c['lnL']) for group in selected for a,c in zip(base,group['rows'])]
        parameter_differences=[abs(a['params'][key]-c['params'][key]) for group in selected
            for a,c in zip(base,group['rows']) for key in a['params']]
        parity.append(dict(states=states,max_lnL_difference=max(differences),
                           max_parameter_difference=max(parameter_differences)))
    write_json(folder/'results.json',dict(status='complete',smoke=args.smoke,
        p0_repetition_requirement_met=not args.smoke and 20<=args.repetitions<=50,
        all_converged=all(r['n_converged']==args.repetitions for r in records),
        scope='One hardware pilot only; formal simulation budgets require both Mac and AMD pilots.',
        groups=records,worker_parity=parity))


if __name__=='__main__':main()
