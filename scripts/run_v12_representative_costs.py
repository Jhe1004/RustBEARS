#!/usr/bin/env python3
"""Small declared cost pilot; not a coverage, size, power, or accuracy experiment.

The default writes a manifest. --execute runs each declared dataset in a fresh
serial subprocess. Every native and BLAS thread budget is one. Parameters and
seeds are fixed before results are seen; numerical failures are retained.
"""
import argparse
from dataclasses import replace
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
from rustbears.uncertainty.guards import clone_model
from rustbears.uncertainty.checkpoint import model_manifest,task_fingerprint,canonical_manifest
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds

CASES=[(n,s) for n in (30,50,150) for s in ('null_e02','alt_j003')]+[(40,'interior')]
TRUTH={'null_e02':dict(d=.015,e=.02,j=0.),'alt_j003':dict(d=.015,e=.02,j=.03),'interior':dict(d=.015,e=.008,j=.05)}
PHASE='initialization'

def safe(x):
    if isinstance(x,np.generic):x=x.item()
    if isinstance(x,float) and not np.isfinite(x):return x.hex()
    if isinstance(x,np.ndarray):return safe(x.tolist())
    if isinstance(x,dict):return {str(k):safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [safe(v) for v in x]
    return x

def write(path,x):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(safe(x),indent=2,sort_keys=True,allow_nan=False)+'\n');tmp.replace(path)

def tree(n,height=10.):
    def visit(labels,age,parent=None):
        edge='' if parent is None else ':'+format(parent-age,'.17g')
        if len(labels)==1:return labels[0]+edge
        m=len(labels)//2;l,r=labels[:m],labels[m:]
        return '('+visit(l,age*.65 if len(l)>1 else 0,age)+','+visit(r,age*.65 if len(r)>1 else 0,age)+')'+edge
    return visit([f'T{i+1}' for i in range(n)],height)+';'

def template(name):
    obj=b.model_template(name).with_parameters(d=dict(min=0,max=2,init=.015),e=dict(min=0,max=2,init=.008))
    return obj.with_parameters(j=dict(min=0,init=.03)) if '+J' in name else obj

def configuration(args):
    return dict(schema='rustbears.representative_costs.v1',mode='cost_pilot_not_formal',seed=args.seed,
        cases=[dict(tips=n,areas=4,states=16,scenario=s,truth=TRUTH[s]) for n,s in CASES],
        tree='deterministic balanced, root_height=10, each internal daughter age .65*parent',
        maxiter=500,workers=1,native_threads=1,blas_threads=1,root_prior='uniform_nonnull',
        fit_contracts=['tip_survival','legacy_unconditioned'],fit_models=['DEC','DEC+J'],
        v1_40tip_initial_fit='DEC+J / tip_survival only',profile_parameter='j',profile_grid=16,
        n_boot_plugin=3,n_boot_lrt=3,
        measured_30tip_alt=['A/C paired model fits','A asymptotic j profile','A plugin3 j profile','A B3 bootstrap_lrt'],
        measured_40tip_interior=['A DEC+J fit','A asymptotic j profile','A plugin3 j profile'],
        conditional_150tip_alt_j_profile='Run once only if initial A/DEC+J fit converges within 120 seconds; otherwise retain an explicit cost-threshold skip. This rule is fixed before execution.',
        failures='Keep every failed fit and phase, including partial point/bootstrap checkpoints; no seed search or replacement dataset.',
        unmeasured='No d/e profiles; no 50tip profile; no 150tip plugin; no formal B199/B399 or outer coverage/power repeats.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),software=b.__version__)

def run_case(args,manifest):
    import resource
    global PHASE
    n,s=args.case.split(':');n=int(n);out=args.output
    out.mkdir(parents=True,exist_ok=True)
    if (out/'case_result.json').exists():raise FileExistsError('Refusing to overwrite a prior cost result')
    text=tree(n);labels=b.read_newick(text).tip_labels
    geo=b.TipRanges(tuple('A'+str(i+1) for i in range(4)),{name:(i%4,) for i,name in enumerate(labels)})
    base=b.BioGeoBEARS(text,geo,template('DEC+J'),likelihood_conditioning='tip_survival')
    roles=['simulate','profile_asymptotic','profile_plugin','bootstrap']+[f'fit_{m}_{c}' for m in ('DEC','DEC+J') for c in ('tip_survival','legacy_unconditioned')]
    analysis='v1' if s=='interior' else 'v2b'
    tasks={r:TaskIdentity(analysis+'_'+s,r.split('_')[1] if r.startswith('fit_') else 'DEC+J','given_tree',1,role=r) for r in roles}
    seeds=derive_task_seeds(args.seed,tasks.values())
    write(out/'case_manifest.json',dict(configuration=manifest,tips=n,scenario=s,truth=TRUTH[s],tree=text,
        input_fingerprint=task_fingerprint(model_manifest(base)),resolved_input=canonical_manifest(model_manifest(base)),
        rng={r:dict(identity=t.to_dict(),seed=seeds[t]) for r,t in tasks.items()},platform=platform.platform(),
        python=sys.version,load_average=os.getloadavg(),source_commit_hint=args.source_commit))
    with b.PreparedSimulation(base,TRUTH[s]) as sampler:batch=sampler.draw(1,start_rep=1,seed=seeds[tasks['simulate']],num_threads=1)
    geo=batch.to_tip_ranges(0);np.savez_compressed(out/'simulated_input.npz',tips=batch.tip_states,nodes=batch.node_states)
    (out/'tree.newick').write_text(text+'\n')
    (out/'geography.data').write_text(f'{n} 4 (A1 A2 A3 A4)\n'+''.join(name+' '+''.join('1' if a in geo.ranges[name] else '0' for a in range(4))+'\n' for name in labels))
    calls=[];phases=[];fits={};fit_durations={};original_fit=b.BioGeoBEARS.fit
    def measured_fit(self,*pos,**kw):
        tick=time.perf_counter();row=dict(id=len(calls)+1,phase=PHASE,model=self.model.name,
            contract=self.likelihood_conditioning,free_parameters=list(self.model.free),maxiter=kw.get('maxiter'),num_threads=kw.get('num_threads'))
        try:
            result=original_fit(self,*pos,**kw)
            row.update(status='converged' if result.converged else 'nonconverged',lnL=result.lnL,n_evals=result.n_evals,params=result.params,message=result.message)
            return result
        except Exception as error:
            row.update(status='error',error=type(error).__name__+': '+str(error));raise
        finally:
            row['elapsed_sec']=time.perf_counter()-tick;calls.append(row)
            with (out/'fit_calls.jsonl').open('a') as f:f.write(json.dumps(safe(row),sort_keys=True,allow_nan=False)+'\n')
    b.BioGeoBEARS.fit=measured_fit
    def phase(name,operation):
        global PHASE
        PHASE=name;tick=time.perf_counter();first=len(calls)
        try:
            value=operation();status='complete';error=None
        except Exception as e:
            value=None;status='failed';error=type(e).__name__+': '+str(e)
        record=dict(name=name,status=status,error=error,elapsed_sec=time.perf_counter()-tick,
            fit_calls=len(calls)-first,fit_evaluations=sum(x.get('n_evals',0) for x in calls[first:]),
            process_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
        phases.append(record);write(out/'phase_progress.json',phases)
        return value,record
    with threadpool_limits(limits=1):
        specs=[('tip_survival','DEC+J')] if s=='interior' else [(c,m) for c in ('tip_survival','legacy_unconditioned') for m in ('DEC','DEC+J')]
        for c,m in specs:
            name='initial_'+c+'_'+m
            def fit(c=c,m=m):
                return clone_model(base,geog=geo,model_object=template(m),likelihood_conditioning=c).fit(maxiter=500,num_threads=1,seed=seeds[tasks[f'fit_{m}_{c}']])
            value,record=phase(name,fit);fits[(c,m)]=value;fit_durations[(c,m)]=record['elapsed_sec']
        af=fits.get(('tip_survival','DEC+J'))
        regular=(n,s) in ((30,'alt_j003'),(40,'interior'))
        conditional=(n,s)==(150,'alt_j003')
        eligible=af is not None and af.converged
        if conditional and (not eligible or fit_durations[('tip_survival','DEC+J')]>120):
            phases.append(dict(name='asymptotic_j',status='skipped_cost_rule',initial_fit_converged=eligible,initial_fit_sec=fit_durations[('tip_survival','DEC+J')],threshold_sec=120))
        elif regular or conditional:
            for method in (['asymptotic_chi2','bootstrap_plugin'] if regular else ['asymptotic_chi2']):
                def profile(method=method):
                    if af is None or not af.converged:raise RuntimeError('Initial A/DEC+J fit failed; dependent profile unavailable')
                    value=b.profile_ci(af,params=['j'],method=method,n_grid=16,n_boot_plugin=3,num_threads=1,
                        seed=seeds[tasks['profile_plugin' if method=='bootstrap_plugin' else 'profile_asymptotic']],
                        fit_options=dict(maxiter=500),checkpoint_file=out/(method+'.jsonl'))
                    value.to_dataframe().to_csv(out/(method+'.csv'),index=False)
                    interval=value.intervals['j'];interval.profile_curve_df.to_csv(out/(method+'_curve.csv'),index=False)
                    if interval.bootstrap_plugin_df is not None:interval.bootstrap_plugin_df.to_csv(out/(method+'_plugin_replicates.csv'),index=False)
                    return value
                result,_=phase(method+'_j',profile)
                if result is not None and result.updated_fit is not None:af=result.updated_fit
        if (n,s)==(30,'alt_j003'):
            def bootstrap():
                nf=fits.get(('tip_survival','DEC'))
                if nf is None or af is None or not nf.converged or not af.converged:raise RuntimeError('Observed paired fits failed; bootstrap unavailable')
                result=b.bootstrap_lrt(nf,af,n_boot=3,num_threads=1,seed=seeds[tasks['bootstrap']],fit_options=dict(maxiter=500),checkpoint_file=out/'bootstrap.jsonl')
                result.replicates_df.to_csv(out/'bootstrap_replicates.csv',index=False)
                write(out/'bootstrap_summary.json',dict(n_planned=result.n_planned,n_failed=result.n_failed,p_failure_bounds=result.pval_failure_bounds,mc_ci95=result.pval_mc_ci95))
                return result
            phase('bootstrap_B3',bootstrap)
    b.BioGeoBEARS.fit=original_fit
    result=dict(status='completed_with_failures' if any(p['status']=='failed' for p in phases) or any(c['status']!='converged' for c in calls) else 'complete',
        mode='cost_pilot_not_formal',tips=n,scenario=s,phases=phases,fit_calls=len(calls),
        nonconverged_fit_calls=sum(c['status']=='nonconverged' for c in calls),errored_fit_calls=sum(c['status']=='error' for c in calls),
        process_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        load_average_end=os.getloadavg(),script_sha256=manifest['script_sha256'])
    write(out/'case_result.json',result);print(json.dumps(safe(result)))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=int,default=20261115);p.add_argument('--source-commit',default='unrecorded')
    p.add_argument('--execute',action='store_true');p.add_argument('--case',help=argparse.SUPPRESS)
    args=p.parse_args();args.output=args.output.resolve();manifest=configuration(args)
    if args.case:
        if tuple([int(args.case.split(':')[0]),args.case.split(':')[1]]) not in CASES:raise ValueError('Undeclared pilot case')
        return run_case(args,manifest)
    args.output.mkdir(parents=True,exist_ok=True);path=args.output/'manifest.json'
    if path.exists() and json.loads(path.read_text())!=safe(manifest):raise ValueError('Manifest changed; use a new directory')
    write(path,manifest)
    if not args.execute:print('Manifest only; no representative cost computations started');return
    rows=[]
    for n,s in CASES:
        folder=args.output/f'{n}_{s}';folder.mkdir(exist_ok=True)
        command=[sys.executable,str(Path(__file__).resolve()),'--output',str(folder),'--case',f'{n}:{s}','--seed',str(args.seed),'--source-commit',args.source_commit]
        tick=time.perf_counter()
        with (folder/'driver.log').open('w') as f:completed=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=False)
        row=dict(tips=n,scenario=s,returncode=completed.returncode,wall_sec=time.perf_counter()-tick)
        if (folder/'case_result.json').exists():row['result']=json.loads((folder/'case_result.json').read_text())
        else:row['result']=dict(status='process_failed',error='See retained driver.log; no replacement dataset.')
        rows.append(row);write(args.output/'progress.json',rows);print(json.dumps(safe(row)),flush=True)
    write(args.output/'results.json',dict(mode='cost_pilot_not_formal',status='complete' if all(r['result']['status']=='complete' for r in rows) else 'completed_with_failures',cases=rows,
        limitations=manifest['unmeasured'],source_commit_hint=args.source_commit))

if __name__=='__main__':main()
