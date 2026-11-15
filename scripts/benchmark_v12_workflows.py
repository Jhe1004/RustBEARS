#!/usr/bin/env python3
"""Plan or measure complete uncertainty APIs in a separate bounded cost study.

No scientific calibration grid is changed. Default: plan only. Formal execution
requires a dedicated locked budget, Gate 2/3 and an idle completed statistical
suite. Each phase is a separately bounded process; all failed/timeout attempts
are retained. Deterministically scaled trees are computational inputs, not a
Bayesian posterior sample. Execution requires POSIX process groups and wait4.
"""
import argparse
from dataclasses import fields
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import threading
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.v12_formal_protocol import load_locked_budget,select_budget_job

PHASES=('six_model_average','profile_chi2','profile_plugin','bootstrap_lrt','run_bsm','fit_over_trees')
MODELS=('DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J')
DEFAULT=dict(states=[16,64,256],workers=[1,8],repeats=3,n_tips=19,n_boot=199,n_boot_plugin=199,
    n_grid=16,bsm_maps=5000,n_trees=20,maps_per_tree=100,maxiter=500,phase_timeout_seconds=14400,seed=20261115)
SMOKE=dict(DEFAULT,states=[4],workers=[1],repeats=1,n_tips=8,n_boot=3,n_boot_plugin=3,n_grid=4,
    bsm_maps=20,n_trees=2,maps_per_tree=12,phase_timeout_seconds=120)
THREAD_ENV={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
SOURCE_FILES=('scripts/benchmark_v12_workflows.py','scripts/benchmark_v12_scaling.py',
    'scripts/run_uncertainty_validation.py','scripts/run_uncertainty_case_study.py',
    'scripts/v12_workflow_timing.py','scripts/v12_formal_protocol.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    import numpy as np
    if isinstance(value,np.generic):value=value.item()
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if isinstance(value,np.ndarray):return clean(value.tolist())
    return value


def write(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(clean(value),sort_keys=True,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)


def now():return datetime.now(timezone.utc).isoformat()


def validate_settings(settings):
    if set(settings)!=set(DEFAULT):raise ValueError('dedicated workflow settings have missing/unknown fields')
    for key in DEFAULT:
        values=settings[key] if key in ('states','workers') else [settings[key]]
        if not values or any(type(v) is not int or v<1 for v in values):raise ValueError('positive integer required: '+key)
        if key in ('states','workers') and len(values)!=len(set(values)):raise ValueError('duplicate grid value')
    if any(n<4 or n>256 or n&(n-1) for n in settings['states']):raise ValueError('states must be powers of two in [4,256]')
    if max(settings['workers'])>40 or settings['n_tips']<5 or settings['n_grid']<4 or min(settings['bsm_maps'],settings['maps_per_tree'])<12:
        raise ValueError('unsupported workflow dimensions or fewer than two maps per model')
    return settings


def cells_for(settings):
    validate_settings(settings)
    return [dict(id=f's{state}_r{rep}_w{workers}',states=state,areas=int(math.log2(state)),replicate=rep,workers=workers)
        for state in settings['states'] for rep in range(1,settings['repeats']+1) for workers in settings['workers']]


def source_manifest():
    from rustbears._ffi import _resolve_native_library_path
    paths=[ROOT/name for name in SOURCE_FILES]+sorted((ROOT/'rustbears').rglob('*.py'))
    return dict(files={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
                native_path=str(_resolve_native_library_path()),native_sha256=sha(_resolve_native_library_path()))


def verify_sources(source):
    if any(sha(ROOT/name)!=digest for name,digest in source['files'].items()) or sha(source['native_path'])!=source['native_sha256']:
        raise ValueError('source/native changed since the prospective plan')


def require_checkout_core():
    import rustbears as b
    expected=ROOT/'rustbears'
    if not (expected/'__init__.py').is_file() or Path(b.__file__).resolve().parent!=expected.resolve():
        raise ValueError('formal workflow benchmark requires the full source checkout and its own imported Python core')


def template(name):
    import rustbears as b
    model=b.model_template(name)
    changes=dict(d=dict(min=0.,max=2.,init=.015),e=dict(min=0.,max=2.,init=.008))
    if model.params['j'].type=='free':changes['j']=dict(min=0.,init=.03)
    return model.with_parameters(**changes)


def seeds_for(settings,cell):
    from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
    # Worker count deliberately absent: paired worker settings use identical streams.
    tasks={role:TaskIdentity('workflow_cost',role,str(cell['states']),cell['replicate'],role=role)
           for role in (*MODELS,*PHASES)}
    seeds=derive_task_seeds(settings['seed'],tasks.values())
    return {role:seeds[identity] for role,identity in tasks.items()}


def prepare_inputs(folder,settings,cell):
    from scripts.benchmark_v12_scaling import make_inputs
    from scripts.run_uncertainty_validation import balanced_tree
    folder.mkdir(parents=True,exist_ok=True)
    make_inputs(folder,settings['n_tips'],cell['areas'])
    scales=[.9+.2*i/max(1,settings['n_trees']-1) for i in range(settings['n_trees'])]
    (folder/'scaled_trees.newick').write_text('\n'.join(balanced_tree(settings['n_tips'],height=10*s) for s in scales)+'\n')
    return dict(files={p.name:sha(p) for p in folder.iterdir() if p.is_file()},tree_scales=scales,
                interpretation='Deterministic branch-length-scaled copies; computational workload only, not posterior samples')


def save_fits(folder,fits):
    """Lossless array cache without pickle or reevaluation of the fitted MLE."""
    import numpy as np
    from rustbears.uncertainty.checkpoint import model_manifest,task_fingerprint
    arrays={};records=[]
    for i,fit in enumerate(fits):
        evaluation={}
        for field in fields(fit.evaluation):
            value=getattr(fit.evaluation,field.name)
            if isinstance(value,np.ndarray):arrays[f'{i}_{field.name}']=value;evaluation[field.name]=dict(array=f'{i}_{field.name}')
            elif field.name not in ('tree','states_list','ranges_list'):evaluation[field.name]=value
        record={f.name:getattr(fit,f.name) for f in fields(fit) if f.name not in ('evaluation','model')}
        records.append(dict(model=fit.model.model.name,record=record,evaluation=evaluation,
                            fingerprint=task_fingerprint(model_manifest(fit))))
    np.savez_compressed(folder/'fits.npz',**arrays)
    write(folder/'fits.json',dict(records=records,array_sha256=sha(folder/'fits.npz')))


def load_fits(folder,input_folder,expected_sha):
    import numpy as np
    import rustbears as b
    from rustbears.model import FitResult,EvaluationResult
    from rustbears.uncertainty.checkpoint import model_manifest,task_fingerprint
    if sha(folder/'fits.json')!=expected_sha:raise ValueError('initial fit cache metadata changed')
    data=json.loads((folder/'fits.json').read_text())
    if sha(folder/'fits.npz')!=data['array_sha256']:raise ValueError('initial fit arrays changed')
    fits=[]
    with np.load(folder/'fits.npz',allow_pickle=False) as arrays:
        for row in data['records']:
            model=b.BioGeoBEARS(input_folder/'tree.newick',input_folder/'geography.data',template(row['model']),likelihood_conditioning='tip_survival')
            kwargs={k:arrays[v['array']].copy() if isinstance(v,dict) and 'array' in v else v for k,v in row['evaluation'].items()}
            evaluation=EvaluationResult(**kwargs,tree=model.tree,states_list=model.states,ranges_list=model.ranges)
            fit=FitResult(**row['record'],evaluation=evaluation,model=model)
            if task_fingerprint(model_manifest(fit))!=row['fingerprint']:raise ValueError('reconstructed fit differs from its exact original cache')
            fits.append(fit)
    if [f.model.model.name for f in fits]!=list(MODELS) or not all(f.converged for f in fits):raise ValueError('six converged initial fits required')
    return fits


def worker(spec_path):
    import numpy as np
    import rustbears as b
    from threadpoolctl import threadpool_limits
    from rustbears.uncertainty.parallel import ParallelConfig,parallel_scope
    from scripts.v12_workflow_timing import WorkflowTimer
    from scripts.run_uncertainty_case_study import write_bsm_outputs
    spec=json.loads(Path(spec_path).read_text())
    if spec['mode']=='formal':require_checkout_core()
    verify_sources(spec['sources'])
    cell,settings=spec['cell'],spec['settings'];phase=spec['phase'];out=Path(spec['attempt']);inputs=Path(spec['inputs'])
    for name,digest in spec['input_sha256'].items():
        if sha(inputs/name)!=digest:raise ValueError('fixed computational input changed')
    checkpoint=Path(spec['checkpoint']) if spec.get('checkpoint') else None
    fits=None
    if phase not in ('six_model_average','fit_over_trees'):
        fits=load_fits(Path(spec['initial_cache']),inputs,spec['initial_cache_sha256'])
    mapping_average=b.model_average(fits) if phase=='run_bsm' else None
    timer=WorkflowTimer(out/'workflow_timings.jsonl',mode=spec['mode'],analysis='workflow_benchmark',
        n_tips=settings['n_tips'],areas=cell['areas'],states=cell['states'],workers=cell['workers'],
        native_threads=cell['workers'] if phase=='run_bsm' else 1,replicate=cell['replicate'],
        resume_requested=spec['resume'],contract='tip_survival',cell_id=cell['id'])
    seeds=seeds_for(settings,cell);lock=threading.Lock();counter=0;original=b.BioGeoBEARS.fit
    def tracked(model,*args,**kwargs):
        nonlocal counter
        with lock:
            counter+=1;identifier=counter
        def event(value):
            with lock:
                with (out/'fit_calls.jsonl').open('a') as stream:stream.write(json.dumps(clean(value),sort_keys=True,allow_nan=False)+'\n')
        base=dict(id=identifier,model=model.model.name,native_threads=kwargs.get('num_threads',1),maxiter=kwargs.get('maxiter',500))
        event(dict(base,event='started',utc=now()));started=time.monotonic();result=None
        try:
            result=original(model,*args,**kwargs)
            finished=dict(status='converged' if result.converged else 'nonconverged',params=result.params,lnL=result.lnL,n_evals=result.n_evals)
            return result
        except BaseException as error:
            finished=dict(status='raised',error=type(error).__name__+': '+str(error));raise
        finally:event(dict(base,event='finished',elapsed_seconds=time.monotonic()-started,**finished))
    b.BioGeoBEARS.fit=tracked
    current_fits=[]
    def call():
        nonlocal current_fits
        if phase=='six_model_average':
            def fit_one(name):
                try:return b.BioGeoBEARS(inputs/'tree.newick',inputs/'geography.data',template(name),likelihood_conditioning='tip_survival').fit(maxiter=settings['maxiter'],num_threads=1,seed=seeds[name])
                except Exception:return None  # Each exception is already retained; finish the other five tasks.
            with parallel_scope(ParallelConfig(cell['workers'])) as dispatch:current_fits=list(dispatch(fit_one,MODELS))
            if any(f is None or not f.converged for f in current_fits):raise RuntimeError('at least one initial model fit failed; six-model averaging unavailable')
            return b.model_average(current_fits,criterion='AICc')
        if phase in ('profile_chi2','profile_plugin'):
            return b.profile_ci(fits[1],params=['d','e','j'],method='asymptotic_chi2' if phase=='profile_chi2' else 'bootstrap_plugin',
                n_boot_plugin=settings['n_boot_plugin'],n_grid=settings['n_grid'],num_threads=cell['workers'],seed=seeds[phase],
                fit_options=dict(maxiter=settings['maxiter']),checkpoint_file=checkpoint,resume=spec['resume'])
        if phase=='bootstrap_lrt':
            return b.bootstrap_lrt(fits[0],fits[1],n_boot=settings['n_boot'],n_alt_sim=0,calc_surrogate_power=False,num_threads=cell['workers'],seed=seeds[phase],
                fit_options=dict(maxiter=settings['maxiter']),checkpoint_file=checkpoint,resume=spec['resume'])
        if phase=='run_bsm':return mapping_average.run_bsm(num_maps=settings['bsm_maps'],num_threads=cell['workers'],seed=seeds[phase])
        return b.fit_over_trees(inputs/'scaled_trees.newick',inputs/'geography.data',models=[template(n) for n in MODELS],
            summary_tree=inputs/'tree.newick',likelihood_conditioning='tip_survival',max_trees=settings['n_trees'],
            bsm_maps_per_tree=settings['maps_per_tree'],on_excess_failure='conservative_bounds',num_threads=cell['workers'],
            seed=seeds[phase],fit_options=dict(maxiter=settings['maxiter']),checkpoint_file=checkpoint,resume=spec['resume'])
    result=None;error=None
    try:
        with threadpool_limits(limits=1):
            write(out/'api_started.json',dict(monotonic=time.monotonic(),utc=now(),pid=os.getpid()))
            result=timer.run(phase,call,checkpoint=checkpoint)
        if phase=='six_model_average':
            save_fits(out,current_fits);result.weights_df.to_csv(out/'models.csv',index=False)
        elif phase.startswith('profile_'):
            result.to_dataframe().to_csv(out/'intervals.csv',index=False)
            for name,interval in result.intervals.items():
                interval.profile_curve_df.to_csv(out/f'curve_{name}.csv',index=False)
                if interval.bootstrap_plugin_df is not None:interval.bootstrap_plugin_df.to_csv(out/f'plugin_{name}.csv',index=False)
            write(out/'profile_result.json',dict(rebased_mle=result.rebased_mle,interval_kind=result.interval_kind,
                updated_fit=None if result.updated_fit is None else dict(params=result.updated_fit.params,lnL=result.updated_fit.lnL)))
        elif phase=='bootstrap_lrt':
            result.replicates_df.to_csv(out/'replicates.csv',index=False)
            write(out/'bootstrap_result.json',dict(n_planned=result.n_planned,n_failed=result.n_failed,pval=result.pval_mc_conservative,
                failure_bounds=result.pval_failure_bounds,mc_ci95=result.pval_mc_ci95,d_obs=result.d_obs))
        elif phase=='run_bsm':write_bsm_outputs(result,out)
        else:
            result.to_dataframe().to_csv(out/'probability_bounds.csv',index=False)
            result.per_tree_df.to_csv(out/'per_tree_models.csv',index=False)
            result.param_tree_sensitivity_df.to_csv(out/'parameter_sensitivity.csv',index=False)
            result.model_winner_freq.to_csv(out/'model_winner_frequency.csv')
            write(out/'multitree_result.json',dict(n_planned=result.n_trees_planned,n_complete=result.n_trees_completed,n_incomplete=result.n_trees_incomplete))
            if result.joint_bsm is None:
                write(out/'no_complete_tree_maps.json',dict(status='no_complete_tree_maps',n_complete=0,n_incomplete=result.n_trees_incomplete))
            else:write_bsm_outputs(result.joint_bsm,out,'joint_')
    except Exception as caught:
        error=type(caught).__name__+': '+str(caught)
        import traceback
        traceback.print_exc()
    finally:
        b.BioGeoBEARS.fit=original
        stable=True
        try:verify_sources(spec['sources'])
        except Exception as changed:
            stable=False;error=error or type(changed).__name__+': '+str(changed)
        write(out/'worker_result.json',dict(status='returned' if error is None else 'failed_retained',error=error,
            source_unchanged=stable,fit_calls_started=counter,phase=phase,formal=spec['mode']=='formal'))
    if error is not None:raise SystemExit(1)


def measured(command,folder,timeout,stop,lock_fds=()):
    """Bound only our process group; report process RSS and censored API time."""
    if not hasattr(os,'wait4'):raise RuntimeError('benchmark execution requires POSIX wait4/process groups')
    folder=Path(folder);cache=folder/'runtime_cache';cache.mkdir(exist_ok=True)
    env=dict(os.environ,**THREAD_ENV,MPLCONFIGDIR=str(cache),TMPDIR=str(cache),PYTHONDONTWRITEBYTECODE='1')
    started=time.monotonic();timed_out=False;interrupted=False;terminated=None
    with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stdout,stderr=stderr,start_new_session=True,pass_fds=lock_fds)
        while True:
            pid,status,usage=os.wait4(process.pid,os.WNOHANG)
            if pid:break
            elapsed=time.monotonic()-started
            if terminated is None and (elapsed>=timeout or stop.is_set()):
                timed_out=elapsed>=timeout;interrupted=stop.is_set();terminated=time.monotonic()
                try:os.killpg(process.pid,signal.SIGTERM)
                except ProcessLookupError:pass
            elif terminated is not None and time.monotonic()-terminated>2:
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
            time.sleep(.025)
        process.returncode=os.waitstatus_to_exitcode(status)
    ended=time.monotonic();api=folder/'api_started.json';lower=None;terminal_api=None
    timing=folder/'workflow_timings.jsonl'
    if timing.exists():
        # A process timeout may occur after a perfectly completed API call,
        # during output export. Such a timeout censors the process, not the API.
        for line in timing.read_bytes().splitlines(keepends=True):
            if line.endswith(b'\n'):
                candidate=json.loads(line)
                if candidate.get('status') in ('returned','raised'):terminal_api=candidate
    if (timed_out or interrupted) and api.exists():
        tick=json.loads(api.read_text())['monotonic']
        if terminal_api is None and terminated is not None and started<=tick<=terminated:lower=max(0,terminated-tick)
    return dict(command=command,returncode=process.returncode,timed_out=timed_out,interrupted=interrupted,
        process_wall_seconds=ended-started,api_right_censored_elapsed_seconds=lower,
        terminal_api_status=None if terminal_api is None else terminal_api['status'],
        terminal_api_wall_seconds=None if terminal_api is None else terminal_api['wall_seconds'],
        timeout_seconds=timeout,process_cpu_seconds=usage.ru_utime+usage.ru_stime,
        process_peak_rss_bytes=int(usage.ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        memory_scope='whole phase process high-water RSS including imports and cache load; not API-exclusive memory',
        stdout_sha256=sha(folder/'stdout.log'),stderr_sha256=sha(folder/'stderr.log'))


def completed_statistical_suite(folder,expected_budget=None,expected_sha=None):
    folder=Path(folder);plan=json.loads((folder/'suite_plan.json').read_text());ledger=json.loads((folder/'suite_ledger.json').read_text())
    ids=[j['id'] for j in plan['jobs']];rows=ledger['outcomes']
    if len(ids)!=len(set(ids)):raise ValueError('duplicate predecessor suite IDs')
    if expected_budget is not None:
        jobs=expected_budget['jobs']
        keys=[(j['analysis'],j.get('scenario'),j.get('n_tips'),j.get('areas'),bool(j.get('finite_b_extrapolation',False))) for j in jobs]
        if (plan.get('budget_sha256')!=expected_sha or len(jobs)!=28 or len(keys)!=len(set(keys))
            or {j['analysis'] for j in jobs}!={'v1','v2b','v2c','v3','v4','d1','d2'}
            or sorted(json.dumps(j['settings'],sort_keys=True) for j in plan['jobs'])!=sorted(json.dumps(j,sort_keys=True) for j in jobs)):
            raise ValueError('predecessor is not the locked 28-task statistical budget')
        if plan.get('script_sha256')!=expected_budget['source_sha256'].get('scripts/run_v12_formal_suite.py'):
            raise ValueError('predecessor launcher source mismatch')
    if (ledger.get('plan_sha256')!=sha(folder/'suite_plan.json') or ledger.get('status') not in ('all_workflows_processed','completed_with_retained_failures')
        or ledger.get('source_unchanged') is not True or len(rows)!=len(ids) or {r['id'] for r in rows}!=set(ids)
        or any(r['status'] not in ('processed','failed_retained') for r in rows)
        or ledger.get('cpu_reservations',{}).get('active_reserved_cpus')!=0):
        raise ValueError('statistical suite is not fully terminal and idle')
    return dict(plan_sha256=sha(folder/'suite_plan.json'),ledger_sha256=sha(folder/'suite_ledger.json'),status=ledger['status'])


def predecessor_budget(budget,budget_path):
    entry=budget.get('predecessor_budget',{});path=Path(entry.get('path',''))
    if not path.is_absolute():path=Path(budget_path).parent/path
    expected=budget.get('predecessor_statistical_budget_sha256')
    if not expected or entry.get('sha256')!=expected or not path.is_file() or sha(path)!=expected:
        raise ValueError('dedicated budget must bind the original statistical budget bytes')
    data=json.loads(path.read_text())
    if (data.get('schema')!='rustbears.formal_budget.v1' or data.get('status')!='locked'
        or data.get('gate_evidence',{}).get('sha256')!=budget['gate_evidence']['sha256']):
        raise ValueError('performance and statistical budgets must use the same passed gates')
    return data,expected


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path);parser.add_argument('--mode',choices=['formal','smoke'],default='formal')
    parser.add_argument('--budget-json',type=Path);parser.add_argument('--after-suite',type=Path)
    parser.add_argument('--execute',action='store_true');parser.add_argument('--resume',action='store_true')
    parser.add_argument('--worker-spec',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker_spec:return worker(args.worker_spec)
    if args.output is None:parser.error('--output required')
    settings=dict(DEFAULT if args.mode=='formal' else SMOKE);budget_record=None
    if args.budget_json:
        if args.mode!='formal':raise ValueError('smoke cannot borrow a formal budget')
        budget,budget_record=load_locked_budget(args.budget_json,source_script=__file__,execute=args.execute)
        selected=select_budget_job(budget,analysis='workflow_benchmark');settings={k:v for k,v in selected.items() if k!='analysis'}
    elif args.execute and args.mode=='formal':raise ValueError('formal workflow costs require their own locked budget')
    if args.execute and args.mode=='formal':require_checkout_core()
    validate_settings(settings);sources=source_manifest()
    if budget_record and any(budget_record['payload']['source_sha256'].get(k)!=v for k,v in sources['files'].items()):
        raise ValueError('dedicated budget must freeze every workflow source')
    if budget_record:
        native=Path(sources['native_path'])
        native_key=native.relative_to(ROOT).as_posix() if native.is_relative_to(ROOT) else None
        frozen=budget_record['payload'].get('native_sha256') or budget_record['payload']['source_sha256'].get(native_key)
        if frozen!=sources['native_sha256']:raise ValueError('dedicated budget must freeze the actual native library SHA')
    output=args.output.resolve()
    if args.execute and args.mode=='formal' and args.after_suite is None:raise ValueError('formal costs may only run after the statistical suite')
    if args.execute and (not hasattr(os,'getpriority') or os.getpriority(os.PRIO_PROCESS,0)!=0):raise ValueError('normal process priority nice=0 required')
    if args.execute and max(settings['workers'])>(os.cpu_count() or 1):raise ValueError('reserved workers exceed available logical CPUs')
    output.mkdir(parents=True,exist_ok=True)
    locks=[]
    if args.execute or os.name=='posix':
        import fcntl
        paths=[output/'benchmark.lock']+([args.after_suite/'suite.lock'] if args.execute and args.mode=='formal' else [])
        for path in paths:
            lock=path.open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(lock)
    after=None
    if args.execute and args.mode=='formal':
        predecessor,predecessor_sha=predecessor_budget(budget_record['payload'],args.budget_json)
        after=completed_statistical_suite(args.after_suite,predecessor,predecessor_sha)
    cells=cells_for(settings)
    plan=dict(schema='rustbears.workflow_benchmark.v1',mode=args.mode,formal=args.mode=='formal',settings=settings,
        host=platform.node(),platform=platform.platform(),
        full_proposal_grid=settings==DEFAULT,cells=cells,phases=list(PHASES),sources=sources,budget=budget_record,
        contract='tip_survival',root='uniform_nonnull',precision='strict_f64',parameter_lower_bounds='physical_zero',
        templates={name:{k:dict(type=p.type,init=p.init,min=p.min,max=p.max) for k,p in template(name).params.items()} for name in MODELS},
        timing='API call wall time includes internal checkpoint I/O; process time additionally includes startup, cache load and final exports',
        tree_interpretation='deterministic scaled copies for computational cost; not posterior samples',
        timeout_policy='per-phase process bound; timeout retained as right-censored cost, no inference or full-cost point',
        after_suite=None if args.after_suite is None else str(args.after_suite.resolve()))
    plan_path=output/'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text())!=plan:raise ValueError('immutable benchmark plan changed; do not reuse smoke/another budget')
    write(plan_path,plan)
    if not args.execute:return
    prior_path=output/'results.json';previous=json.loads(prior_path.read_text()) if prior_path.exists() else None
    if previous and (not args.resume or previous.get('plan_sha256')!=sha(plan_path)):raise ValueError('resume requires matching original plan')
    planned={(c['id'],p) for c in cells for p in PHASES};records={} if previous is None else {(r['cell_id'],r['phase']):r for r in previous['records']}
    if previous and (len(records)!=len(previous['records']) or set(records)-planned):raise ValueError('duplicate or unplanned benchmark records')
    terminal={'complete','timeout_retained','failed_retained','blocked_initial_fit_failure'}
    if any(r.get('status') not in terminal|{'running','interrupted'} for r in records.values()):raise ValueError('unknown benchmark task status')
    if previous is None and any((output/c['id']).exists() for c in cells):raise ValueError('unowned existing benchmark cell')
    stop=threading.Event();signal.signal(signal.SIGINT,lambda *_:stop.set());signal.signal(signal.SIGTERM,lambda *_:stop.set())
    first_started_utc=now() if previous is None else previous['first_started_utc']
    def save(status):write(prior_path,dict(schema=plan['schema'],formal=plan['formal'],plan_sha256=sha(plan_path),
        status=status,after_suite=after,records=list(records.values()),source_unchanged=True,first_started_utc=first_started_utc,updated_utc=now()))
    save('running')
    for cell in cells:
        if stop.is_set():break
        cell_dir=output/cell['id'];cell_dir.mkdir(exist_ok=True);inputs=cell_dir/'inputs'
        input_path=cell_dir/'input_manifest.json'
        if input_path.exists():
            input_manifest=json.loads(input_path.read_text())
            if any(sha(inputs/n)!=h for n,h in input_manifest['files'].items()):raise ValueError('benchmark input changed')
        else:input_manifest=prepare_inputs(inputs,settings,cell);write(input_path,input_manifest)
        for phase in PHASES:
            if stop.is_set():break
            key=(cell['id'],phase);old=records.get(key)
            if old and any(sha(path)!=digest for path,digest in old.get('files',{}).items()):raise ValueError('terminal or previously observed evidence changed')
            if old and old['status'] not in ('running','interrupted'):
                continue  # Failed/timeout tasks never automatically rerun.
            if old and old['status']=='running' and old['attempts'][-1]['status']=='running':
                # The inherited lock proves that no old child remains active.
                # Its parent died before recording wait4/clock observations.
                old['attempts'][-1].update(status='unobserved_process_exit',returncode=None,timed_out=None,interrupted=True,
                    process_wall_seconds=None,process_cpu_seconds=None,process_peak_rss_bytes=None,
                    api_right_censored_elapsed_seconds=None,terminal_api_status=None,terminal_api_wall_seconds=None,
                    process_observation='parent exit left no process termination measurement; API rows retained separately')
            initial=records.get((cell['id'],'six_model_average'))
            if phase not in ('six_model_average','fit_over_trees') and (not initial or initial['status']!='complete'):
                records[key]=dict(cell_id=cell['id'],phase=phase,status='blocked_initial_fit_failure',attempts=[]);save('running');continue
            folder=cell_dir/phase;folder.mkdir(exist_ok=True);attempts=[] if old is None else old['attempts']
            attempt=folder/f'attempt_{len(attempts)+1:03d}'
            if attempt.exists():raise ValueError('unowned pre-existing phase attempt')
            attempt.mkdir();checkpoint=folder/'checkpoint.jsonl' if phase in ('profile_chi2','profile_plugin','bootstrap_lrt','fit_over_trees') else None
            spec=dict(mode=args.mode,cell=cell,settings=settings,phase=phase,attempt=str(attempt),inputs=str(inputs),started_utc=now(),
                input_sha256=input_manifest['files'],sources=sources,checkpoint=None if checkpoint is None else str(checkpoint),
                resume=checkpoint is not None and checkpoint.exists())
            if phase not in ('six_model_average','fit_over_trees'):
                cache=Path(initial['attempts'][-1]['directory']);spec.update(initial_cache=str(cache),initial_cache_sha256=sha(cache/'fits.json'))
            write(attempt/'spec.json',spec)
            attempts.append(dict(directory=str(attempt),status='running'))
            # Prior attempt artifacts are immutable. The shared checkpoint can
            # advance during this attempt, so its old digest must not masquerade
            # as an immutable observation after a parent interruption.
            immutable={name:digest for name,digest in (old or {}).get('files',{}).items()
                       if any(Path(name).is_relative_to(Path(a['directory'])) for a in attempts[:-1])}
            immutable[str(attempt/'spec.json')]=sha(attempt/'spec.json')
            records[key]=dict(cell_id=cell['id'],phase=phase,status='running',attempts=attempts,files=immutable);save('running')
            execution=measured([sys.executable,str(Path(__file__).resolve()),'--worker-spec',str(attempt/'spec.json')],attempt,
                settings['phase_timeout_seconds'],stop,tuple(lock.fileno() for lock in locks))
            status='interrupted' if execution['interrupted'] else 'timeout_retained' if execution['timed_out'] else 'complete' if execution['returncode']==0 else 'failed_retained'
            worker_path=attempt/'worker_result.json'
            if status=='complete' and (not worker_path.exists() or json.loads(worker_path.read_text())['status']!='returned'):status='failed_retained'
            attempts[-1].update(execution,status=status)
            files={str(p):sha(p) for p in folder.rglob('*') if p.is_file() and 'runtime_cache' not in p.parts}
            records[key].update(status=status,files=files);verify_sources(sources);save('running')
    verify_sources(sources)
    save('interrupted' if stop.is_set() else 'all_planned_phases_processed')
    if stop.is_set():raise SystemExit(1)


if __name__=='__main__':main()
