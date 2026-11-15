#!/usr/bin/env python3
"""Prospective formal-suite launcher with a bounded CPU allocation and ledger.

Default: write commands only. Execution requires the hashed locked budget and
passed Gate 2/3. Statistical failures are retained by each adapter, never rerun
until they pass. This launcher controls only its own process groups.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

from v12_formal_protocol import load_locked_budget

ROOT=Path(__file__).resolve().parents[1]
THREAD_ENV=dict(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
VALIDATOR_SOURCES=('scripts/run_v12_formal_suite.py','scripts/v12_formal_protocol.py',
                   'scripts/plot_v12_calibration.py','scripts/v12_case_evidence.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    path=Path(path);temporary=path.with_name(path.name+'.partial')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    os.replace(temporary,path)


def now():return datetime.now(timezone.utc).isoformat()


def is_large_v1(job):
    return job['analysis']=='v1' and job.get('n_tips')==150 and job.get('areas')==5


def allocated_workers(job,regular,large):
    return large if is_large_v1(job) else regular


def launch_priority(job):
    analysis=job['analysis']
    ranks={'d1':0,'d2':1,'v3':2,'v4':3,'v2c':5,'v2b':6,'v1':7}
    return (4 if is_large_v1(job['settings']) else ranks[analysis],job['id'])


class CPUReservations:
    """FIFO admission using the immutable per-job CPU requests.

    A large waiting task can leave capacity idle until its full allocation is
    available. Later small tasks cannot postpone that preregistered priority.
    These are reserved worker counts, not CPU affinity or measured utilization.
    """
    def __init__(self,capacity,jobs):
        self.capacity=capacity
        self.requests={job['id']:job['reserved_cpus'] for job in jobs}
        if len(self.requests)!=len(jobs) or any(type(n) is not int or not 1<=n<=capacity for n in self.requests.values()):
            raise ValueError('invalid per-job CPU reservation')
        self.waiting=[job['id'] for job in jobs];self.active={};self.used=0;self.peak=0
        self.condition=threading.Condition()

    def acquire(self,job_id,stop):
        with self.condition:
            if job_id not in self.waiting:raise ValueError('job is not awaiting CPU admission')
            amount=self.requests[job_id]
            while True:
                if stop.is_set():
                    self.waiting.remove(job_id);self.condition.notify_all();return None
                if self.waiting[0]==job_id and self.used+amount<=self.capacity:
                    self.waiting.pop(0);self.active[job_id]=amount;self.used+=amount
                    self.peak=max(self.peak,self.used);self.condition.notify_all();return self.used
                # Signals set the shared Event; bounded waits also release jobs
                # when no running process is available to notify this condition.
                self.condition.wait(timeout=.1)

    def release(self,job_id):
        with self.condition:
            self.used-=self.active.pop(job_id);self.condition.notify_all()

    def snapshot(self):
        with self.condition:
            return dict(capacity=self.capacity,active_reserved_cpus=self.used,
                peak_reserved_cpus_this_launch=self.peak,active_allocations=dict(self.active),
                waiting_jobs=list(self.waiting))


def command_for(job,index,output,budget_path,revbayes,d2,workers):
    analysis=job['analysis'];name=f'{index:02d}_{analysis}'
    if analysis in ('v1','v2b','v2c'):
        name+=f'_n{job["n_tips"]}_a{job["areas"]}_{job["scenario"]}'
        if job.get('finite_b_extrapolation'):name+='_finiteB'
    directory=output/name
    if Path(name).name!=name or directory.parent.resolve()!=output.resolve():raise ValueError('job ID escapes its suite directory')
    prefix=[sys.executable];budget=['--budget-json',str(budget_path)]
    if analysis in ('v1','v2b','v2c'):
        argv=prefix+[str(ROOT/'scripts/run_uncertainty_validation.py'),'--analysis',analysis,'--mode','formal',
            *budget,'--scenario',job['scenario'],'--n-tips',str(job['n_tips']),'--areas',str(job['areas']),
            '--workers',str(workers),'--seed',str(job['seed']),'--output',str(directory),'--execute','--resume']
        if job.get('finite_b_extrapolation'):argv.append('--finite-b-extrapolation')
        result=directory/'summary.json'
    elif analysis=='v3':
        argv=prefix+[str(ROOT/'scripts/run_v12_v3_study.py'),*budget,'--output',str(directory),
            '--executable',str(revbayes),'--workers',str(workers),'--execute','--resume'];result=directory/'results.json'
    elif analysis=='v4':
        argv=prefix+[str(ROOT/'scripts/validate_v12_v4_protocol.py'),'--mode','formal',*budget,
            '--output-dir',str(directory),'--seed',str(job['seed']),'--resume'];result=directory/'results.json'
    elif analysis in ('d1','d2'):
        argv=prefix+[str(ROOT/'scripts/run_uncertainty_case_study.py'),'--mode','formal',*budget,
            '--workers',str(workers),'--seed',str(job['seed']),'--out',str(directory),'--execute',
            '--case-id',job['case_id']];result=directory/'results.json'
        if analysis=='d2':argv+=['--tree',str(d2/'raw/reference_biogeo.tre'),'--geog',str(d2/'raw/geography.data'),
            '--posterior-trees',str(d2/'raw/posterior.trees'),'--max-range-size','2']
    else:raise ValueError('unknown analysis '+analysis)
    return dict(id=name,analysis=analysis,settings=job,argv=argv,output=str(directory),result=str(result),
        reserved_cpus=workers,cwd=str(ROOT),environment=THREAD_ENV)


def check_result(job,budget_hash):
    from plot_v12_calibration import load_run,job_key,verify_workflow_source,read_json
    path=Path(job['result']);analysis=job['analysis'];settings=job['settings']
    argv=job['argv'];budget_path=Path(argv[argv.index('--budget-json')+1]);budget=read_json(budget_path)
    if sha(budget_path)!=budget_hash or budget['jobs'].count(settings)!=1:raise ValueError('result job is not the locked budget task')
    if analysis in ('v1','v2b','v2c','v3','v4'):
        run=load_run(dict(analysis=analysis,path=str(path.parent)),ROOT,budget_hash)
        manifest=run['manifest'];metadata=manifest['budget'] if analysis=='v3' else manifest['locked_budget']
        if metadata['payload']!=budget:raise ValueError('result stores a different budget payload')
        verify_workflow_source(manifest,analysis,budget)
        if run['key']!=job_key(settings):raise ValueError('result belongs to a different planned job')
        actual=manifest['settings'] if analysis=='v3' else manifest['execution_settings'] if analysis=='v4' else manifest
        if analysis=='v4':
            if actual!=settings['settings'] or manifest['seed']!=settings['seed']:raise ValueError('V4 settings differ from planned job')
        elif any(actual.get(k)!=v for k,v in settings.items()):raise ValueError('result settings differ from planned job')
        evidence=run['files']
    else:
        from v12_case_evidence import verify_case
        evidence=verify_case(path.parent,settings,budget,budget_hash)['files']
    return dict(path=str(path),sha256=sha(path),evidence_files=evidence,
                meaning='all planned tasks processed; retained statistical failures are not claimed successful')


def validate_ledger(prior,jobs,plan_hash):
    if prior.get('plan_sha256')!=plan_hash:raise ValueError('ledger plan mismatch')
    outcomes=prior['outcomes'];ids=[row['id'] for row in outcomes];planned={job['id'] for job in jobs}
    if len(ids)!=len(set(ids)) or set(ids)-planned:raise ValueError('duplicate or unknown ledger task IDs')
    statuses={'not_started','running','interrupted','processed','failed_retained'}
    if any(row.get('status') not in statuses for row in outcomes):raise ValueError('unknown ledger terminal/continuation status')
    if any(row['status']=='processed' and not isinstance(row.get('result'),dict) for row in outcomes):raise ValueError('processed task lacks verified result evidence')
    return {row['id']:row for row in outcomes}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--budget-json',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--revbayes',type=Path,required=True);p.add_argument('--d2-input-dir',type=Path,required=True)
    p.add_argument('--cpu-budget',type=int,default=32);p.add_argument('--workers-per-job',type=int,default=2)
    p.add_argument('--large-v1-workers',type=int,default=None,
                   help='Fixed workers for V1 with 150 tips and 5 areas; defaults to --workers-per-job')
    p.add_argument('--execute',action='store_true');p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    large_workers=args.workers_per_job if args.large_v1_workers is None else args.large_v1_workers
    if (not 1<=args.workers_per_job<=args.cpu_budget or not 1<=large_workers<=args.cpu_budget
        or args.cpu_budget>(os.cpu_count() or 1)):
        p.error('require 1 <= each worker allocation <= cpu-budget <= available logical CPUs')
    if args.execute and os.getpriority(os.PRIO_PROCESS,0)!=0:
        raise ValueError('formal suite must be launched at the authorized normal process priority (nice=0)')
    budget_path=args.budget_json.resolve();budget=json.loads(budget_path.read_text())
    if args.execute:budget,_=load_locked_budget(budget_path,source_script=__file__,execute=True)
    if not budget.get('jobs'):raise ValueError('empty budget')
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    # Protect the immutable plan as well as execution. A plan-only invocation
    # cannot overwrite a live owner's files, and two first launches cannot race.
    lock=(output/'suite.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('another suite launcher owns this output')
    jobs=[command_for(job,index,output,budget_path,args.revbayes.resolve(),args.d2_input_dir.resolve(),
                      allocated_workers(job,args.workers_per_job,large_workers))
          for index,job in enumerate(budget['jobs'],1)]
    # Pilot cost evidence fixes this priority before the formal run. Scheduling
    # changes no budget settings, random seed, or scientific task identity.
    jobs.sort(key=launch_priority)
    if len({json.dumps(j['settings'],sort_keys=True) for j in jobs})!=len(jobs):raise ValueError('duplicate budget job')
    plan=dict(schema='rustbears.formal_suite.v1',budget_sha256=sha(budget_path),budget_status=budget.get('status'),
        cpu_budget=args.cpu_budget,workers_per_job=args.workers_per_job,large_v1_workers=large_workers,
        max_concurrent_jobs=min(len(jobs),args.cpu_budget//min(j['reserved_cpus'] for j in jobs)),
        executor_threads=min(len(jobs),args.cpu_budget),admission_policy='fixed_priority_fifo_weighted_cpu_reservation',
        launch_order=['d1','d2','v3','v4','v1_150tips_5areas','v2c','v2b','remaining_v1'],
        worker_allocation_rule='150-tip/5-area V1 receives large_v1_workers; every other job receives workers_per_job; immutable across resume',
        jobs=jobs,script_sha256=sha(__file__),normal_priority=True,statistical_failure_policy='retain_without_outcome_driven_rerun')
    planpath=output/'suite_plan.json'
    if planpath.exists() and json.loads(planpath.read_text())!=plan:raise ValueError('existing suite has a different immutable plan; use a new output directory')
    write(planpath,plan)
    if not args.execute:
        print(json.dumps(dict(status='commands_written_no_execution',jobs=len(jobs))));return
    for name in VALIDATOR_SOURCES:
        if budget.get('source_sha256',{}).get(name)!=sha(ROOT/name):raise ValueError('terminal validator source is not frozen in budget: '+name)
    for name,digest in budget['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('budget source mismatch: '+name)
    if not args.revbayes.is_file():raise FileNotFoundError(args.revbayes)
    ledgerpath=output/'suite_ledger.json';existing={}
    if ledgerpath.exists():
        if not args.resume:raise ValueError('existing execution requires --resume')
        prior=json.loads(ledgerpath.read_text())
        existing=validate_ledger(prior,jobs,sha(planpath))
    if any(Path(job['output']).exists() and job['id'] not in existing for job in jobs):
        raise ValueError('unowned pre-existing job output; cannot silently adopt another run')
    if any(Path(job['output']).is_symlink() for job in jobs):raise ValueError('job output must not be a symlink to another run')
    stop=threading.Event();running={};guard=threading.Lock();rows=dict(existing)
    resources=CPUReservations(args.cpu_budget,jobs)
    def interrupt(_signal,_frame):stop.set()
    signal.signal(signal.SIGTERM,interrupt);signal.signal(signal.SIGINT,interrupt)
    def save():
        write(ledgerpath,dict(plan_sha256=sha(planpath),updated_utc=now(),pid=os.getpid(),
            status='interrupted' if stop.is_set() else 'running',cpu_reservations=resources.snapshot(),
            outcomes=sorted(rows.values(),key=lambda r:r['id'])))
    def execute_one(job,admitted_total):
        previous=existing.get(job['id'])
        if previous and previous['status'] in ('processed','failed_retained'):
            if previous['status']=='processed' and check_result(job,plan['budget_sha256'])!=previous['result']:
                raise ValueError('completed result changed since prior execution')
            return previous
        if stop.is_set():return dict(id=job['id'],status='not_started')
        cache=Path(job['output'])/'runtime_cache'
        env=dict(os.environ,**THREAD_ENV,TMPDIR=str(cache),MPLCONFIGDIR=str(cache),XDG_CACHE_HOME=str(cache))
        row=dict(id=job['id'],status='running',started_utc=now(),reserved_cpus=job['reserved_cpus'],
                 total_reserved_cpus_at_admission=admitted_total)
        with guard:rows[job['id']]=row.copy();save()
        cache.mkdir(parents=True,exist_ok=True)
        logdir=output/'logs';logdir.mkdir(exist_ok=True)
        # Append across infrastructure resume; never overwrite prior failures.
        with (logdir/(job['id']+'.stdout.log')).open('a') as stdout,(logdir/(job['id']+'.stderr.log')).open('a') as stderr:
            started=time.monotonic()
            # Children retain the flock descriptor. If the launcher itself is
            # killed abruptly, a restart cannot overlap still-running children.
            process=subprocess.Popen(job['argv'],cwd=job['cwd'],env=env,stdout=stdout,stderr=stderr,
                start_new_session=True,pass_fds=(lock.fileno(),))
            with guard:running[job['id']]=process
            term_sent=None
            while process.poll() is None:
                if stop.is_set():
                    if term_sent is None:
                        try:os.killpg(process.pid,signal.SIGTERM)
                        except ProcessLookupError:pass
                        term_sent=time.monotonic()
                    elif time.monotonic()-term_sent>15:
                        try:os.killpg(process.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                time.sleep(.25)
            with guard:running.pop(job['id'],None)
            row.update(returncode=process.returncode,wall_seconds=time.monotonic()-started,finished_utc=now())
        row['status']='interrupted' if stop.is_set() else 'failed_retained'
        if process.returncode==0 and not stop.is_set():
            try:row.update(result=check_result(job,plan['budget_sha256']),status='processed')
            except (ValueError,KeyError,OSError) as error:row['error']=f'{type(error).__name__}: {error}'
        return row
    def one(job):
        total=resources.acquire(job['id'],stop)
        if total is None:
            return existing.get(job['id'],dict(id=job['id'],status='not_started',reserved_cpus=job['reserved_cpus']))
        try:return execute_one(job,total)
        finally:resources.release(job['id'])
    with guard:save()
    with ThreadPoolExecutor(max_workers=plan['executor_threads']) as pool:
        pending={pool.submit(one,job):job for job in jobs}
        for future in as_completed(pending):
            job=pending[future]
            try:row=future.result()
            except Exception as error:row=dict(id=job['id'],status='failed_retained',error=f'{type(error).__name__}: {error}')
            with guard:rows[job['id']]=row;save()
            print(json.dumps(row,ensure_ascii=False),flush=True)
    stable=sha(budget_path)==plan['budget_sha256'] and all(sha(ROOT/name)==digest for name,digest in budget['source_sha256'].items())
    all_processed=set(rows)=={j['id'] for j in jobs} and all(r['status']=='processed' for r in rows.values())
    final=dict(plan_sha256=sha(planpath),updated_utc=now(),source_unchanged=stable,
        cpu_reservations=resources.snapshot(),
        status='interrupted' if stop.is_set() else 'all_workflows_processed' if stable and all_processed else 'completed_with_retained_failures',
        outcomes=sorted(rows.values(),key=lambda r:r['id']))
    write(ledgerpath,final)
    if final['status']!='all_workflows_processed':raise SystemExit(1)


if __name__=='__main__':main()
