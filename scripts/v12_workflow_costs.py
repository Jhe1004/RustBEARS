"""Read source-bound formal API timings without importing the POSIX launcher.

Timing rows are attempts, not independent statistical replicates. No timing is
imputed for a killed/unstarted call, and cached work is not a full-cost point.
"""
from collections import defaultdict
from datetime import datetime, timedelta
import csv
import json
import math
from pathlib import Path

if __package__:
    from .plot_v12_calibration import (budget_info, decode, job_key, load_run,
        read_json, reject_constant, resolve, sha, verify_workflow_source)
    from .v12_case_evidence import verify_case
else:
    from plot_v12_calibration import (budget_info, decode, job_key, load_run,
        read_json, reject_constant, resolve, sha, verify_workflow_source)
    from v12_case_evidence import verify_case

TIMER='scripts/v12_workflow_timing.py'
MODELS=('DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J')
WORKLOADS=('profile_chi2','profile_plugin','lrt_199_199','six_model_average','maps_5000','trees_200_maps_500')
LABELS=('χ² profiles: d, e, j','Plug-in profiles: d, e, j','LRT: 199 + 199','Six models + average','5000 single-tree maps','200 trees × 500 maps')
PHASES={'fit','profile_asymptotic_chi2','profile_bootstrap_plugin','bootstrap_lrt','model_average','run_bsm','fit_over_trees'}


def timestamp(value):
    result=datetime.fromisoformat(value.replace('Z','+00:00'))
    if result.tzinfo is None:raise ValueError('timezone-aware timing required')
    return result


def verify_plan(plan,ledger,budget,budget_hash,plan_hash):
    """Validate every suite ID and setting; a summary counter is insufficient."""
    if plan.get('schema')!='rustbears.formal_suite.v1' or plan.get('budget_status')!='locked' or plan.get('budget_sha256')!=budget_hash:
        raise ValueError('suite plan is not this locked formal budget')
    if plan.get('script_sha256')!=budget['source_sha256'].get('scripts/run_v12_formal_suite.py'):
        raise ValueError('suite launcher source mismatch')
    if ledger.get('plan_sha256')!=plan_hash or ledger.get('source_unchanged') is not True:
        raise ValueError('suite ledger source/plan binding missing')
    if ledger.get('status')!='all_workflows_processed':
        raise ValueError('every required workflow must be processed before formal rendering')
    jobs=plan['jobs'];ids=[j['id'] for j in jobs];outcomes=ledger['outcomes']
    if len(ids)!=len(set(ids)) or len(jobs)!=len(budget['jobs']):raise ValueError('duplicate or missing suite jobs')
    if sorted(json.dumps(j['settings'],sort_keys=True) for j in jobs)!=sorted(json.dumps(j,sort_keys=True) for j in budget['jobs']):
        raise ValueError('suite settings differ from the locked budget')
    if len(outcomes)!=len(ids) or {r['id'] for r in outcomes}!=set(ids) or any(r.get('status')!='processed' for r in outcomes):
        raise ValueError('suite ledger lacks one processed outcome per required ID')
    for job in jobs:
        if job['analysis']!=job['settings']['analysis']:raise ValueError('suite analysis mismatch')
        if type(job['reserved_cpus']) is not int or not 1<=job['reserved_cpus']<=plan['cpu_budget']:raise ValueError('invalid worker reservation')
        argv=job['argv']
        if job['analysis']!='v4' and int(argv[argv.index('--workers')+1])!=job['reserved_cpus']:
            raise ValueError('actual worker arguments differ from immutable plan')
    return {row['id']:row for row in outcomes}


def verify_processed_job(job,budget,budget_hash,outcome):
    """Re-read the scientific checkpoints, not just the launcher's counters."""
    folder=Path(job['output']);settings=job['settings'];analysis=job['analysis']
    if Path(job['result']).parent.resolve()!=folder.resolve():raise ValueError('result is outside its job output')
    previous=outcome['result'];files=previous['evidence_files']
    if Path(previous['path']).resolve()!=Path(job['result']).resolve() or sha(previous['path'])!=previous['sha256']:
        raise ValueError('suite result path/digest changed')
    if not files or any(sha(path)!=digest for path,digest in files.items()):raise ValueError('processed evidence changed')
    if analysis in ('d1','d2'):
        checked=verify_case(folder,settings,budget,budget_hash)
        manifest=read_json(folder/'manifest.json')
        return manifest,checked['files']
    run=load_run(dict(analysis=analysis,path=str(folder)),folder,budget_hash)
    manifest=run['manifest'];metadata=manifest['budget'] if analysis=='v3' else manifest['locked_budget']
    if metadata['payload']!=budget or run['key']!=job_key(settings):raise ValueError('run budget/task mismatch')
    verify_workflow_source(manifest,analysis,budget)
    actual=manifest['settings'] if analysis=='v3' else manifest['execution_settings'] if analysis=='v4' else manifest
    if analysis=='v4':
        if actual!=settings['settings'] or manifest['seed']!=settings['seed']:raise ValueError('V4 settings mismatch')
    elif any(actual.get(k)!=v for k,v in settings.items()):raise ValueError('run settings mismatch')
    return manifest,run['files']


def read_timings(path,job,manifest,helper_hash,earliest):
    """Validate rows from a separately verified formal job, retaining failures."""
    path=Path(path)
    if not path.exists():return []  # Missing observations stay missing, never zero.
    raw=path.read_bytes()
    if raw and not raw.endswith(b'\n'):raise ValueError('truncated timing ledger; preserve and investigate its tail')
    rows=[];seen=set();invocations={}
    resolved=manifest.get('resolved_input',decode(manifest.get('resolved_reference_model',{})))
    expected=dict(analysis=job['analysis'],workers=job['reserved_cpus'],
                  n_tips=manifest.get('n_tips',manifest.get('tip_count')),
                  areas=manifest['areas'] if isinstance(manifest['areas'],int) else len(manifest['areas']),
                  states=len(resolved['states']))
    for line,encoded in enumerate(raw.splitlines(),1):
        row=json.loads(encoded,parse_constant=reject_constant)
        if row.get('mode')!='formal':raise ValueError('smoke/pilot timing cannot become a formal cost')
        if row.get('telemetry_source_sha256')!=helper_hash:raise ValueError('timer helper source mismatch')
        if any(row.get(k)!=v for k,v in expected.items()):raise ValueError('timing dimensions/workers differ from the verified job')
        native_threads=job['reserved_cpus'] if row.get('phase')=='run_bsm' else 1
        if row.get('native_threads')!=native_threads:raise ValueError('native thread metadata differs from the API phase')
        if any(row.get('thread_environment',{}).get(k)!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')):
            raise ValueError('formal timing lacks fixed single-thread BLAS environment')
        if row.get('status') not in ('returned','raised') or row.get('phase') not in PHASES:raise ValueError('unknown API timing phase/status')
        if row['phase']=='fit' and row.get('model') not in MODELS:raise ValueError('unregistered initial model')
        if row['phase']=='model_average' and row.get('n_models')!=6:raise ValueError('model-average workload is not six models')
        if job['analysis']=='v1' and row['phase'].startswith('profile_') and row.get('params')!=['d','e','j']:
            raise ValueError('V1 timing must cover all three parameters')
        if job['analysis']=='v2b' and row['phase']=='bootstrap_lrt' and row.get('route') not in ('A','C','C_to_A'):
            raise ValueError('unregistered bootstrap route')
        if any(not isinstance(row.get(k),(int,float)) or isinstance(row[k],bool) or not math.isfinite(row[k]) or row[k]<0 for k in ('wall_seconds','process_cpu_seconds')):
            raise ValueError('invalid elapsed time')
        if row['wall_seconds']==0:raise ValueError('zero wall time cannot supply a measured cost')
        before,after=row.get('checkpoint_rows_before'),row.get('checkpoint_rows_after')
        if type(before) is not int or type(after) is not int or not 0<=before<=after:raise ValueError('invalid checkpoint row counts')
        start=timestamp(row['started_utc'])
        if start<earliest:raise ValueError('timing predates formal execution')
        invocation=row['invocation'];signature=tuple(row.get(k) for k in ('host','platform','pid','workers'))
        if not row.get('host') or not invocation:raise ValueError('timing host/invocation missing')
        if invocation in invocations and invocations[invocation]!=signature:raise ValueError('mixed host/process within one invocation')
        invocations[invocation]=signature
        if job['analysis'] in ('v1','v2b','v2c'):
            if row.get('scenario')!=manifest['scenario'] or type(row.get('replicate')) is not int or not 1<=row['replicate']<=manifest['outer_replicates']:
                raise ValueError('timing scenario/replicate not in the planned job')
        elif row.get('case_id')!=manifest['case_id']:raise ValueError('timing case ID mismatch')
        identity=(invocation,row['phase'],row.get('replicate'),row.get('model'),row.get('contract'),row.get('route'),row.get('null'),row.get('alt'))
        if identity in seen:raise ValueError('duplicated phase in one invocation/task')
        seen.add(identity)
        if row['phase']=='fit' and row.get('call_workers')!=1:raise ValueError('initial fit worker scope differs')
        if row['phase']=='profile_bootstrap_plugin' and row.get('plugin_bootstrap')!=manifest['plugin_bootstrap']:raise ValueError('plugin B mismatch')
        if row['phase']=='bootstrap_lrt':
            nboot=manifest['inner_bootstrap'] if job['analysis']=='v2b' else manifest['n_boot']
            nalt=manifest['surrogate_alt_replicates'] if job['analysis']=='v2b' else 0
            if row.get('n_boot')!=nboot or row.get('n_alt_sim',0)!=nalt:raise ValueError('bootstrap workload mismatch')
            if row['status']=='returned' and (after!=nboot+nalt or row.get('n_planned')!=nboot):raise ValueError('returned bootstrap is incomplete')
        if row['phase']=='run_bsm' and (row.get('planned_maps')!=manifest['bsm_maps'] or (row['status']=='returned' and row.get('num_maps')!=manifest['bsm_maps'])):
            raise ValueError('single-tree mapping workload mismatch')
        if row['phase']=='fit_over_trees':
            if row.get('planned_trees')!=manifest['max_trees'] or row.get('maps_per_tree')!=manifest['bsm_maps_per_tree']:raise ValueError('tree/map workload mismatch')
            if row['status']=='returned' and (after!=manifest['max_trees'] or row.get('n_trees_planned')!=after or row.get('n_trees_completed',0)+row.get('n_trees_incomplete',0)!=after):raise ValueError('returned tree workflow is incomplete')
        rows.append(dict(row,job_id=job['id'],source_path=str(path),source_line=line,
                         row_id=f'{job["id"]}:{line}',api_complete_fresh=row['status']=='returned' and before==0,
                         classification='raised' if row['status']=='raised' else 'partial_resume' if before else 'complete_fresh'))
    # API calls in this recorder are sequential within a process invocation.
    # This also prevents summing overlapping calls into a composite V2c cost.
    for invocation in invocations:
        group=sorted((r for r in rows if r['invocation']==invocation),key=lambda r:r['started_utc'])
        for left,right in zip(group,group[1:]):
            if timestamp(left['started_utc'])+timedelta(seconds=left['wall_seconds'])>timestamp(right['started_utc'])+timedelta(seconds=.01):
                raise ValueError('overlapping API calls within one invocation')
    return rows


def workload(row):
    analysis,phase=row['analysis'],row['phase']
    if analysis=='v1' and row.get('params')==['d','e','j']:
        if phase=='profile_asymptotic_chi2':return 'profile_chi2'
        if phase=='profile_bootstrap_plugin':return 'profile_plugin'
    if analysis=='v2b' and phase=='bootstrap_lrt' and row.get('n_boot')==199 and row.get('n_alt_sim')==199:return 'lrt_199_199'
    if analysis in ('d1','d2') and phase=='run_bsm' and row.get('planned_maps')==5000:return 'maps_5000'
    if analysis=='d2' and phase=='fit_over_trees' and row.get('planned_trees')==200 and row.get('maps_per_tree')==500:return 'trees_200_maps_500'
    return None


def cost_task_id(row,category):
    return json.dumps([row['job_id'],category,row.get('replicate'),row.get('route'),
                       row.get('contract') if row['analysis'] in ('d1','d2') else None],separators=(',',':'))


def planned_cost_tasks(job,manifest):
    analysis=job['analysis'];result=[]
    def add(category,**identity):
        row=dict(job_id=job['id'],analysis=analysis,**identity)
        result.append(dict(workload=category,task_id=cost_task_id(row,category),job_id=job['id'],**identity))
    if analysis in ('v1','v2b','v2c'):
        for rep in range(1,manifest['outer_replicates']+1):
            if analysis=='v1':
                for category in WORKLOADS[:2]:add(category,replicate=rep)
            elif analysis=='v2c':add('six_model_average',replicate=rep)
            elif manifest['inner_bootstrap']==199 and manifest['surrogate_alt_replicates']==199:
                for route in ('A','C','C_to_A'):add('lrt_199_199',replicate=rep,route=route)
    elif analysis in ('d1','d2'):
        for contract in manifest['likelihood_contracts']:
            if manifest['bsm_maps']==5000:add('maps_5000',contract=contract)
            if analysis=='d2' and manifest['max_trees']==200 and manifest['bsm_maps_per_tree']==500:add('trees_200_maps_500',contract=contract)
    return result


def cost_points(rows):
    """Return all eligible full attempts; never choose a fastest or final attempt."""
    points=[];groups=defaultdict(list)
    for row in rows:
        category=workload(row)
        if category and row['api_complete_fresh']:
            points.append(dict(row,workload=category,component_row_ids=[row['row_id']],component_count=1))
        if row['analysis']=='v2c':groups[(row['job_id'],row['invocation'],row['replicate'])].append(row)
    for group in groups.values():
        fits=[r for r in group if r['phase']=='fit'];average=[r for r in group if r['phase']=='model_average']
        if len(fits)!=6 or len(average)!=1 or {r.get('model') for r in fits}!=set(MODELS):continue
        components=fits+average
        if not all(r['api_complete_fresh'] for r in components) or not all(r.get('converged') is True for r in fits):continue
        if any(r.get('contract')!='tip_survival' for r in fits) or average[0].get('n_models')!=6:continue
        base=average[0]
        points.append(dict(base,workload='six_model_average',phase='six_fits_plus_aggregation',contract='tip_survival',
            wall_seconds=sum(r['wall_seconds'] for r in components),process_cpu_seconds=sum(r['process_cpu_seconds'] for r in components),
            component_row_ids=[r['row_id'] for r in components],component_count=7))
    return [dict(point,task_id=cost_task_id(point,point['workload'])) for point in points]


def inspect_suite(inventory_path):
    """Only complete, locked, source-bound full suites enter the render path."""
    path=Path(inventory_path).resolve();inventory=read_json(path)
    if inventory.get('schema')!='rustbears.workflow_costs.inputs.v1' or inventory.get('formal') is not True:
        raise ValueError('formal=true workflow inventory required; no smoke promotion')
    envelope=dict(inventory,schema='rustbears.figure2.inputs.v1')
    budget,budget_path,_=budget_info(envelope,path.parent);budget_hash=sha(budget_path)
    if set(j['analysis'] for j in budget['jobs'])!={'v1','v2b','v2c','v3','v4','d1','d2'}:raise ValueError('full seven-analysis budget required')
    suite=resolve(path.parent,inventory['suite']);plan_path=suite/'suite_plan.json';ledger_path=suite/'suite_ledger.json'
    plan,ledger=read_json(plan_path),read_json(ledger_path)
    outcomes=verify_plan(plan,ledger,budget,budget_hash,sha(plan_path))
    if not budget['source_sha256'].get(TIMER):raise ValueError('timer source is not frozen in budget')
    files={str(p):sha(p) for p in (path,budget_path,plan_path,ledger_path)}
    for name in ('gate_evidence','design_evidence'):
        if name in budget:
            p=resolve(budget_path.parent,budget[name]['path']);files[str(p)]=sha(p)
    roots={Path(j['cwd']) for j in plan['jobs']}
    if len(roots)!=1:raise ValueError('suite mixes execution source roots')
    root=roots.pop()
    for name,digest in budget['source_sha256'].items():
        source=root/name
        if sha(source)!=digest:raise ValueError('frozen execution source changed: '+name)
        files[str(source)]=digest
    rows=[];jobs=[];planned=[]
    for job in plan['jobs']:
        folder=Path(job['output'])
        if folder.parent.resolve()!=suite or folder.name!=job['id'] or folder.is_symlink():raise ValueError('suite output ownership mismatch')
        manifest,evidence=verify_processed_job(job,budget,budget_hash,outcomes[job['id']]);files.update(evidence)
        timed=[]
        if job['analysis'] not in ('v3','v4'):
            started=read_json(folder/'execution_started.json')
            timed=read_timings(folder/'workflow_timings.jsonl',job,manifest,budget['source_sha256'][TIMER],timestamp(started['first_started_utc']))
            if (folder/'workflow_timings.jsonl').exists():files[str(folder/'workflow_timings.jsonl')]=sha(folder/'workflow_timings.jsonl')
            rows.extend(timed)
            planned.extend(planned_cost_tasks(job,manifest))
        jobs.append(dict(job_id=job['id'],analysis=job['analysis'],settings=job['settings'],timing_rows=len(timed),
            timing_status='not_instrumented_for_panel_c' if job['analysis'] in ('v3','v4') else 'observed' if timed else 'unavailable'))
    points=cost_points(rows)
    available={p['task_id'] for p in points}
    if len({p['task_id'] for p in planned})!=len(planned) or available-{p['task_id'] for p in planned}:raise ValueError('unregistered cost task')
    task_rows=[dict(p,full_cost_attempts=sum(x['task_id']==p['task_id'] for x in points),
                    status='full_cost_observed' if p['task_id'] in available else 'unavailable_full_cost') for p in planned]
    return dict(status='eligible_formal_workflow_costs',formal=True,source_files=files,budget_sha256=budget_hash,
        rows=rows,points=points,jobs=jobs,cost_tasks=task_rows,missingness='Absent or interrupted calls have no imputed elapsed time; incomplete resumes are retained outside full-cost points.',
        suite_wall_time='Not computed from API phase durations; simultaneous jobs overlap.',
        category_counts={key:dict(full_attempts=sum(p['workload']==key for p in points),
            unique_tasks_observed=sum(p['workload']==key and p['status']=='full_cost_observed' for p in task_rows),
            planned_tasks=sum(p['workload']==key for p in task_rows)) for key in WORKLOADS})


def export_costs(audit,folder):
    folder=Path(folder)
    if audit.get('status')!='eligible_formal_workflow_costs' or audit.get('formal') is not True:raise ValueError('formal eligibility required')
    for key,name in (('rows','workflow_all_attempts.jsonl'),('points','workflow_complete_costs.jsonl'),('jobs','workflow_job_denominators.jsonl'),('cost_tasks','workflow_cost_task_denominators.jsonl')):
        (folder/name).write_text(''.join(json.dumps(r,sort_keys=True,ensure_ascii=False,allow_nan=False)+'\n' for r in audit[key]))
    columns=['workload','job_id','analysis','scenario','case_id','replicate','invocation','host','n_tips','areas','states','workers','contract','route','plugin_bootstrap','n_boot','n_alt_sim','wall_seconds','process_cpu_seconds','n_failed','n_trees_incomplete','component_count','component_row_ids']
    with (folder/'workflow_complete_costs.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        for point in audit['points']:
            writer.writerow({k:json.dumps(point[k]) if isinstance(point.get(k),(list,dict)) else point.get(k) for k in columns})
    (folder/'workflow_input_audit.json').write_text(json.dumps({k:v for k,v in audit.items() if k not in ('rows','points')},indent=2,ensure_ascii=False)+'\n')
