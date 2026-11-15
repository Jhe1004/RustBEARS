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
    from .plot_v12_calibration import (budget_info, checkpoint, decode, job_key, load_run,
        read_json, reject_constant, resolve, sha, verify_workflow_source)
    from .v12_case_evidence import verify_case
else:
    from plot_v12_calibration import (budget_info, checkpoint, decode, job_key, load_run,
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
    columns=['workload','job_id','task_id','row_id','analysis','scenario','case_id','replicate','invocation','host','n_tips','areas','states','workers','contract','route','plugin_bootstrap','n_boot','n_alt_sim','status','process_status','wall_seconds','process_cpu_seconds','n_failed','n_trees_incomplete','profile_endpoint_diagnostics','profile_rebased_mle','profile_box_censored_endpoints','profile_disconnected_parameters','profile_missing_se_parameters','profile_artifacts','component_count','component_row_ids']
    with (folder/'workflow_complete_costs.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        for point in audit['points']:
            writer.writerow({k:json.dumps(point[k]) if isinstance(point.get(k),(list,dict)) else point.get(k) for k in columns})
    if 'stratum_counts' in audit:
        with (folder/'workflow_stratum_counts.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=['workload','states','workers','unique_tasks_observed','planned_tasks','full_attempts','unavailable_tasks'])
            writer.writeheader();writer.writerows(audit['stratum_counts'])
        with (folder/'workflow_profile_diagnostics.csv').open('w',newline='') as stream:
            columns=['task_id','row_id','workload','states','workers','replicate','status','process_status','api_complete_fresh','profile_rebased_mle','param','interval_kind','mle','lower','upper','status_lower','status_upper','onesided_upper_from_zero','critical_value','disconnected_support_detected','endpoint_deviance_residual_max','se_hessian_local','se_hessian_missing']
            writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
            for row in audit['rows']:
                for diagnostic in row.get('profile_endpoint_diagnostics',[]):
                    values=dict(row,**diagnostic)
                    writer.writerow({k:values.get(k) for k in columns})
    (folder/'workflow_input_audit.json').write_text(json.dumps({k:v for k,v in audit.items() if k not in ('rows','points')},indent=2,ensure_ascii=False)+'\n')


def profile_cost_diagnostics(directory,evidence_files,api_status):
    """Read returned-profile diagnostics only from the already checked artifacts.

    Box censoring, disconnected support and an unavailable Hessian SE do not
    turn an observed API duration into a numerical failure. In this producer,
    numerical certification failures raise instead of returning an interval.
    """
    if api_status!='returned':raise ValueError('profile diagnostics require a returned API call')
    directory=Path(directory);artifacts={}
    for name in ('intervals.csv','profile_result.json'):
        path=directory/name;digest=evidence_files.get(str(path))
        if digest is None or sha(path)!=digest:raise ValueError('profile diagnostic artifact is not hash-bound: '+name)
        artifacts[str(path)]=digest
    with (directory/'intervals.csv').open(newline='') as stream:intervals=list(csv.DictReader(stream))
    if len(intervals)!=3 or {r['param'] for r in intervals}!={'d','e','j'}:raise ValueError('profile diagnostics omit a parameter')
    result=read_json(directory/'profile_result.json')
    if type(result.get('rebased_mle')) is not bool:raise ValueError('profile rebase diagnostic must be boolean')
    def number(value,nullable=False):
        if nullable and value=='':return None
        try:parsed=float(value)
        except (ValueError,TypeError):raise ValueError('invalid profile diagnostic number') from None
        if not math.isfinite(parsed):raise ValueError('nonfinite profile diagnostic number')
        return parsed
    diagnostics=[]
    for record in sorted(intervals,key=lambda r:('d','e','j').index(r['param'])):
        diagnostic={'param':record['param'],'interval_kind':record['interval_kind'],'mle':number(record['mle']),
            'onesided_upper_from_zero':number(record['onesided_upper_from_zero'],nullable=True),
            'critical_value':number(record['critical_value'])}
        for side in ('lower','upper'):
            status=record['status_'+side]
            if status not in ('interior_bracketed','infeasible_zero_approached','at_physical_zero','censored_at_box_'+side):
                raise ValueError('unknown or uncertified returned-profile endpoint status: '+status)
            value=number(record[side],nullable=True)
            if (value is None)!=(status=='censored_at_box_'+side):raise ValueError('profile endpoint/status mismatch')
            diagnostic.update({side:value,'status_'+side:status})
        if record['disconnected_support_detected'] not in ('True','False'):raise ValueError('invalid disconnected-support flag')
        diagnostic['disconnected_support_detected']=record['disconnected_support_detected']=='True'
        diagnostic['endpoint_deviance_residual_max']=number(record['endpoint_deviance_residual_max'])
        if diagnostic['endpoint_deviance_residual_max']<0:raise ValueError('negative profile endpoint residual')
        diagnostic['se_hessian_local']=number(record['se_hessian_local'],nullable=True)
        if diagnostic['se_hessian_local'] is not None and diagnostic['se_hessian_local']<=0:raise ValueError('invalid Hessian SE')
        diagnostic['se_hessian_missing']=diagnostic['se_hessian_local'] is None
        diagnostics.append(diagnostic)
    return dict(profile_endpoint_diagnostics=diagnostics,profile_rebased_mle=result['rebased_mle'],
        profile_box_censored_endpoints=sum(r['status_'+side]=='censored_at_box_'+side for r in diagnostics for side in ('lower','upper')),
        profile_disconnected_parameters=sum(r['disconnected_support_detected'] for r in diagnostics),
        profile_missing_se_parameters=sum(r['se_hessian_missing'] for r in diagnostics),profile_artifacts=artifacts)


def dedicated_stratum_counts(settings,tasks,points,workloads):
    """Count unique planned tasks on the locked grid, including wholly absent strata."""
    expected={(phase,state,workers,rep) for phase in workloads for state in settings['states']
              for workers in settings['workers'] for rep in range(1,settings['repeats']+1)}
    keys=[(r['workload'],r['states'],r['workers'],r['replicate']) for r in tasks]
    task_by_id={r['task_id']:r for r in tasks}
    if len(keys)!=len(expected) or set(keys)!=expected or len(task_by_id)!=len(tasks):raise ValueError('incomplete or duplicate planned cost strata')
    attempts=defaultdict(int)
    for point in points:
        task=task_by_id.get(point['task_id'])
        if task is None or any(point[k]!=task[k] for k in ('workload','states','workers','replicate')):
            raise ValueError('cost point is outside its planned stratum')
        attempts[point['task_id']]+=1
    if any(r['full_cost_attempts']!=attempts[r['task_id']] for r in tasks):raise ValueError('task cost-attempt count disagrees with observations')
    rows=[]
    for phase in workloads:
        for state in settings['states']:
            for workers in settings['workers']:
                selected=[r for r in tasks if (r['workload'],r['states'],r['workers'])==(phase,state,workers)]
                observed=sum(attempts[r['task_id']]>0 for r in selected)
                rows.append(dict(workload=phase,states=state,workers=workers,unique_tasks_observed=observed,
                    planned_tasks=len(selected),full_attempts=sum(attempts[r['task_id']] for r in selected),
                    unavailable_tasks=len(selected)-observed))
    return rows


def inspect_dedicated_benchmark(inventory_path):
    """Validate the independent 16/64/256-state cost protocol, not calibration."""
    if __package__:
        from . import benchmark_v12_workflows as bench
    else:
        import benchmark_v12_workflows as bench
    path=Path(inventory_path).resolve();inventory=read_json(path)
    core=None
    if inventory.get('schema')=='rustbears.workflow_benchmark.inputs.v2':
        if __package__:
            from .v12_core_reporting import core_cost_context
        else:
            from v12_core_reporting import core_cost_context
        core=core_cost_context(path,inventory,bench)
    if not core and (inventory.get('schema')!='rustbears.workflow_benchmark.inputs.v1' or inventory.get('formal') is not True):
        raise ValueError('formal dedicated benchmark inventory required; no smoke promotion')
    budget_path=resolve(path.parent,inventory['budget']['path'])
    if sha(budget_path)!=inventory['budget']['sha256']:raise ValueError('dedicated budget SHA mismatch')
    if core:budget=core['budget']
    else:budget,_=bench.load_locked_budget(budget_path,source_script=bench.__file__,execute=True)
    folder=resolve(path.parent,inventory['benchmark']);plan_path=folder/'plan.json';result_path=folder/'results.json'
    plan=read_json(plan_path);result=read_json(result_path)
    if plan.get('mode')!='formal' or plan.get('formal') is not True or result.get('formal') is not True:
        raise ValueError('smoke or pilot benchmarks cannot supply formal workflow costs')
    if result.get('status')!='all_planned_phases_processed' or result.get('source_unchanged') is not True or result.get('plan_sha256')!=sha(plan_path):
        raise ValueError('dedicated benchmark is not completely processed with unchanged sources')
    if timestamp(result['first_started_utc'])<timestamp(budget['locked_utc']):raise ValueError('benchmark execution predates its budget lock')
    if core and timestamp(result['first_started_utc'])<timestamp(core['specified_utc']):raise ValueError('core design postdates benchmark execution')
    if plan.get('budget',{}).get('sha256')!=sha(budget_path) or plan['budget']['payload']!=budget:raise ValueError('dedicated benchmark budget binding mismatch')
    selected=bench.select_budget_job(budget,analysis='workflow_benchmark');settings={k:v for k,v in selected.items() if k!='analysis'}
    if plan['settings']!=settings or plan['cells']!=bench.cells_for(settings) or plan['phases']!=list(bench.PHASES):raise ValueError('dedicated benchmark dimensions differ from locked design')
    predecessor,digest=bench.predecessor_budget(budget,budget_path)
    predecessor_checked=core['predecessor_checked'] if core else bench.completed_statistical_suite(plan['after_suite'],predecessor,digest)
    if result.get('after_suite')!=predecessor_checked:raise ValueError('predecessor completion evidence changed')
    if not core:bench.verify_sources(plan['sources'])
    if any(budget['source_sha256'].get(k)!=v for k,v in plan['sources']['files'].items()):raise ValueError('benchmark source not frozen in budget')
    source_root=core['root'] if core else bench.ROOT
    if core:
        expected_sources=set(bench.SOURCE_FILES)|{p.relative_to(source_root).as_posix() for p in (source_root/'rustbears').rglob('*.py')}
        if set(plan['sources']['files'])!=expected_sources:raise ValueError('dedicated source manifest is incomplete')
    native=Path(plan['sources']['native_path']);native_key=native.relative_to(source_root).as_posix() if native.is_relative_to(source_root) else None
    if (budget.get('native_sha256') or budget['source_sha256'].get(native_key))!=sha(native):raise ValueError('native library not frozen in dedicated budget')
    if core and (native_key is None or plan['sources']['native_sha256']!=sha(native)):raise ValueError('dedicated native path or recorded SHA differs from execution source')
    expected={(c['id'],p) for c in plan['cells'] for p in bench.PHASES};records=result['records']
    keys=[(r['cell_id'],r['phase']) for r in records]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:raise ValueError('missing or duplicated dedicated workload task')
    by_key=dict(zip(keys,records));files={str(p):sha(p) for p in (path,budget_path,plan_path,result_path,native)}
    if core:files.update(core['files'])
    for name,digest in plan['sources']['files'].items():files[str(source_root/name)]=digest
    for item in (budget['gate_evidence'],budget['predecessor_budget']):
        p=resolve(budget_path.parent,item['path']);files[str(p)]=sha(p)
    for name in ('suite_plan.json','suite_ledger.json'):
        p=Path(plan['after_suite'])/name;files[str(p)]=sha(p)
    rows=[];points=[];tasks=[]
    for cell in plan['cells']:
        cell_dir=folder/cell['id'];input_path=cell_dir/'input_manifest.json';inputs=read_json(input_path)
        files[str(input_path)]=sha(input_path)
        for name,digest in inputs['files'].items():
            p=cell_dir/'inputs'/name
            if sha(p)!=digest:raise ValueError('dedicated input changed')
            files[str(p)]=digest
        for phase in bench.PHASES:
            record=by_key[(cell['id'],phase)];task_id=f'{cell["id"]}:{phase}';available=0
            if record['status'] not in ('complete','timeout_retained','failed_retained','blocked_initial_fit_failure'):
                raise ValueError('unprocessed dedicated phase')
            if record['status']=='blocked_initial_fit_failure':
                if phase in ('six_model_average','fit_over_trees') or by_key[(cell['id'],'six_model_average')]['status']=='complete':raise ValueError('unjustified blocked phase')
            else:
                if not record.get('attempts') or not record.get('files'):raise ValueError('phase lacks attempt evidence')
                if record['status']!=record['attempts'][-1]['status']:raise ValueError('phase and final attempt status disagree')
                actual={str(p):sha(p) for p in (cell_dir/phase).rglob('*') if p.is_file() and 'runtime_cache' not in p.parts}
                if actual!=record['files']:raise ValueError('dedicated attempt evidence changed or includes unhashed artifacts')
                files.update(actual)
                for index,attempt in enumerate(record['attempts'],1):
                    directory=Path(attempt['directory'])
                    if directory.resolve()!=cell_dir/phase/f'attempt_{index:03d}' or directory.is_symlink():raise ValueError('attempt ownership mismatch')
                    spec=read_json(directory/'spec.json')
                    if spec['mode']!='formal' or spec['cell']!=cell or spec['settings']!=settings or spec['sources']!=plan['sources'] or spec['input_sha256']!=inputs['files']:
                        raise ValueError('phase spec differs from its locked inputs')
                    if timestamp(spec['started_utc'])<timestamp(result['first_started_utc']):raise ValueError('phase predates benchmark execution')
                    timing=directory/'workflow_timings.jsonl';observed=[]
                    if timing.exists():
                        raw=timing.read_bytes()
                        if not raw.endswith(b'\n'):raise ValueError('truncated dedicated timing')
                        observed=[json.loads(line,parse_constant=reject_constant) for line in raw.splitlines()]
                        if len(observed)!=1:raise ValueError('each phase process must record one API call')
                    complete_attempt=attempt['status']=='complete'
                    if complete_attempt and (attempt['returncode']!=0 or attempt['timed_out'] or attempt['interrupted'] or not observed):raise ValueError('false completed phase process')
                    if complete_attempt:
                        returned=read_json(directory/'worker_result.json')
                        if returned['status']!='returned' or returned['formal'] is not True or returned.get('source_unchanged') is not True:raise ValueError('missing successful formal worker result')
                        count=returned['fit_calls_started'];calls_path=directory/'fit_calls.jsonl'
                        calls=[] if not calls_path.exists() else [json.loads(line,parse_constant=reject_constant) for line in calls_path.read_text().splitlines()]
                        starts=[r['id'] for r in calls if r['event']=='started'];finishes=[r['id'] for r in calls if r['event']=='finished']
                        if sorted(starts)!=list(range(1,count+1)) or sorted(finishes)!=list(range(1,count+1)):
                            raise ValueError('complete cost is missing actual fit-call audit events')
                        if any(r['native_threads']!=1 or r['maxiter']!=settings['maxiter'] for r in calls):raise ValueError('actual refit thread/optimizer budget differs')
                    for timing_row in observed:
                        if (timing_row.get('mode')!='formal' or timing_row['phase']!=phase or timing_row.get('cell_id')!=cell['id']
                            or timing_row.get('status') not in ('returned','raised') or timing_row['host']!=plan['host']
                            or timing_row['telemetry_source_sha256']!=plan['sources']['files'][TIMER]
                            or any(timing_row.get(k)!=cell[k] for k in ('states','workers','areas','replicate'))
                            or timing_row['n_tips']!=settings['n_tips']
                            or timing_row['native_threads']!=(cell['workers'] if phase=='run_bsm' else 1)
                            or any(timing_row['thread_environment'].get(k)!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))):
                            raise ValueError('dedicated timing provenance mismatch')
                        if not math.isfinite(timing_row['wall_seconds']) or timing_row['wall_seconds']<=0:raise ValueError('invalid API cost')
                        before=timing_row['checkpoint_rows_before'];fresh=complete_attempt and timing_row['status']=='returned' and before==0
                        row=dict(timing_row,job_id=cell['id'],task_id=task_id,workload=phase,row_id=f'{task_id}:{index}',
                            source_path=str(timing),source_line=1,process_status=attempt['status'],api_complete_fresh=fresh,
                            classification='complete_fresh' if fresh else 'partial_resume' if before else attempt['status'],
                            process_wall_seconds=attempt['process_wall_seconds'],process_peak_rss_bytes=attempt['process_peak_rss_bytes'])
                        if complete_attempt and phase.startswith('profile_'):
                            row.update(profile_cost_diagnostics(directory,actual,timing_row['status']))
                        rows.append(row)
                        if fresh:
                            verify_dedicated_phase(directory,cell_dir/phase,phase,settings,timing_row)
                            points.append(dict(row,component_count=7 if phase=='six_model_average' else 1,component_row_ids=[row['row_id']]))
                            available+=1
                    if not observed:
                        rows.append(dict(job_id=cell['id'],task_id=task_id,workload=phase,process_status=attempt['status'],
                            host=plan['host'],analysis='workflow_benchmark',states=cell['states'],workers=cell['workers'],
                            areas=cell['areas'],replicate=cell['replicate'],n_tips=settings['n_tips'],api_complete_fresh=False,
                            classification='unavailable_terminal_API_timing',wall_seconds=None,
                            process_wall_seconds=attempt['process_wall_seconds'],api_right_censored_elapsed_seconds=attempt['api_right_censored_elapsed_seconds']))
            tasks.append(dict(job_id=cell['id'],task_id=task_id,workload=phase,states=cell['states'],workers=cell['workers'],replicate=cell['replicate'],full_cost_attempts=available,
                status='full_cost_observed' if available else 'unavailable_full_cost',scientific_process_status=record['status']))
    order=('profile_chi2','profile_plugin','bootstrap_lrt','six_model_average','run_bsm','fit_over_trees')
    labels=['χ² profiles: d, e, j',f'Plug-in profiles: B={settings["n_boot_plugin"]}',f'LRT: B={settings["n_boot"]}',
        'Six models + average',f'{settings["bsm_maps"]} single-tree maps',f'{settings["n_trees"]} trees × {settings["maps_per_tree"]} maps']
    audit=dict(status='eligible_formal_workflow_costs',formal=True,source_files=files,budget_sha256=sha(budget_path),rows=rows,points=points,
        jobs=records,cost_tasks=tasks,workload_order=order,workload_labels=labels,settings=settings,protocol='dedicated_computational_workloads',
        stratum_counts=dedicated_stratum_counts(settings,tasks,points,order),
        full_proposal_grid=plan['full_proposal_grid'],category_counts={p:dict(full_attempts=sum(x['workload']==p for x in points),
            unique_tasks_observed=sum(x['workload']==p and x['status']=='full_cost_observed' for x in tasks),
            planned_tasks=sum(x['workload']==p for x in tasks)) for p in order})
    if core:
        audit.update(scope='core_v1',repetition_interpretation=core['repetition_interpretation'],
            independent_repetitions=1,between_repeat_variability_estimated=False)
    return audit


def verify_dedicated_phase(directory,shared,phase,settings,timing):
    """Verify task IDs and primary outputs before calling an API cost complete."""
    def table(name):
        with (directory/name).open(newline='') as stream:return list(csv.DictReader(stream))
    if phase=='six_model_average':
        models=table('models.csv');cache=read_json(directory/'fits.json')
        if len(models)!=6 or {r['model'] for r in models}!=set(MODELS) or not math.isclose(sum(float(r['weight']) for r in models),1,abs_tol=1e-10):raise ValueError('incomplete six-model cost')
        if [r['model'] for r in cache['records']]!=list(MODELS) or not all(r['record']['converged'] for r in cache['records']):raise ValueError('six-model initial fit cache incomplete')
        if sha(directory/'fits.npz')!=cache['array_sha256']:raise ValueError('fit cache array mismatch')
    elif phase.startswith('profile_'):
        intervals=table('intervals.csv')
        if len(intervals)!=3 or {r['param'] for r in intervals}!={'d','e','j'}:raise ValueError('profile cost omits a parameter')
        for name in ('d','e','j'):
            if not table(f'curve_{name}.csv'):raise ValueError('profile curve missing')
            if phase=='profile_plugin' and sorted(int(r['rep_id']) for r in table(f'plugin_{name}.csv'))!=list(range(1,settings['n_boot_plugin']+1)):
                raise ValueError('plug-in profile omitted a calibration replicate')
    elif phase=='bootstrap_lrt':
        records,_=checkpoint(shared/'checkpoint.jsonl',settings['n_boot'])
        reps=table('replicates.csv')
        if sorted(int(r['rep_id']) for r in reps)!=list(range(1,settings['n_boot']+1)):raise ValueError('bootstrap cost discarded a replicate')
        if timing['checkpoint_rows_after']!=len(records):raise ValueError('bootstrap timing checkpoint count mismatch')
    elif phase=='fit_over_trees':
        records,_=checkpoint(shared/'checkpoint.jsonl',settings['n_trees'])
        result=read_json(directory/'multitree_result.json')
        if result['n_planned']!=len(records) or result['n_complete']+result['n_incomplete']!=len(records):raise ValueError('multitree cost discarded a tree')
        if len(table('per_tree_models.csv'))!=settings['n_trees']*6:raise ValueError('multitree model grid incomplete')
        if timing['checkpoint_rows_after']!=len(records):raise ValueError('tree timing checkpoint count mismatch')
        _verify_cost_maps(directory,'joint_',settings['n_tips'],result['n_complete']*settings['maps_per_tree'])
    else:_verify_cost_maps(directory,'',settings['n_tips'],settings['bsm_maps'])


def _verify_cost_maps(directory,prefix,tips,planned):
    from collections import Counter
    if planned==0:
        if (directory/(prefix+'map_metadata.csv')).exists():raise ValueError('maps present despite zero complete trees')
        marker=read_json(directory/'no_complete_tree_maps.json')
        if marker.get('status')!='no_complete_tree_maps' or marker.get('n_complete')!=0:raise ValueError('missing all-failed tree mapping marker')
        return
    def table(name):
        with (directory/(prefix+name)).open(newline='') as stream:return list(csv.DictReader(stream))
    metadata=table('map_metadata.csv');ids=[int(r['map']) for r in metadata]
    tuples={(r['tree_id'],r['model_id'],r['map_id']) for r in metadata}
    if sorted(ids)!=list(range(1,planned+1)) or len(tuples)!=planned:raise ValueError('map cost has duplicated or missing IDs')
    keys=('tree_id','tree_topology_hash','model_id','map_id','stratum_weight','precision_mode','likelihood_conditioning')
    by_id={int(r['map']):r for r in metadata}
    for r in metadata:
        if any(not r.get(k) for k in keys) or r['precision_mode']!='strict_f64' or r['likelihood_conditioning']!='tip_survival':raise ValueError('map provenance missing')
    clado=table('cladogenetic_events.csv');counts=Counter(int(r['map']) for r in clado)
    if counts!={i:tips-1 for i in ids}:raise ValueError('cladogenetic event count incomplete')
    for row in clado+table('anagenetic_events.csv'):
        metadata_row=by_id.get(int(row['map']))
        if metadata_row is None or any(row[k]!=metadata_row[k] for k in keys):raise ValueError('event provenance differs from its map')
