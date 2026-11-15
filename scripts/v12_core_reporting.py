"""Explicit core-scope reporting over immutable, independently bound runs.

This module reads evidence only. It never relabels an old checkpoint with a new
budget, resumes a calculation, or treats a prospective design as a passed gate.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import copy
import json


def plot_module():
    if __package__:
        from . import plot_v12_calibration as plot
    else:
        import plot_v12_calibration as plot
    return plot


def need(value,message):
    if not value:raise ValueError(message)


def timestamp(value):
    result=datetime.fromisoformat(value.replace('Z','+00:00'))
    need(result.tzinfo is not None,'timezone-aware design/execution time required')
    return result


def bound_json(ref,base,files):
    p=plot_module();path=p.resolve(base,ref['path'])
    need(p.sha(path)==ref['sha256'],'referenced evidence SHA mismatch: '+str(path))
    files[str(path)]=ref['sha256']
    return p.read_json(path),path


def bound_budget(ref,base,files):
    """Validate each original budget without requiring unrelated jobs to finish."""
    budget,path=bound_json(ref,base,files)
    need(budget.get('schema')=='rustbears.formal_budget.v1' and budget.get('status')=='locked','a locked formal budget is required')
    need(timestamp(budget['locked_utc'])<=datetime.now(timezone.utc),'future budget lock')
    gate,_=bound_json(budget['gate_evidence'],path.parent,files)
    need(gate.get('gate2_passed') is True and gate.get('gate3_passed') is True,'actual Gate 2 and Gate 3 are required')
    need(gate.get('source_commit')==budget.get('source_commit') and gate.get('native_sha256')==budget.get('native_sha256')
         and bool(gate.get('native_sha256')),'gate source/native identity differs from locked budget')
    if 'source_sha256' in gate:
        need(bool(gate['source_sha256']) and all(budget.get('source_sha256',{}).get(k)==v for k,v in gate['source_sha256'].items()),
             'gate source inventory differs from locked budget')
    design,design_path=bound_json(budget['design_evidence'],path.parent,files)
    need(design['jobs']==budget['jobs'],'budget settings differ from its prospective design')
    if design.get('specified_utc'):
        need(timestamp(design['specified_utc'])<=timestamp(budget['locked_utc']),'design postdates budget lock')
    need(isinstance(budget.get('source_commit'),str) and len(budget['source_commit'])==40,'locked budget lacks its execution source identity')
    if design.get('source_commit') is not None:
        need(design['source_commit']==budget['source_commit'],'budget source differs from the source-bound design')
    return budget,path,design,design_path


def source_binding(item,budget,base,files):
    """The execution checkout is separate from this postprocessing checkout."""
    p=plot_module();source=item['source'];root=p.resolve(base,source['root'])
    need(source['commit']==budget['source_commit'],'run source commit differs from its own budget')
    need(source['sha256']==budget['source_sha256'] and bool(source['sha256']),'run source inventory differs from its own budget')
    for name,digest in source['sha256'].items():
        relative=Path(name)
        need(not relative.is_absolute() and '..' not in relative.parts,'source inventory escapes execution root')
        path=root/relative
        need(p.sha(path)==digest,'execution source changed: '+str(path))
        files[str(path)]=digest
    return root


def core_panels():
    p=plot_module();panels=copy.deepcopy(p.CONTRACT['main_panels'])
    panels['a']['tips']=[40]
    panels['b']['scenarios']=['null_e02','alt_j003']
    return panels


def planned_ids(job):
    analysis=job['analysis']
    if analysis in ('v1','v2b','v2c'):
        return dict(kind='outer_replicate',first=1,last=job['outer_replicates'],count=job['outer_replicates'])
    if analysis=='v3':return dict(kind='dataset',ids=list(range(1,job['replicates']+1)),count=job['replicates'])
    if analysis=='v4':return {name:list(range(1,settings['replicates']+1)) for name,settings in job['settings'].items()}
    return dict(kind='complete_case_workflow',case_id=job['case_id'],
        models=['DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J'],
        contracts=['tip_survival','legacy_unconditioned'])


def verify_registry(design):
    p=plot_module();jobs=design['jobs'];registry=design['job_registry']
    need(design.get('schema')=='rustbears.core_statistical_design.v1' and design.get('scope')=='core_v1','explicit core_v1 scientific design required')
    need(len(jobs)==11 and len(registry)==11 and len({r['id'] for r in registry})==11,'core registry must declare eleven distinct fresh jobs')
    need(sorted(r['settings_index'] for r in registry)==list(range(11)),'registry settings indices missing or duplicated')
    need(Counter(j['analysis'] for j in jobs)==Counter(v1=4,v2b=2,v2c=2,v3=1,d1=1,d2=1),'core analysis inventory differs from the authorized design')
    need({(j['n_tips'],j['areas'],j['scenario']) for j in jobs if j['analysis']=='v1'}=={(40,4,s) for s in p.V1_SCENARIOS},'core V1 requires all four 40-tip, four-area boundary cells')
    need({j['scenario'] for j in jobs if j['analysis']=='v2b'}=={'null_e02','alt_j003'},'core V2b scenarios differ from the prospective choice')
    need({j['scenario'] for j in jobs if j['analysis']=='v2c'}=={'higher_jump','lower_jump'},'core V2c requires both jump scenarios')
    need(next(j for j in jobs if j['analysis']=='v3')['replicates']==10 and design['sampling_policy']['v3_dataset_ids']==list(range(1,11)),'core V3 requires the first ten planned datasets')
    result={}
    for row in registry:
        job=jobs[row['settings_index']]
        need(row['analysis']==job['analysis'] and row['expected_task_ids']==planned_ids(job),'registry expected task IDs differ from the complete planned job')
        result[row['id']]=dict(job=job,origin='fresh',expected_task_ids=row['expected_task_ids'])
    inherited=design['inherited_evidence']
    need(len(inherited)==1 and inherited[0]['analysis']=='v4' and inherited[0]['fresh_execution'] is False,'only complete original V4 may be inherited')
    old=inherited[0]
    need(old['id'] not in result,'inherited and fresh job identities collide')
    result[old['id']]=dict(job=old['settings'],origin='inherited',expected_task_ids=planned_ids(old['settings']),inheritance=old)
    return result


def actual_ids(run,expected):
    analysis=run['analysis'];rows=run['records']
    if analysis=='v4':
        actual={section:sorted(r['replicate'] for r in rows if r['section']==section) for section in expected}
        need(set(r['section'] for r in rows)==set(expected) and actual==expected,'inherited V4 task IDs are incomplete or duplicated')
    else:
        key='dataset' if analysis=='v3' else 'id'
        ids=expected['ids'] if key=='dataset' else list(range(expected['first'],expected['last']+1))
        need(sorted(r[key] for r in rows)==ids and len(rows)==expected['count'],'planned outcome IDs are missing or duplicated')


def recorded_sources(manifest,analysis,budget):
    p=plot_module();p.verify_workflow_source(manifest,analysis,budget)
    recorded=manifest.get('production_source_sha256') if analysis=='v3' else manifest.get('source_sha256')
    need(bool(recorded),'run lacks recorded production source hashes')
    for name,digest in recorded.items():
        key='rustbears/'+name if analysis=='v3' else name
        need(budget['source_sha256'].get(key)==digest,'recorded production source differs from original budget')
    if analysis=='v3':need(manifest['native_sha256']==budget['native_sha256'],'V3 native differs from original budget')


def inspect_collection(path):
    p=plot_module();path=Path(path).resolve();inventory=p.read_json(path)
    need(inventory.get('schema')=='rustbears.figure2.inputs.v2' and inventory.get('formal') is True and inventory.get('scope')=='core_v1','explicit formal core collection v2 required; no smoke promotion')
    files={str(path):p.sha(path)}
    design,design_path=bound_json(inventory['core_design'],path.parent,files)
    registry=verify_registry(design);specified=timestamp(design['specified_utc'])
    need(inventory['panel_selection']==core_panels(),'main panels differ from the prospective core selection')
    items=inventory['runs'];ids=[i['id'] for i in items]
    need(len(ids)==len(registry) and len(ids)==len(set(ids)) and set(ids)==set(registry),'every declared fresh and inherited run must occur exactly once')
    runs=[];cases=[];origins=[];paths=set()
    for item in items:
        registered=registry[item['id']];job=registered['job'];analysis=job['analysis']
        need(item['analysis']==analysis and item['origin']==registered['origin'] and item['expected_task_ids']==registered['expected_task_ids'],'run declaration differs from prospective task IDs or origin')
        folder=p.resolve(path.parent,item['path']);need(folder not in paths,'a run folder cannot be counted twice');paths.add(folder)
        budget,budget_path,run_design,run_design_path=bound_budget(item['budget'],path.parent,files)
        need(item['design']['sha256']==p.sha(run_design_path) and p.resolve(path.parent,item['design']['path'])==run_design_path,'run design reference differs from its original budget')
        need(sum(j==job for j in budget['jobs'])==1,'required job not uniquely registered in original budget')
        if registered['origin']=='fresh':
            need(run_design_path==design_path and p.sha(budget_path)==item['budget']['sha256'],'fresh run must use the new core design')
            # The immutable scientific design precedes source freeze and may
            # have source_commit=null; the final budget binds the actual source.
            if design.get('source_commit') is not None:
                need(budget['source_commit']==design['source_commit'],'fresh source differs from locked core source')
        else:
            old=registered['inheritance']
            need(p.sha(budget_path)==old['original_budget']['sha256'] and budget['source_commit']==old['original_source_commit'] and budget['native_sha256']==old['original_native_sha256'],'inherited V4 source/budget/native identity changed')
        source_binding(item,budget,path.parent,files)
        if analysis in ('d1','d2'):
            if __package__:
                from .v12_case_evidence import verify_case
            else:
                from v12_case_evidence import verify_case
            checked=verify_case(folder,job,budget,p.sha(budget_path));files.update(checked['files'])
            manifest=p.read_json(folder/'manifest.json')
            recorded_sources(manifest,analysis,budget)
            cases.append(dict(id=item['id'],analysis=analysis,path=str(folder),status=checked['status'],expected_task_ids=item['expected_task_ids'],budget_sha256=p.sha(budget_path)))
        else:
            run=p.load_run(dict(analysis=analysis,path=str(folder)),path.parent,p.sha(budget_path));manifest=run['manifest']
            metadata=manifest['budget'] if analysis=='v3' else manifest['locked_budget']
            need(metadata['payload']==budget,'run payload differs from its original budget bytes')
            settings=manifest['settings'] if analysis=='v3' else manifest['execution_settings'] if analysis=='v4' else manifest
            need(settings==job['settings'] if analysis=='v4' else all(settings.get(k)==v for k,v in job.items()),'run settings differ from its planned job')
            if analysis=='v4':need(manifest['seed']==job['seed'],'inherited V4 seed changed')
            actual_ids(run,item['expected_task_ids']);recorded_sources(manifest,analysis,budget)
            run.update(job=job,collection_id=item['id'],source_budget_sha256=p.sha(budget_path),origin=item['origin'])
            files.update(run['files']);runs.append(run)
        if registered['origin']=='fresh':
            if analysis=='v3':
                coordinator=item['coordinator_evidence']
                plan,plan_path=bound_json(coordinator['plan'],path.parent,files)
                ledger,_=bound_json(coordinator['ledger'],path.parent,files)
                need(plan['budget_sha256']==p.sha(budget_path) and ledger['plan_sha256']==p.sha(plan_path),'V3 coordinator budget/plan mismatch')
                scheduled=[r for r in plan['jobs'] if r['id']==item['id']]
                completed=[r for r in ledger['outcomes'] if r['id']==item['id']]
                need(len(scheduled)==len(completed)==1 and scheduled[0]['settings']==job
                     and Path(scheduled[0]['output']).resolve()==folder and completed[0]['status']=='processed','V3 coordinator lacks the exact processed task')
                first=completed[0]['started_utc']
            else:first=p.read_json(folder/'execution_started.json')['first_started_utc']
            need(timestamp(first)>=specified,'core design was specified after execution began')
            need(timestamp(first)>=timestamp(budget['locked_utc']),'core execution began before its budget lock')
        origins.append(dict(id=item['id'],analysis=analysis,origin=item['origin'],path=str(folder),budget_sha256=p.sha(budget_path),source_commit=budget['source_commit'],expected_task_ids=item['expected_task_ids']))
    return runs,dict(status='eligible_formal_inputs',formal=True,scope='core_v1',inventory_sha256=p.sha(path),
        core_design_sha256=p.sha(design_path),main_panels=inventory['panel_selection'],jobs=len(registry),
        calibration_jobs=len(runs),planned_outcomes=sum(len(r['records']) for r in runs),
        source_files=files,run_origins=origins,cases=cases,
        exclusions='All twelve declared core collection jobs retained. Expanded original jobs are outside the prospective amendment; the original full suite is not claimed complete.')


def core_predecessor(folder,budget,budget_hash,registry,files):
    """Read the declared core suite terminal ledger, including retained failures."""
    p=plot_module();folder=Path(folder);plan_path=folder/'suite_plan.json';ledger_path=folder/'suite_ledger.json'
    plan=p.read_json(plan_path);ledger=p.read_json(ledger_path)
    files.update({str(x):p.sha(x) for x in (plan_path,ledger_path)})
    expected={key:value['job'] for key,value in registry.items() if value['origin']=='fresh'}
    jobs=plan['jobs'];rows=ledger['outcomes'];ids=[j['id'] for j in jobs]
    need(plan['budget_sha256']==budget_hash and plan['script_sha256']==budget['source_sha256'].get('scripts/run_v12_formal_suite.py'),'core predecessor source/budget mismatch')
    need(len(ids)==len(expected) and len(set(ids))==len(ids) and set(ids)==set(expected),'missing or duplicated core predecessor task')
    need(all(j['settings']==expected[j['id']] for j in jobs),'core predecessor settings mismatch')
    need(sorted(json.dumps(j,sort_keys=True) for j in expected.values())==sorted(json.dumps(j,sort_keys=True) for j in budget['jobs']),'predecessor budget differs from core design')
    need(ledger['plan_sha256']==p.sha(plan_path) and ledger.get('source_unchanged') is True
         and ledger.get('status') in ('all_workflows_processed','completed_with_retained_failures')
         and len(rows)==len(ids) and len({r['id'] for r in rows})==len(rows)
         and {r['id'] for r in rows}==set(ids) and all(r['status'] in ('processed','failed_retained') for r in rows)
         and ledger.get('cpu_reservations',{}).get('active_reserved_cpus')==0,'core predecessor is not fully terminal and idle')
    for row in rows:
        if row['status']=='processed':
            result=row['result'];need(bool(result.get('evidence_files')),'processed predecessor lacks bound evidence')
            for name,digest in {result['path']:result['sha256'],**result['evidence_files']}.items():
                need(p.sha(name)==digest,'predecessor evidence changed');files[str(Path(name))]=digest
    return dict(plan_sha256=p.sha(plan_path),ledger_sha256=p.sha(ledger_path),status=ledger['status'])


def core_cost_context(path,inventory,bench):
    p=plot_module();path=Path(path).resolve();files={str(path):p.sha(path)}
    need(inventory.get('schema')=='rustbears.workflow_benchmark.inputs.v2' and inventory.get('formal') is True
         and inventory.get('scope')=='core_v1','explicit formal core workflow inventory v2 required')
    budget,budget_path,workflow,workflow_path=bound_budget(inventory['budget'],path.parent,files)
    declared,declared_path=bound_json(inventory['workflow_design'],path.parent,files)
    need(declared_path==workflow_path and declared==workflow and workflow.get('schema')=='rustbears.core_workflow_cost_design.v1'
         and workflow.get('scope')=='core_v1','core workflow design mismatch')
    core,core_path=bound_json(inventory['core_design'],path.parent,files);registry=verify_registry(core)
    predecessor,predecessor_path,pred_design,pred_design_path=bound_budget(budget['predecessor_budget'],budget_path.parent,files)
    need(p.sha(predecessor_path)==budget['predecessor_statistical_budget_sha256'] and pred_design_path==core_path
         and predecessor['gate_evidence']['sha256']==budget['gate_evidence']['sha256']
         and predecessor['source_commit']==budget['source_commit']
         and predecessor['source_sha256']==budget['source_sha256'],'core predecessor is not the declared statistical budget/gate/source')
    root=source_binding(inventory,budget,path.parent,files)
    settings={k:v for k,v in bench.select_budget_job(budget,analysis='workflow_benchmark').items() if k!='analysis'}
    need(settings['states']==[16,64,256] and settings['workers']==[1,8] and settings['repeats']==1,'core cost scope requires one repetition of six planned cells')
    cells=bench.cells_for(settings)
    need(workflow['expected_cell_ids']==[c['id'] for c in cells] and workflow['expected_phases']==list(bench.PHASES)
         and workflow['planned_tasks']==36 and inventory['expected_task_ids']==[c['id']+':'+phase for c in cells for phase in bench.PHASES],
         'core cost task IDs differ from the prospective 36-task grid')
    labels=['χ² profiles: d, e, j',f'Plug-in profiles: B={settings["n_boot_plugin"]}',f'LRT: B={settings["n_boot"]}',
            'Six models + average',f'{settings["bsm_maps"]} single-tree maps',f'{settings["n_trees"]} trees × {settings["maps_per_tree"]} maps']
    need(inventory['workload_labels']==labels,'core workload labels differ from the planned payload')
    # Inherited V4 keeps its original source, budget and all planned records.
    old_item=inventory['inherited_v4'];old=core['inherited_evidence'][0]
    need(old_item['id']==old['id'] and old_item['analysis']=='v4' and old_item['origin']=='inherited'
         and old_item['expected_task_ids']==planned_ids(old['settings']),'inherited V4 declaration differs from core design')
    old_budget,old_path,old_design,old_design_path=bound_budget(old_item['budget'],path.parent,files)
    need(p.sha(old_path)==old['original_budget']['sha256'] and old_budget['source_commit']==old['original_source_commit']
         and old_budget['native_sha256']==old['original_native_sha256'],'inherited V4 identity changed')
    need(sum(job==old['settings'] for job in old_budget['jobs'])==1,'inherited V4 job missing or duplicated in original budget')
    need(old_item['design']['sha256']==p.sha(old_design_path) and p.resolve(path.parent,old_item['design']['path'])==old_design_path,'inherited V4 design reference mismatch')
    source_binding(old_item,old_budget,path.parent,files)
    old_run=p.load_run(old_item,path.parent,p.sha(old_path));m=old_run['manifest']
    need(m['locked_budget']['payload']==old_budget and m['execution_settings']==old['settings']['settings']
         and m['seed']==old['settings']['seed'],'inherited V4 settings changed')
    actual_ids(old_run,planned_ids(old['settings']));recorded_sources(m,'v4',old_budget);files.update(old_run['files'])
    folder=p.resolve(path.parent,inventory['benchmark']);plan=p.read_json(folder/'plan.json')
    checked=core_predecessor(plan['after_suite'],predecessor,p.sha(predecessor_path),registry,files)
    return dict(budget=budget,budget_path=budget_path,root=root,files=files,predecessor_checked=checked,
        specified_utc=workflow['specified_utc'],scope='core_v1',labels=labels,
        repetition_interpretation='One planned repetition per state/worker/workflow; descriptive observed costs only. Attempts after interruption are not independent repetitions; no between-repeat variation or confidence interval is estimated.')
