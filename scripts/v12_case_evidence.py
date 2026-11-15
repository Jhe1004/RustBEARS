"""Read-only verification of processed empirical workflows, including failures.

Counters in results.json are insufficient. Model identities, bootstrap/tree IDs,
map metadata and complete cladogenesis records must agree with source artifacts.
This module never fits, simulates, retries, or edits a case study.
"""
from collections import Counter
import csv
from datetime import datetime
import json
import math
from pathlib import Path

if __package__:
    from .plot_v12_calibration import checkpoint,read_json,sha
else:
    from plot_v12_calibration import checkpoint,read_json,sha

MODELS=['DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J']
CONTRACTS=['tip_survival','legacy_unconditioned']
PROCESSED_CASE_STATUSES={'complete','incomplete_profiles','complete_with_tree_failures'}
E4=('tree_id','tree_topology_hash','model_id','map_id','stratum_weight','precision_mode','likelihood_conditioning')


def _integer(value):
    number=float(value)
    if not math.isfinite(number) or number!=int(number):raise ValueError('invalid integer identifier')
    return int(number)


def verify_case_content(folder,manifest,report):
    """Check numerical artifacts; the formal/budget gate is in verify_case."""
    folder=Path(folder);files={};statistics={}
    def track(path):
        path=Path(path)
        if str(path) not in files:files[str(path)]=sha(path)
        return path
    def get(path):return read_json(track(path))
    def csv_rows(path):
        with track(path).open(newline='') as handle:
            yield from csv.DictReader(handle)
    def checked_store(path,n):
        track(path);rows,header=checkpoint(path,n)
        # Check every stored array shard without loading arrays into memory.
        with path.open() as handle:
            next(handle)
            for line in handle:
                item=json.loads(line)
                if item.get('shard') is not None:
                    name=item['shard']
                    if Path(name).name!=name:raise ValueError('checkpoint shard escapes its directory')
                    shard=path.with_name(path.name+'.shards')/name
                    track(shard)
                    if files[str(shard)]!=item['shard_sha256']:raise ValueError('case checkpoint shard checksum mismatch')
        return rows,header
    def model_definition(model,name,contract):
        if model['model']['name']!=name or model['config']['likelihood_conditioning']!=contract or model['config']['precision']!='strict_f64':
            raise ValueError('checkpoint model identity, contract or precision differs from case')
        bounds={k:[v['min'],v['max']] for k,v in model['model']['params'].items() if v['type']=='free'}
        if bounds!=manifest['model_parameter_bounds'][name]:raise ValueError('checkpoint physical parameter domain differs from reference case')
    def mapping(directory,prefix,expected,summary,tree_ids=None):
        if expected==0:
            if summary.get('status')!='no_complete_tree_maps':raise ValueError('zero complete trees require explicit unavailable mapping')
            return dict(maps=0,cladogenetic_events=0,anagenetic_events=0,status='unavailable')
        if summary.get('status')!='complete' or summary.get('maps')!=expected:raise ValueError('map summary count differs from planned count')
        observed=get(directory/(prefix+'bsm_summary.json'))
        if observed!=summary:raise ValueError('mapping summary differs from final report')
        metadata=list(csv_rows(directory/(prefix+'map_metadata.csv')))
        by_map={_integer(r['map']):r for r in metadata}
        tuples={tuple(r[k] for k in ('tree_id','model_id','map_id')) for r in metadata}
        if len(metadata)!=expected or len(tuples)!=expected or set(by_map)!=set(range(1,expected+1)):
            raise ValueError('missing or duplicated map ID/stratum tuple')
        for row in metadata:
            if any(row.get(k) in ('',None) for k in E4):raise ValueError('map lacks required event metadata')
            if row['likelihood_conditioning']!=directory.name or row['precision_mode']!='strict_f64':raise ValueError('map contract/precision differs from its case')
            if _integer(row['model_id']) not in range(1,7):raise ValueError('map has unknown model ID')
        if tree_ids is not None:
            counts=Counter(_integer(r['tree_id']) for r in metadata)
            if counts!={i:manifest['bsm_maps_per_tree'] for i in tree_ids}:raise ValueError('mapping discarded or invented a completed tree')
        clado=Counter();seen_nodes={};ana_count=0
        for kind in ('cladogenetic','anagenetic'):
            for row in csv_rows(directory/(prefix+kind+'_events.csv')):
                rid=_integer(row['map'])
                if rid not in by_map or any(row.get(k)!=by_map[rid][k] for k in E4):raise ValueError('event and map metadata disagree')
                if kind=='cladogenetic':
                    node=_integer(row['node'])-manifest['tip_count']-1
                    if not 0<=node<manifest['tip_count']-1:raise ValueError('cladogenesis node is not an internal node')
                    bit=1<<node
                    if seen_nodes.get(rid,0)&bit:raise ValueError('duplicated internal node/map event')
                    seen_nodes[rid]=seen_nodes.get(rid,0)|bit;clado[rid]+=1
                else:ana_count+=1
        if clado!={i:manifest['tip_count']-1 for i in by_map}:raise ValueError('missing cladogenetic event at an internal node/map')
        # Retain diagnostic tables even when their precision target is unmet.
        for suffix in ('event_precision','route_precision','route_comparisons','event_through_time'):
            track(directory/(prefix+suffix+'.csv'))
        return dict(maps=expected,cladogenetic_events=sum(clado.values()),anagenetic_events=ana_count,status='complete')

    if manifest.get('models')!=MODELS or manifest.get('likelihood_contracts')!=CONTRACTS or set(report.get('contracts',{}))!=set(CONTRACTS):
        raise ValueError('case omits a planned contract/model family')
    profiles_complete=True;trees_complete=True
    for contract in CONTRACTS:
        directory=folder/contract;summary=report['contracts'][contract]
        initial=get(directory/'initial_fitted_parameters.json');final=get(directory/'fitted_parameters.json')
        for fits in (initial,final):
            if [r['model'] for r in fits]!=MODELS or any(r.get('converged') is not True for r in fits):raise ValueError('case lacks six converged reference fits')
            for fit in fits:
                if fit['parameter_bounds']!=manifest['model_parameter_bounds'][fit['model']]:raise ValueError('physical parameter domain changed within case')
                if not math.isfinite(fit['lnL']):raise ValueError('nonfinite reference fit likelihood')
        for i,row in enumerate(initial,1):
            if get(directory/f'initial_fit_{i}.json')!=row:raise ValueError('initial fit ledger mismatch')
        weights=list(csv_rows(directory/'models.csv'))
        if [r['model'] for r in weights]!=MODELS or any(r['likelihood_conditioning']!=contract or r['precision_mode']!='strict_f64' for r in weights):
            raise ValueError('model-average table omits or mixes planned models/contracts')
        probabilities=[float(r['weight']) for r in weights]
        if any(not math.isfinite(v) or not 0<=v<=1 for v in probabilities) or abs(sum(probabilities)-1)>1e-10:
            raise ValueError('invalid model-average weights')
        track(directory/'ancestral_arrays.npz')
        profiles=get(directory/'profile_status.json')
        if profiles!=summary['profiles'] or [r['model'] for r in profiles]!=MODELS:raise ValueError('case profile inventory differs from report')
        for i,row in enumerate(profiles,1):
            if row['status']=='complete':
                intervals=list(csv_rows(directory/f'profile_{i}.csv'))
                expected=set(final[i-1]['parameter_bounds'])
                if len(intervals)!=len(expected) or {r['param'] for r in intervals}!=expected:raise ValueError('completed profile omits a free parameter')
                for parameter in expected:track(directory/f'profile_{i}_{parameter}_curve.csv')
            elif row['status']=='numerically_uncertified':profiles_complete=False
            else:raise ValueError('unprocessed profile status')
        if len(summary['bootstrap'])!=3:raise ValueError('missing bootstrap model comparison')
        bootstrap_failures=[]
        for pair in range(1,4):
            rows,header=checked_store(directory/f'bootstrap_{pair}.jsonl',manifest['n_boot'])
            null,alt=MODELS[2*(pair-1):2*pair];cell=summary['bootstrap'][pair-1]
            if header['n_boot']!=manifest['n_boot'] or header['n_alt_sim']!=0 or header['seed']!=manifest['seed']+pair-1:
                raise ValueError('bootstrap checkpoint configuration differs from case')
            for key,name in (('null',null),('alt',alt)):
                model_definition(header[key],name,contract)
            if header['simulation_conditioning']!='tip_survival' or cell['null']!=null or cell['alt']!=alt:
                raise ValueError('bootstrap generating contract or model pair mismatch')
            table=list(csv_rows(directory/f'bootstrap_{pair}.csv'))
            if sorted(_integer(r['rep_id']) for r in table)!=list(range(1,manifest['n_boot']+1)):raise ValueError('bootstrap CSV dropped a planned replicate')
            by_id={r['rep_id']:r for r in rows}
            if set(by_id)!=set(range(1,manifest['n_boot']+1)) or any(r['phase']!='null' for r in rows):raise ValueError('bootstrap checkpoint replicate identity mismatch')
            if any(r['status']!=by_id[_integer(r['rep_id'])]['status'] for r in table):raise ValueError('bootstrap CSV and checkpoint statuses differ')
            failures=sum(r['status']!='converged' for r in rows)
            if failures!=cell['n_failed']:raise ValueError('bootstrap summary erased retained failures')
            bootstrap_failures.append(failures)
        local=dict(bootstrap_failed_counts=bootstrap_failures,
            profiles_uncertified=sum(r['status']=='numerically_uncertified' for r in profiles),
            single_tree=mapping(directory,'',manifest['bsm_maps'],summary['single_tree_mapping']))
        if manifest['original_posterior_required']:
            rows,header=checked_store(directory/'multitree.jsonl',manifest['max_trees'])
            if (header['bsm_maps_per_tree']!=manifest['bsm_maps_per_tree'] or header['seed']!=manifest['seed']
                or len(header['trees'])!=manifest['max_trees'] or header['max_trees']!=manifest['max_trees'] or header['match_mode']!='monophyletic'):
                raise ValueError('multitree checkpoint configuration differs from case')
            if len(header['models'])!=6:raise ValueError('multitree checkpoint lacks six model definitions')
            for model,name in zip(header['models'],MODELS):model_definition(model,name,contract)
            complete=[]
            for row in rows:
                models=row['models']
                if [r['model'] for r in models]!=MODELS or [r['model_id'] for r in models]!=list(range(1,7)):
                    raise ValueError('tree record omits a candidate model')
                if any(r['tree_id']!=row['id'] or r['likelihood_conditioning']!=contract or r['fit_status'] not in ('converged','failed') for r in models):
                    raise ValueError('tree/model identity, fit status or contract mismatch')
                expected='complete' if all(r['fit_status']=='converged' for r in models) else 'incomplete_model_set'
                if row['tree_status']!=expected:raise ValueError('tree status hides incomplete model set')
                if expected=='complete':complete.append(row['id'])
            multi=summary['multitree'];n=len(complete);trees_complete &= n==manifest['max_trees']
            if (multi['n_planned'],multi['n_complete'],multi['n_incomplete'])!=(manifest['max_trees'],n,manifest['max_trees']-n):
                raise ValueError('multitree report differs from per-tree evidence')
            table=list(csv_rows(directory/'per_tree_models.csv'))
            if sorted((_integer(r['tree_id']),_integer(r['model_id'])) for r in table)!=[(i,j) for i in range(1,manifest['max_trees']+1) for j in range(1,7)]:
                raise ValueError('multitree CSV omitted failed model/tree rows')
            for name in ('multitree_probability_bounds.csv','tree_parameter_sensitivity.csv','tree_model_winner_frequency.csv'):
                track(directory/name)
            local['trees_complete']=n;local['trees_incomplete']=manifest['max_trees']-n
            local['joint']=mapping(directory,'joint_',n*manifest['bsm_maps_per_tree'],multi['joint_mapping'],complete)
        statistics[contract]=local
    expected='complete' if profiles_complete and trees_complete else 'incomplete_profiles' if not profiles_complete else 'complete_with_tree_failures'
    if report.get('status')!=expected:raise ValueError('case terminal status disagrees with retained profile/tree failures')
    for path in folder.rglob('*'):
        if path.suffix in ('.json','.csv','.jsonl','.npz') and 'runtime_cache' not in path.parts:track(path)
    return dict(files=files,status=expected,contracts=statistics)


def verify_case(folder,job,budget,budget_hash):
    folder=Path(folder);manifest=read_json(folder/'manifest.json');report=read_json(folder/'results.json')
    if manifest.get('mode')!='formal' or manifest.get('formal') is not True or report.get('formal') is not True or report.get('status') not in PROCESSED_CASE_STATUSES:
        raise ValueError('case is not a processed formal workflow')
    metadata=manifest['locked_budget']
    if metadata['sha256']!=budget_hash or metadata['payload']!=budget:raise ValueError('case budget mismatch')
    if any(manifest.get(k)!=v for k,v in job.items() if k!='analysis'):raise ValueError('case settings differ from planned job')
    if manifest['original_posterior_required']!=(job['analysis']=='d2'):raise ValueError('case posterior source differs from job')
    if manifest['workflow_script_sha256']!=budget['source_sha256'].get('scripts/run_uncertainty_case_study.py'):
        raise ValueError('case workflow source differs from frozen budget')
    started=read_json(folder/'execution_started.json')
    first=datetime.fromisoformat(started['first_started_utc'].replace('Z','+00:00'))
    locked=datetime.fromisoformat(metadata['locked_utc'].replace('Z','+00:00'))
    if started['budget_sha256']!=budget_hash or started.get('gate2_gate3_verified') is not True or first.tzinfo is None or first<locked:
        raise ValueError('case execution predates its locked budget or gate')
    checked=verify_case_content(folder,manifest,report)
    checked['files'].update({str(folder/name):sha(folder/name) for name in ('manifest.json','results.json','execution_started.json')})
    return checked
