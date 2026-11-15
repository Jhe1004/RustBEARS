#!/usr/bin/env python3
"""Plan/check Figure 2; only verified formal inventories can export tables or render.

No smoke promotion, mock rendering, seed selection or missing-value imputation.
The default writes a contract and checks provenance without opening a graphics
backend. The final-size visual review remains a human gate after actual export.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
ANALYSES=('v1','v2b','v2c','v3','v4')
V1_SCENARIOS=('interior','e_zero','j_zero','double_zero')
V2_SCENARIOS=('null_e0','null_e02','alt_j001','alt_j003','alt_j008')
METHODS=('asymptotic_chi2','onesided_zero_reference','wald_local','bootstrap_plugin')
ROUTES=('A','C','C_to_A')
CONTRACT=dict(schema='rustbears.figure2.contract.v1',status='awaiting_formal_data',backend='python',
    question='How do uncertainty statements behave under the declared boundary, contract, tree and stopping conditions?',
    result_claim='Pending formal data; no universal coverage or anytime-validity claim.',
    main_panels=dict(a=dict(analysis='v1',parameter='j',areas=4,tips=[40,150],scenarios=list(V1_SCENARIOS),methods=list(METHODS)),
        b=dict(analysis='v2b',scenarios=list(V2_SCENARIOS),routes=list(ROUTES),finite_b_extrapolation=False),
        c=dict(analysis='v3',support_bins=['<0.7','0.7-0.95','>0.95'],contrast='dataset-paired Brier: multiple trees minus MCC'),
        d=dict(analysis='v4',sections=['fixed','adaptive'],events=['d','e'])),
    supplemental='All registered jobs, parameters, methods, predictions, finite-B/surrogate summaries, tree-MCSE and route/stopping results; raw numerical records retained.',
    missingness='NA remains NA. Failures remain planned outcomes. Failure bounds are distinct from confidence intervals.',
    rate_display='Available-outcome estimate; error bars use exact-binomial planned-outcome CP envelopes. Thick segments show deterministic failure bounds.',
    dimensions_mm=[183,230],minimum_glyph_pt=8,primary_font_pt=8,
    alignment_tolerance_pt=1.5,qa='Source, alignment, PDF glyph and collision audits; final-size human inspection remains pending.')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reject_constant(value):raise ValueError('non-standard JSON number: '+value)


def read_json(path):return json.loads(Path(path).read_text(),parse_constant=reject_constant)


def write_json(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def resolve(base,value):
    value=Path(value)
    return value.resolve() if value.is_absolute() else (Path(base)/value).resolve()


def decode(value):
    """Read typed checkpoint manifests without reconstructing executable objects."""
    if isinstance(value,list):return [decode(x) for x in value]
    if not isinstance(value,dict):return value
    if '__mapping__' in value:return {decode(k):decode(v) for k,v in value['__mapping__']}
    if '__sequence__' in value:
        result=[decode(x) for x in value['items']]
        return tuple(result) if value['__sequence__']=='tuple' else result
    if '__dataclass__' in value:return {k:decode(v) for k,v in value['fields'].items()}
    if '__float64__' in value:return struct.unpack('d',bytes.fromhex(value['__float64__']))[0]
    if '__nonfinite_float__' in value:return value['__nonfinite_float__']
    return {k:decode(v) for k,v in value.items()}


def job_key(job):
    analysis=job['analysis']
    if analysis in ('v1','v2b','v2c'):
        return (analysis,job['scenario'],job.get('n_tips'),job.get('areas'),bool(job.get('finite_b_extrapolation',False)))
    return (analysis,)


def budget_info(inventory,base):
    if inventory.get('schema')!='rustbears.figure2.inputs.v1' or inventory.get('formal') is not True:
        raise ValueError('Figure 2 requires a formal=true inventory; smoke/pilot promotion is forbidden')
    item=inventory['budget'];path=resolve(base,item['path'])
    if sha(path)!=item['sha256']:raise ValueError('budget SHA differs from the inventory')
    budget=read_json(path)
    if budget.get('schema')!='rustbears.formal_budget.v1' or budget.get('status')!='locked':
        raise ValueError('only a prospectively locked formal budget is eligible')
    locked=datetime.fromisoformat(budget['locked_utc'].replace('Z','+00:00'))
    if locked.tzinfo is None or locked>datetime.now(timezone.utc):raise ValueError('invalid prospective lock time')
    gate=budget['gate_evidence'];gate_path=resolve(path.parent,gate['path'])
    if sha(gate_path)!=gate['sha256']:raise ValueError('gate evidence SHA mismatch')
    gates=read_json(gate_path)
    if gates.get('gate2_passed') is not True or gates.get('gate3_passed') is not True:
        raise ValueError('actual Gate 2 and Gate 3 passes are required')
    if 'design_evidence' in budget:
        item=budget['design_evidence'];design_path=resolve(path.parent,item['path'])
        if sha(design_path)!=item['sha256'] or read_json(design_path)['jobs']!=budget['jobs']:
            raise ValueError('final settings differ from the prospectively frozen design')
    jobs=[j for j in budget['jobs'] if j['analysis'] in ANALYSES]
    if set(j['analysis'] for j in jobs)!=set(ANALYSES):raise ValueError('all V1/V2b/V2c/V3/V4 analyses must be registered')
    if len({job_key(j) for j in jobs})!=len(jobs):raise ValueError('ambiguous registered job keys')
    return budget,path,jobs


def checkpoint(path,planned,*,expected_fingerprint=None,expected_manifest=None):
    raw=Path(path).read_bytes()
    if not raw.endswith(b'\n'):raise ValueError('uncommitted/truncated checkpoint tail')
    lines=[json.loads(line,parse_constant=reject_constant) for line in raw.splitlines()]
    if not lines or lines[0].get('schema')!='rustbears.checkpoint.v1':raise ValueError('invalid checkpoint header')
    if expected_fingerprint is not None and lines[0]['task_fingerprint']!=expected_fingerprint:
        raise ValueError('checkpoint task fingerprint mismatch')
    actual=hashlib.sha256(json.dumps(lines[0]['manifest'],sort_keys=True,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    if actual!=lines[0]['task_fingerprint']:raise ValueError('checkpoint manifest digest does not match its contents')
    if expected_manifest is not None and lines[0]['manifest']!=expected_manifest:
        raise ValueError('checkpoint differs from the preregistered complete configuration')
    records=[decode(row['record']) for row in lines[1:]]
    if sorted(row['id'] for row in records)!=list(range(1,planned+1)):
        raise ValueError('every planned outcome must be processed exactly once, including failures')
    return records,decode(lines[0]['manifest'])


def verify_lock(metadata,budget_hash):
    if not isinstance(metadata,dict) or metadata.get('sha256')!=budget_hash:
        raise ValueError('run is not bound to this locked budget; no smoke promotion')
    if metadata.get('payload',{}).get('status')!='locked':raise ValueError('run records an unlocked budget')


def verify_workflow_source(manifest,analysis,budget):
    if analysis in ('v1','v2b','v2c'):
        workflow='scripts/run_uncertainty_validation.py';recorded=manifest.get('script_sha256')
    elif analysis=='v4':
        workflow='scripts/validate_v12_v4_protocol.py';recorded=manifest['source_sha256'].get(workflow)
    else:
        return 'V3 coordinator source recorded through budget; no separate execution-time workflow digest'
    if not recorded or recorded!=budget['source_sha256'].get(workflow):
        raise ValueError('workflow source differs from the locked budget or lacks its digest')
    return 'recorded workflow digest matches locked budget'


def load_run(item,base,budget_hash):
    folder=resolve(base,item['path']);analysis=item['analysis'];files=[]
    def get(name):
        path=folder/name;value=read_json(path);files.append(path);return value
    if analysis in ('v1','v2b','v2c'):
        prereg=get('preregistered_manifest.json');manifest=decode(prereg['manifest']);summary=get('summary.json')
        if manifest.get('mode')!='formal' or summary.get('mode')!='formal':raise ValueError('smoke/pilot validation runs cannot be rendered')
        verify_lock(manifest.get('locked_budget'),budget_hash)
        if summary.get('locked_budget_sha256')!=budget_hash:raise ValueError('summary budget mismatch')
        started=get('execution_started.json')
        if started.get('budget_sha256')!=budget_hash or started.get('gate2_gate3_verified') is not True:raise ValueError('formal execution-start evidence missing')
        first=datetime.fromisoformat(started['first_started_utc'].replace('Z','+00:00'))
        locked=datetime.fromisoformat(manifest['locked_budget']['locked_utc'].replace('Z','+00:00'))
        if first.tzinfo is None or locked.tzinfo is None or first<locked:raise ValueError('budget was not locked before execution')
        records,header=checkpoint(folder/'outer_records.jsonl',manifest['outer_replicates'],expected_fingerprint=prereg['task_fingerprint'],expected_manifest=prereg['manifest'])
        files.append(folder/'outer_records.jsonl')
        if summary.get('n_records')!=len(records) or summary.get('n_planned')!=manifest['outer_replicates']:
            raise ValueError('summary is incomplete or has a different planned denominator')
        if header.get('locked_budget',{}).get('sha256')!=budget_hash or job_key(header)!=job_key(manifest):raise ValueError('checkpoint configuration changed')
        key=job_key(manifest)
    elif analysis=='v3':
        manifest=get('execution_manifest.json');summary=get('results.json')
        if summary.get('formal') is not True or summary.get('status')!='all_planned_datasets_processed':raise ValueError('V3 is not a processed formal study')
        verify_lock(manifest.get('budget'),budget_hash)
        if summary.get('budget_sha256')!=budget_hash:raise ValueError('V3 summary budget mismatch')
        records=[dict(r) for r in summary['outcomes']];planned=manifest['settings']['replicates']
        if sorted(r['dataset'] for r in records)!=list(range(1,planned+1)):raise ValueError('V3 discarded a planned dataset')
        if len(records)!=summary['n_planned']:raise ValueError('V3 planned denominator mismatch')
        sys.path.insert(0,str(ROOT))
        from scripts.run_v12_v3_study import verify_prepared_study,verify_completed_dataset
        verify_prepared_study(folder,manifest['settings'],budget_hash)
        for row in records:
            ledger=get(f'dataset_{row["dataset"]:03d}/study_outcome.json')
            if ledger!=row:raise ValueError('V3 outcome differs from its retained ledger')
            if row['status']=='complete':
                verify_completed_dataset(folder/f'dataset_{row["dataset"]:03d}',row)
                row['score']=get(f'dataset_{row["dataset"]:03d}/v3_scores.json')
            elif row['status']!='failed_retained':raise ValueError('V3 has an unprocessed outcome')
        key=('v3',)
    elif analysis=='v4':
        raw=get('manifest.json');manifest=decode(raw['manifest']);summary=get('results.json')
        if manifest.get('mode')!='formal' or summary.get('formal') is not True:raise ValueError('V4 smoke/subset is not a complete formal protocol')
        verify_lock(manifest.get('locked_budget'),budget_hash)
        if summary.get('status')!='completed' or summary.get('source_unchanged') is not True:raise ValueError('V4 run did not finish with stable sources')
        settings=manifest['execution_settings'];planned=sum(s['replicates'] for s in settings.values())
        records,header=checkpoint(folder/'records.jsonl',planned,expected_fingerprint=raw['task_fingerprint'],expected_manifest=raw['manifest']);files.append(folder/'records.jsonl')
        if set(manifest['sections'])!=set(settings):raise ValueError('V4 omits registered sections')
        for section,config in settings.items():
            if sorted(r['replicate'] for r in records if r['section']==section)!=list(range(1,config['replicates']+1)):
                raise ValueError('V4 section has missing or duplicated outcomes')
        key=('v4',)
    else:raise ValueError('unknown analysis in input inventory')
    return dict(analysis=analysis,key=key,path=folder,manifest=manifest,summary=summary,records=records,
                files={str(p):sha(p) for p in files})


def inspect_inputs(path):
    path=Path(path).resolve();inventory=read_json(path);budget,budget_path,jobs=budget_info(inventory,path.parent)
    budget_hash=sha(budget_path);runs=[load_run(item,path.parent,budget_hash) for item in inventory['runs']]
    expected={job_key(j):j for j in jobs}
    if len(runs)!=len(expected) or {r['key'] for r in runs}!=set(expected):raise ValueError('every required registered job must be present exactly once')
    for run in runs:
        job=expected[run['key']];run['job']=job
        metadata=run['manifest'].get('budget') if run['analysis']=='v3' else run['manifest'].get('locked_budget')
        if metadata['payload']!=budget:raise ValueError('stored budget payload differs from the locked budget bytes')
        run['workflow_source_check']=verify_workflow_source(run['manifest'],run['analysis'],budget)
        settings=run['manifest']['settings'] if run['analysis']=='v3' else run['manifest']['execution_settings'] if run['analysis']=='v4' else run['manifest']
        if run['analysis']=='v4':
            if settings!=job['settings']:raise ValueError('V4 settings differ from the budget')
        elif run['analysis']=='v3':
            if settings!=job:raise ValueError('V3 settings differ from the budget')
        else:
            for key in ('outer_replicates','inner_bootstrap','plugin_bootstrap','surrogate_alt_replicates','seed'):
                if settings[key]!=job[key]:raise ValueError('validation settings differ from the budget: '+key)
    main_a={(j['n_tips'],j['scenario']) for j in jobs if j['analysis']=='v1' and j['areas']==4}
    if main_a!={(n,s) for n in (40,150) for s in V1_SCENARIOS}:raise ValueError('predeclared panel a requires all eight 4-area cells')
    if {j['scenario'] for j in jobs if j['analysis']=='v2b' and not j['finite_b_extrapolation']}!=set(V2_SCENARIOS):raise ValueError('panel b requires all five scenarios')
    supporting_files=[resolve(budget_path.parent,budget['gate_evidence']['path'])]
    if 'design_evidence' in budget:supporting_files.append(resolve(budget_path.parent,budget['design_evidence']['path']))
    return runs,dict(status='eligible_formal_inputs',formal=True,inventory_sha256=sha(path),budget_sha256=budget_hash,
        jobs=len(runs),planned_outcomes=sum(len(r['records']) for r in runs),
        workflow_source_checks={str(r['path']):r['workflow_source_check'] for r in runs},
        source_files={str(path):sha(path),str(budget_path):budget_hash,**{str(p):sha(p) for p in supporting_files},**{p:h for r in runs for p,h in r['files'].items()}},
        exclusions='No registered analysis removed; main-panel subset fixed by the figure contract; all other results retained in supplementary tables.')


def rate_row(values,planned,*,possible=None):
    """Keep available rates, deterministic failure bounds and CP envelopes distinct."""
    from scipy.stats import beta
    values=list(values)
    if any(v is not None and not isinstance(v,bool) for v in values):raise ValueError('outcomes must be booleans or explicit null')
    valid=[v for v in values if v is not None];hit=sum(valid);n=len(valid);missing=planned-n
    if not 0<=hit<=n<=planned:raise ValueError('invalid planned/observed outcome counts')
    lo,hi=(hit,hit+missing) if possible is None else possible
    if not 0<=lo<=hi<=planned:raise ValueError('invalid missing-outcome counts')
    cp=None if n==0 else [0. if hit==0 else float(beta.ppf(.025,hit,n-hit+1)),1. if hit==n else float(beta.ppf(.975,hit+1,n-hit))]
    envelope=None if planned==0 else [0. if lo==0 else float(beta.ppf(.025,lo,planned-lo+1)),1. if hi==planned else float(beta.ppf(.975,hi+1,planned-hi))]
    return dict(n_planned=planned,n_observed=n,n_missing=missing,n_true=hit,
        estimate=None if n==0 else hit/n,cp95_available=cp,cp95_planned_envelope=envelope,
        failure_bounds=None if planned==0 else [lo/planned,hi/planned],
        status='unavailable' if n==0 else 'complete' if missing==0 and lo==hi else 'missing_or_failure_bounded')


def main_tables(runs):
    tables={name:[] for name in ('a','b','c','d')}
    for run in runs:
        analysis=run['analysis'];job=run['job'];rows=run['records'];summary=run['summary']
        common=dict(source_run=str(run['path']),source_analysis=analysis)
        if analysis=='v1' and job['areas']==4:
            for method in METHODS:
                key=('bootstrap_plugin' if method=='bootstrap_plugin' else 'asymptotic_chi2')+'_j'
                field={'onesided_zero_reference':'one_sided_reference_covers','wald_local':'wald_covers'}.get(method,'covered')
                values=[r['outputs'].get(key,{}).get(field) for r in rows]
                tables['a'].append(dict(common,n_tips=job['n_tips'],areas=job['areas'],scenario=job['scenario'],method=method,
                    j_truth=run['manifest']['truth']['j'],**rate_row(values,job['outer_replicates'])))
        elif analysis=='v2b' and not job['finite_b_extrapolation']:
            for route in ROUTES:
                available=[r['outputs'][route] for r in rows if route in r['outputs']]
                lower=sum(r['p_failure_bounds'][1]<=.05 for r in available)
                upper=sum(r['p_failure_bounds'][0]<=.05 for r in available)+job['outer_replicates']-len(available)
                tables['b'].append(dict(common,scenario=job['scenario'],route=route,inner_bootstrap=job['inner_bootstrap'],
                    n_outer_with_inner_failure=sum(r['n_failed']>0 for r in available),
                    **rate_row([r['reject'] for r in available],job['outer_replicates'],possible=(lower,upper))))
        elif analysis=='v3':
            for support in ('<0.7','0.7-0.95','>0.95'):
                cell=summary['support_bins'][support];n=cell['n_datasets_with_scored_clades']
                tables['c'].append(dict(common,support_bin=support,n_observed=n,n_planned=summary['n_planned'],
                    n_failed_datasets=summary['n_missing'],n_clades=cell['n_clades'],
                    estimate=cell['paired_dataset_mean_brier_difference'] if n else None,
                    ci95=cell['paired_dataset_bootstrap_ci95'] if n else None,
                    status='unavailable' if n==0 else 'single_dataset_no_interval' if n==1 else 'complete_case_dataset_pair'))
        elif analysis=='v4':
            for section in ('fixed','adaptive'):
                for event in ('d','e'):
                    cell=summary['summaries'][section][event]
                    available=[r for r in rows if r['section']==section and event in r.get('estimates',{})]
                    truth=cell['truth'];values=[r['estimates'][event] for r in available]
                    hits=[v['mean']-1.96*v['se']<=truth<=v['mean']+1.96*v['se']
                          if all(isinstance(v[k],(int,float)) and not isinstance(v[k],bool) and math.isfinite(v[k]) for k in ('mean','se')) and v['se']>=0
                          else None for v in values]
                    tables['d'].append(dict(common,section=section,event=event,truth=truth,
                        **rate_row(hits,job['settings'][section]['replicates'])))
    if [len(tables[x]) for x in ('a','b','c','d')]!=[32,15,3,4]:raise ValueError('predeclared main-panel inventory is incomplete')
    return tables


def leaf_rows(value,prefix=''):
    if isinstance(value,dict):
        for key,child in value.items():yield from leaf_rows(child,prefix+'/'+str(key))
    elif isinstance(value,(list,tuple)):
        for i,child in enumerate(value):yield from leaf_rows(child,prefix+'/'+str(i))
    else:yield dict(metric_path=prefix,value=value,value_type=type(value).__name__)


def write_csv(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields or ['status']);writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,ensure_ascii=False,allow_nan=False) if isinstance(v,(dict,list,tuple)) else v for k,v in row.items()})


def export_tables(runs,tables,out):
    folder=out/'source_data';folder.mkdir(exist_ok=True)
    for panel,rows in tables.items():write_csv(folder/f'figure2_{panel}.csv',rows)
    leaves=[];intervals=[];lrt=[];predictions=[];weights=[];clades=[];datasets=[];raw=[];indices=[]
    for run in runs:
        source=dict(analysis=run['analysis'],source_run=str(run['path']),job=json.dumps(run['job'],sort_keys=True))
        leaves.extend(dict(source,**leaf) for leaf in leaf_rows(run['summary']))
        for row in run['records']:
            identifier=row.get('rep_id',row.get('dataset',row.get('id')))
            raw.append(dict(source,record_id=identifier,record=row))
            if run['analysis']=='v1':
                for key,value in row['outputs'].items():intervals.append(dict(source,replicate=identifier,interval=key,**value))
            elif run['analysis']=='v2b':
                for route,value in row['outputs'].items():lrt.append(dict(source,replicate=identifier,route=route,**value))
            elif run['analysis']=='v2c':
                for name,value in row['outputs'].items():
                    if name=='weights':weights.extend(dict(source,replicate=identifier,**v) for v in value)
                    else:predictions.append(dict(source,replicate=identifier,method=name,**value))
            elif run['analysis']=='v3' and row['status']=='complete':
                for clade in row['score']['paired_clades']:clades.append(dict(source,dataset=identifier,**clade))
        if run['analysis']=='v3':
            for support,cell in run['summary']['support_bins'].items():
                datasets.extend(dict(source,support_bin=support,**r) for r in cell['dataset_scores'])
        # The checkpoint index binds all inner/bootstrap/profile numerical records
        # and their recorded array-shard checksums, without duplicating huge arrays.
        for path in sorted(run['path'].rglob('*.jsonl')):
            indices.append(dict(analysis=run['analysis'],path=str(path),sha256=sha(path),bytes=path.stat().st_size))
    for name,rows in [('all_summary_values',leaves),('v1_all_parameter_intervals',intervals),('v2b_all_outer_route_results',lrt),
                      ('v2c_all_predictions',predictions),('v2c_all_model_weights',weights),('v3_all_shared_clade_scores',clades),('v3_all_dataset_scores',datasets)]:
        write_csv(folder/(name+'.csv'),rows)
    with (folder/'all_outer_and_mapping_records.jsonl').open('w') as handle:
        for row in raw:handle.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n')
    with (folder/'all_run_summaries.jsonl').open('w') as handle:
        for run in runs:handle.write(json.dumps(dict(source_run=str(run['path']),summary=run['summary']),ensure_ascii=False,allow_nan=False)+'\n')
    write_json(folder/'checkpoint_index.json',indices)
    write_json(folder/'source_data_scope.json',dict(main_selection=CONTRACT['main_panels'],
        all_registered_summaries_and_outcomes_retained=True,unavailable_values='empty CSV field or JSON null; never zero-imputed',
        checkpoint_arrays='Original checkpoints and array files remain untouched; their index and stored shard digests are retained, not flattened into CSV.',
        interval_definitions=CONTRACT['rate_display']))
    return {p.relative_to(out).as_posix():sha(p) for p in folder.iterdir() if p.is_file()}


def audit_command(command,log):
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
    Path(log).write_text(result.stdout)
    if result.returncode:raise RuntimeError('Figure QA blocked export: '+Path(command[1]).name+'; see '+str(log))


def render(tables,out,qa_tools,*,eligibility=None):
    """Called only after every registered formal outcome passes provenance checks."""
    if not isinstance(eligibility,dict) or eligibility.get('status')!='eligible_formal_inputs' or eligibility.get('formal') is not True:
        raise ValueError('render requires verified formal input eligibility; no synthetic or smoke image')
    if not eligibility.get('source_files') or any(sha(path)!=digest for path,digest in eligibility['source_files'].items()):
        raise ValueError('formal sources changed after verification')
    qa_tools=Path(qa_tools).resolve()
    required=('validate_figure.py','audit_panel_alignment.py','audit_pdf_text.py','audit_figure_collisions.py')
    if any(not (qa_tools/name).is_file() for name in required):raise ValueError('complete QA tools directory is required before rendering')
    audit_command([sys.executable,str(qa_tools/'validate_figure.py'),__file__],out/'figure2.source-audit.log')
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sys.path.insert(0,str(qa_tools))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
        'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
        'legend.fontsize':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.right':False,
        'axes.spines.top':False,'axes.linewidth':.65,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(7.204724409,9.055118110),gridspec_kw={'height_ratios':[1.65,1]})
    fig.subplots_adjust(left=.20,right=.92,bottom=.10,top=.91,wspace=.70,hspace=.62)
    axs=dict(zip('abcd',axes.ravel()));headers=[]
    method_style={name:(color,marker,label) for name,color,marker,label in zip(METHODS,
        ['#365E7D','#A26C35','#6E6E6E','#57896D'],['o','^','s','D'],['Profile','Zero upper','Wald','Plug-in'])}
    route_style=dict(A=('#365E7D','o','A'),C=('#8B6E91','s','C'),C_to_A=('#57896D','D','C → A'))
    map_style=dict(fixed=('#365E7D','o','Fixed'),adaptive=('#A26C35','D','Adaptive'))
    def header(ax,key,styles=None):
        box=ax.get_position();h=fig.add_axes([box.x0,box.y1+.012,box.width,.045]);h.set_axis_off();headers.append(h)
        if styles:
            handles=[Line2D([],[],color=c,marker=m,linestyle='none',markersize=3.5,label=label) for c,m,label in styles]
            h.legend(handles=handles,loc='center',ncol=2 if len(styles)>3 else len(styles),frameon=False,
                     handletextpad=.4,columnspacing=.9,borderaxespad=0)
        ax.annotate(key+')',(0,1),xytext=(-24,15),textcoords='offset points',fontweight='bold',fontsize=8,ha='left',va='bottom')
    def point(ax,row,y,color,marker):
        value=row['estimate']
        if value is None:
            ax.text(1.02,y,'NA',transform=ax.get_yaxis_transform(),ha='left',va='center',fontsize=8,color=color,clip_on=False)
            return
        bounds=row.get('failure_bounds')
        if bounds is not None and bounds[0]!=bounds[1]:ax.plot(bounds,[y,y],color=color,lw=3,alpha=.25,solid_capstyle='butt')
        interval=row.get('cp95_planned_envelope',row.get('ci95'))
        if interval is not None:
            if not interval[0]<=value<=interval[1]:raise ValueError('estimate lies outside its declared interval')
            ax.errorbar(value,y,xerr=[[value-interval[0]],[interval[1]-value]],fmt=marker,
                        color=color,markersize=3.5,lw=.75,capsize=1.5)
        else:ax.plot(value,y,marker=marker,color=color,markersize=3.5,linestyle='none',fillstyle='none')
    def rate_axes(ax,labels):
        ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.set_yticks(range(len(labels)),labels)
        ax.set_ylim(len(labels)-.5,-.5);ax.grid(axis='x',color='#E7E7E7',lw=.5,zorder=0)
        ax.tick_params(axis='y',length=0,pad=5)
    titles=dict(interior='Interior',e_zero='e = 0',j_zero='j = 0',double_zero='e = j = 0')
    ax=axs['a'];labels=[]
    for i,(tips,scenario) in enumerate((n,s) for n in (40,150) for s in V1_SCENARIOS):
        local=[next(r for r in tables['a'] if r['n_tips']==tips and r['scenario']==scenario and r['method']==method) for method in METHODS]
        labels.append(f'{tips} tips: {titles[scenario]}\nN = '+ '/'.join(str(r['n_observed']) for r in local))
        for offset,row,method in zip([-.27,-.09,.09,.27],local,METHODS):point(ax,row,i+offset,*method_style[method][:2])
    rate_axes(ax,labels);ax.axvline(.95,color='#777777',linestyle='--',lw=.6);ax.set_xlabel('j interval coverage')
    header(ax,'a',list(method_style.values()))
    labels=[];ax=axs['b'];names=dict(null_e0='Null: e = 0',null_e02='Null: e = 0.02',alt_j001='j = 0.01',alt_j003='j = 0.03',alt_j008='j = 0.08')
    for i,scenario in enumerate(V2_SCENARIOS):
        local=[next(r for r in tables['b'] if r['scenario']==scenario and r['route']==route) for route in ROUTES]
        labels.append(names[scenario]+'\nN = '+'/'.join(str(r['n_observed']) for r in local))
        for offset,row,route in zip([-.20,0,.20],local,ROUTES):point(ax,row,i+offset,*route_style[route][:2])
    rate_axes(ax,labels);ax.axvline(.05,color='#777777',linestyle='--',lw=.6);ax.set_xlabel('Bootstrap rejection rate')
    header(ax,'b',list(route_style.values()))
    ax=axs['c'];labels=[];extent=.01
    for i,row in enumerate(tables['c']):
        labels.append(row['support_bin']+'\nN = '+str(row['n_observed']))
        point(ax,row,i,'#365E7D','o')
        if row['estimate'] is not None:extent=max(extent,abs(row['estimate']),*(abs(x) for x in (row['ci95'] or [])))
    ax.set_xlim(-1.12*extent,1.12*extent);ax.axvline(0,color='#777777',linestyle='--',lw=.6)
    ax.set_yticks(range(3),labels);ax.set_ylim(2.5,-.5);ax.tick_params(axis='y',length=0,pad=5)
    ax.set_xlabel('Brier difference:\nMultiple trees − MCC');ax.set_ylabel('Clade support',fontsize=8,labelpad=3)
    header(ax,'c')
    labels=[];ax=axs['d']
    for i,event in enumerate(('d','e')):
        local=[next(r for r in tables['d'] if r['event']==event and r['section']==section) for section in ('fixed','adaptive')]
        labels.append(('Dispersal (d)' if event=='d' else 'Extinction (e)')+'\nN = '+'/'.join(str(r['n_observed']) for r in local))
        for offset,row,section in zip([-.12,.12],local,('fixed','adaptive')):point(ax,row,i+offset,*map_style[section][:2])
    rate_axes(ax,labels);ax.axvline(.95,color='#777777',linestyle='--',lw=.6);ax.set_xlabel('Mapping interval coverage')
    header(ax,'d',list(map_style.values()))
    candidate=out/'qa_candidate';candidate.mkdir(exist_ok=True);stem=candidate/'figure2'
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,axes=list(axs.values()),panel_ids=list('abcd'),
        row_groups=[['a','b'],['c','d']],column_groups=[['a','c'],['b','d']],exclude_axes=headers,
        json_out=out/'figure2.alignment.json',overlay_svg=out/'figure2.alignment.svg',
        tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=True,strict=True)
    fig.savefig(stem.with_suffix('.pdf'),dpi=600)
    fig.savefig(stem.with_suffix('.svg'),dpi=600)
    fig.savefig(stem.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(stem.with_suffix('.png'),dpi=600)
    plt.close(fig)
    audit_command([sys.executable,str(qa_tools/'audit_pdf_text.py'),str(stem.with_suffix('.pdf')),'--min-pt','8','--json'],out/'figure2.glyph-audit.json')
    audit_command([sys.executable,str(qa_tools/'audit_figure_collisions.py'),str(stem.with_suffix('.pdf')),
                   '--json-out',str(out/'figure2.collision-audit.json'),'--overlay-pdf',str(out/'figure2.collision-overlay.pdf'),'--strict'],out/'figure2.collision-audit.log')
    for suffix in ('pdf','svg','tiff','png'):stem.with_suffix('.'+suffix).replace(out/('figure2.'+suffix))
    return dict(status='rendered_pending_final_size_visual_review',automatic_checks='passed',
                human_panel_review='pending',publication_ready=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--export-tables',action='store_true');p.add_argument('--render',action='store_true')
    p.add_argument('--qa-tools-dir',type=Path)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    write_json(out/'figure2_contract.json',CONTRACT)
    write_json(out/'input_inventory_template.json',dict(schema='rustbears.figure2.inputs.v1',formal=False,
        budget=dict(path='formal_budget.locked.json',sha256='PENDING'),runs=[],status='template_only_not_formal_data'))
    if args.inputs is None:
        audit=dict(status='awaiting_formal_inputs',formal=False,rendered=False)
        write_json(out/'input_audit.json',audit)
        if args.render or args.export_tables:p.error('formal input inventory is required; no figure or source table produced')
        print(json.dumps(audit));return
    try:runs,audit=inspect_inputs(args.inputs)
    except (ValueError,KeyError,FileNotFoundError,TypeError) as error:
        write_json(out/'input_audit.json',dict(status='ineligible_inputs',formal=False,rendered=False,error=str(error)))
        p.error(str(error))
    write_json(out/'input_audit.json',audit)
    if not args.render and not args.export_tables:
        print(json.dumps(dict(status='formal_inventory_checked_no_render',jobs=len(runs))));return
    tables=main_tables(runs)
    if not all(sha(path)==digest for path,digest in audit['source_files'].items()):raise ValueError('source inputs changed during inspection')
    artifacts=export_tables(runs,tables,out)
    state=dict(status='formal_source_tables_exported',formal=True,rendered=False,artifacts=artifacts,
        plotting_source_sha256=sha(__file__),source_inventory_sha256=audit['inventory_sha256'])
    if args.render:
        if args.qa_tools_dir is None:p.error('--qa-tools-dir is required for render; no graphic backend opened')
        try:state.update(render(tables,out,args.qa_tools_dir,eligibility=audit),rendered=True)
        except Exception as error:
            state.update(status='render_or_qa_blocked',rendered=False,error=type(error).__name__+': '+str(error))
            write_json(out/'figure2_status.json',state);raise
    write_json(out/'figure2_status.json',state)
    print(json.dumps(dict(status=state['status'],formal=True,rendered=state['rendered'])))


if __name__=='__main__':main()
