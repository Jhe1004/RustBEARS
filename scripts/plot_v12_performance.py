#!/usr/bin/env python3
"""Plan Figure 1; render formal benchmarks and source-bound API workflow costs.

The workflow-only option draws a methods diagram without inventing data panels.
Matched process ratios include startup and are never called kernel speedups.
"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
if __package__:
    from .v12_figure_utils import (audit_alignment, configure_plotting, copy_evidence,
                                  export_rows, save_figure, sha, write_json)
    from .v12_workflow_costs import WORKLOADS, LABELS, inspect_suite, export_costs
else:
    from v12_figure_utils import (audit_alignment, configure_plotting, copy_evidence,
                                  export_rows, save_figure, sha, write_json)
    from v12_workflow_costs import WORKLOADS, LABELS, inspect_suite, export_costs

CONTRACT = dict(figure='1', backend='python', final_size_mm=[180, 225], minimum_font_pt=8,
    question='How do explicit probability conventions connect analysis and measured computational cost?',
    panels=dict(a='Methods workflow; no quantitative claim.',
        b='All matched fresh-process external/Rust time ratios; startup included; all planned failures retained.',
        c='Actual formal API workflow wall time; complete fresh attempts only; failures and partial resumes remain in source tables.'),
    selection='All formal benchmark cells. All V1 three-parameter profiles, V2b B199/surrogate199, V2c six fits plus aggregation, case 5000-map and D2 200-tree/500-map workflows. P0 is supplementary only.',
    uncertainty='Raw timing attempts; repeated attempts are not independent datasets. No pooled confidence interval across tree sizes, hosts or workers.',
    scope='Actual calibration/case state counts only. Original 64/256-state complete-workflow benchmarks remain pending a dedicated protocol. Panel b memory remains supplementary. No standard-DEC comparison for Lagrange-ng: two-area v=2 stratum only. No optimizer speedup from fixed-likelihood timings.',
    legends=dict(b='software color; host marker; tips recorded in source and point horizontal offset', c='range-state color; reserved-CPU marker; counts show unique tasks observed/planned'),
    missingness='N=0 remains NA; failed calls and partial checkpoint resumes cannot become full-cost points. Returned calls with retained inner failures are marked.',
    raster_dpi=600,
    outputs=['figure1.pdf', 'figure1.svg', 'figure1.png', 'source_data', 'source_inventory', 'caption', 'QA'])


def workflow(ax):
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    ax.set_axis_off()
    ax.set_xlim(-.015, 1.015)
    ax.set_ylim(0, 1)
    boxes = {
        'input': (.005, .40, .17, .23, 'Tree + ranges'),
        'contract': (.21, .40, .20, .23, 'Probability rule\n+ root weights'),
        'fit': (.45, .40, .10, .23, 'Fit'),
        'parameter': (.60, .72, .21, .20, 'Profile / bootstrap'),
        'trees': (.60, .40, .21, .20, 'Model / tree averages\nExternal trees supported'),
        'maps': (.60, .08, .21, .20, 'Maps + precision'),
        'output': (.855, .40, .14, .23, 'Results\n+ failures')}
    for key, (x, y, w, h, label) in boxes.items():
        color = '#E7ECEF' if key == 'input' else '#DCEDEE' if key != 'output' else '#E3ECDC'
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.006,rounding_size=0.012',
            facecolor=color, edgecolor='#4B5960', linewidth=.7))
        ax.text(x+w/2, y+h/2, label, ha='center', va='center', fontsize=8)
    for source, target in [('input','contract'),('contract','fit'),('fit','parameter'),
                           ('fit','trees'),('fit','maps'),('parameter','output'),('trees','output'),('maps','output')]:
        x,y,w,h,_ = boxes[source]; tx,ty,tw,th,_ = boxes[target]
        ax.add_patch(FancyArrowPatch((x+w+.008,y+h/2),(tx-.008,ty+th/2),
            arrowstyle='-|>',mutation_scale=8,linewidth=.7,color='#39474D'))
    ax.text(.38,.10,'Bootstrap: simulate → refit',fontsize=8,ha='center')
    ax.text(0,1.11,'a',fontweight='bold',fontsize=10)


def benchmark_rows(inputs):
    rows=[]; grid=None
    for host,path in inputs:
        result=json.loads(path.read_text()); manifest=result['manifest']
        if result.get('formal') is not True or manifest.get('mode')!='formal':
            raise ValueError('Figure 1 requires formal benchmark results, never relabelled smoke data')
        if not result.get('source_stable') or not result.get('original_oracle_stable'):
            raise ValueError('benchmark changed source or oracle during measurement')
        this_grid=(manifest['tips'],manifest['states'],manifest['repeats'])
        if grid is not None and grid!=this_grid:
            raise ValueError('host benchmark grids must match')
        grid=this_grid
        expected=set(itertools.product(manifest['tips'],manifest['states'],range(manifest['repeats'])))
        observed=[]
        for record in result['records']:
            if record['preflight']: continue
            observed.append((record['tips'],record['states'],record['replicate']))
            expected_names={'original_biogeobears','revbayes'}|({'lagrange_ng'} if record['states']==4 else set())
            if set(record['comparisons'])!=expected_names:
                raise ValueError('missing planned software comparison')
            for software,comparison in record['comparisons'].items():
                rust=comparison['rust_reference']; left=record['execution'][software];right=record['execution'][rust]
                passed=comparison['alignment_passed'] and left['returncode']==right['returncode']==0
                ratio=left['wall_seconds']/right['wall_seconds'] if passed else None
                if passed and (not np.isfinite(ratio) or ratio<=0): raise ValueError('timing ratio must be finite and strictly positive')
                rows.append(dict(host=host,tips=record['tips'],states=record['states'],replicate=record['replicate'],
                    software=software,aligned=passed,time_ratio=ratio,
                    external_seconds=left['wall_seconds'],rust_seconds=right['wall_seconds'],
                    external_peak_rss_bytes=left['peak_rss_bytes'],rust_peak_rss_bytes=right['peak_rss_bytes'],
                    maximum_lnL_error=comparison['maximum_lnL_error'],
                    status='aligned' if passed else 'failed_or_unsupported',
                    input_tree_sha256=record['input_sha256']['tree.newick']))
        if len(observed)!=len(set(observed)) or set(observed)!=expected:
            raise ValueError('benchmark grid is incomplete or duplicated; retain failed cells as records')
    return rows, grid


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--benchmark',action='append',default=[],metavar='HOST=RESULTS_JSON')
    p.add_argument('--p0',action='append',default=[],metavar='HOST=RESULTS_JSON')
    p.add_argument('--workflow-input',type=Path,help='Locked full-suite timing inventory; never a P0 or smoke replacement')
    p.add_argument('--render',action='store_true');p.add_argument('--workflow-only',action='store_true')
    p.add_argument('--audit-scripts',type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    write_json(args.output/'figure_contract.json',CONTRACT)
    write_json(args.output/'workflow_input_template.json',dict(schema='rustbears.workflow_costs.inputs.v1',formal=False,
        budget=dict(path='LOCKED_BUDGET.json',sha256='PENDING'),suite='FORMAL_SUITE_DIRECTORY'))
    costs=inspect_suite(args.workflow_input) if args.workflow_input else None
    if costs is not None:export_costs(costs,args.output)
    if not args.render and not args.workflow_only:return
    if args.audit_scripts is None:p.error('rendering requires --audit-scripts for the final alignment check')
    if args.workflow_only:
        configure_plotting()
        import matplotlib.pyplot as plt
        fig=plt.figure(figsize=(180/25.4,55/25.4));ax=fig.add_axes([.025,.12,.95,.68]);workflow(ax)
        audit_alignment(fig,args.output,'workflow',args.audit_scripts)
        save_figure(fig,args.output,'workflow');plt.close(fig)
        write_json(args.output/'script_binding.json',{'sha256':sha(__file__),'helper_sha256':sha(Path(__file__).with_name('v12_figure_utils.py'))})
        (args.output/'caption.md').write_text('Figure 1a. An explicit probability convention and root specification connect fitting with parameter, model/tree, and mapping summaries. Bootstrap repeats simulation and fitting; external tree samples retain their provenance. Each workflow retains failed tasks and its numerical diagnostics. This methods diagram contains no measured performance or calibration result.\n')
        return
    parse=lambda items:[(item.split('=',1)[0],Path(item.split('=',1)[1])) for item in items]
    benchmarks=parse(args.benchmark);p0=parse(args.p0)
    if len(benchmarks)!=2 or len({x[0] for x in benchmarks})!=2:
        raise ValueError('provide two distinct matched benchmark hosts')
    if costs is None:raise ValueError('Figure 1c requires a verified formal --workflow-input; P0 cannot substitute')
    if len({r['host'] for r in costs['rows']})>1:raise ValueError('panel c requires one execution host; multi-host attempts need explicit separate panels')
    rows,grid=benchmark_rows(benchmarks);fitrows=[]
    for host,path in p0:
        data=json.loads(path.read_text())
        if data.get('status')!='complete' or data.get('smoke') or not data.get('p0_repetition_requirement_met'):
            raise ValueError('P0 inputs must be completed measurements')
        for group in data['groups']:
            for row in group['rows']:
                fitrows.append(dict(host=host,states=group['states'],workers=group['workers'],contract=group['contract'],
                    replicate=row['replicate'],seconds=row['elapsed_sec'],converged=row['converged'],
                    plotted=False))
    # Eligibility checks finish before opening any quantitative drawing backend.
    if any(sha(path)!=digest for path,digest in costs['source_files'].items()):raise ValueError('workflow evidence changed after inspection')
    configure_plotting()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    archive=[Path(path) for path in costs['source_files'] if Path(path).name in ('workflow_timings.jsonl','suite_plan.json','suite_ledger.json')]
    copy_evidence([path for _,path in benchmarks+p0]+[args.workflow_input]+archive,args.output)
    export_rows(args.output/'benchmark_source.csv',rows)
    if fitrows:export_rows(args.output/'supplementary_p0_source.csv',fitrows)
    fig_width_mm=180
    fig=plt.figure(figsize=(fig_width_mm/25.4,225/25.4))
    axa=fig.add_axes([.025,.80,.95,.155]);workflow(axa)
    axb=fig.add_axes([.31,.45,.63,.245]);axc=fig.add_axes([.31,.08,.63,.245])
    colors={'original_biogeobears':'#257C87','revbayes':'#A56529','lagrange_ng':'#767676'}
    names={'original_biogeobears':'BioGeoBEARS','revbayes':'RevBayes','lagrange_ng':'Lagrange-ng (v=2)'}
    states=grid[1];markers=['o','s'];tips=grid[0]
    for row in rows:
        if not row['aligned']:continue
        host=[x[0] for x in benchmarks].index(row['host']);software=list(colors).index(row['software'])
        x=states.index(row['states'])+(software-1)*.20+(host-.5)*.07+(tips.index(row['tips'])-(len(tips)-1)/2)*.018
        axb.scatter(x,row['time_ratio'],s=13,marker=markers[host],color=colors[row['software']],alpha=.65,linewidths=0)
    axb.axhline(1,color='#333333',linewidth=.7,linestyle=':')
    axb.set_yscale('log');axb.set_xticks(range(len(states)),list(map(str,states)))
    axb.set_xlabel('Range states');axb.set_ylabel('External / Rust process time')
    axb.yaxis.set_major_formatter(plt.FuncFormatter(lambda value,_:format(value,'g')))
    handles=[Line2D([],[],color=color,marker='o',linestyle='',label=names[name],markersize=4) for name,color in colors.items()]
    handles += [Line2D([],[],color='black',marker=markers[i],linestyle='',label=host,markersize=4) for i,(host,_) in enumerate(benchmarks)]
    axb.legend(handles=handles,loc='lower left',bbox_to_anchor=(0,1.015),frameon=False,ncol=3,columnspacing=.8,handletextpad=.3)
    failed=sum(not row['aligned'] for row in rows)
    axb.text(0,-.30,f'{failed}/{len(rows)} comparisons failed alignment',transform=axb.transAxes,fontsize=8)
    state_values=sorted({r['states'] for r in costs['points']});cpu_values=sorted({r['workers'] for r in costs['points']})
    palette=['#257C87','#A56529','#767676','#7C568D'];symbols=['o','s','^','D']
    if len(state_values)>len(palette) or len(cpu_values)>len(symbols):raise ValueError('new timing strata require a prospective legend revision')
    state_colors=dict(zip(state_values,palette));cpu_markers=dict(zip(cpu_values,symbols));labels=[]
    for i,(category,label) in enumerate(zip(WORKLOADS,LABELS)):
        chosen=[r for r in costs['points'] if r['workload']==category];counts=costs['category_counts'][category]
        labels.append(f'{label}\n{counts["unique_tasks_observed"]}/{counts["planned_tasks"]} tasks')
        for row,offset in zip(chosen,np.linspace(-.15,.15,len(chosen))):
            axc.scatter(row['wall_seconds'],i+offset,s=13,color=state_colors[row['states']],marker=cpu_markers[row['workers']],alpha=.55,linewidths=0)
            if row.get('n_failed',0) or row.get('n_trees_incomplete',0):axc.scatter(row['wall_seconds'],i+offset,marker='x',s=24,color='black')
        if not chosen:axc.text(.02,i,'NA',transform=axc.get_yaxis_transform(),va='center',fontsize=8)
    axc.set_xscale('log');axc.set_ylim(5.5,-.5);axc.set_yticks(range(6),labels)
    axc.xaxis.set_major_formatter(plt.FuncFormatter(lambda value,_:format(value,'g')))
    axc.set_xlabel('API workflow wall time (s)')
    handles=[Line2D([],[],color=color,marker='o',linestyle='',label=f'{state} states',markersize=4) for state,color in state_colors.items()]
    handles += [Line2D([],[],color='black',marker=marker,linestyle='',label=f'{cpu} CPUs reserved',markersize=4) for cpu,marker in cpu_markers.items()]
    if any(r.get('n_failed',0) or r.get('n_trees_incomplete',0) for r in costs['points']):handles.append(Line2D([],[],color='black',marker='x',linestyle='',label='Inner failures',markersize=4))
    if handles:axc.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.22,1.03),frameon=False,ncol=3,columnspacing=.8,handletextpad=.3)
    fig.text(.025,.745,'b',fontweight='bold',fontsize=10);fig.text(.025,.385,'c',fontweight='bold',fontsize=10)
    audit_alignment(fig,args.output,'figure1',args.audit_scripts)
    save_figure(fig,args.output,'figure1');plt.close(fig)
    (args.output/'caption.md').write_text('Figure 1. Analysis flow and measured computational workloads. (a) The probability convention travels with each analysis. (b) Each point is a mathematically aligned external/Rust fresh-process time ratio, including startup, imports, construction and one fixed likelihood evaluation; it is not an optimizer speedup. Lagrange-ng is restricted to the two-area v=2 stratum. Failed comparisons remain in source tables. (c) All fresh complete API attempts for the predeclared workflows. Profiles include d, e and j together; the LRT includes 199 null and 199 surrogate replicates; model averaging sums six sequential fits and aggregation only within the same invocation/dataset. Tree mapping denotes the planned 200-tree/500-map workload; crosses mark retained inner failures. Colors identify range states and symbols reserved CPUs, not measured CPU utilization. Counts are unique tasks with a full-cost observation divided by planned tasks; attempts can repeat after interruption. Timing includes internal checkpoint I/O and excludes later table exports. Partial resumes, raised calls, missing observations and every unplotted phase are retained, without imputation. Points span the registered tree sizes, scenarios and probability routes; no pooled mean, confidence interval or suite elapsed time is inferred. P0 fits appear only in supplementary evidence.\n')
    write_json(args.output/'script_binding.json',{'sha256':sha(__file__),'helper_sha256':sha(Path(__file__).with_name('v12_figure_utils.py')),'workflow_reader_sha256':sha(Path(__file__).with_name('v12_workflow_costs.py'))})


if __name__=='__main__':main()
