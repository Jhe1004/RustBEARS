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
    from .v12_workflow_costs import WORKLOADS, LABELS, inspect_suite, inspect_dedicated_benchmark, export_costs
else:
    from v12_figure_utils import (audit_alignment, configure_plotting, copy_evidence,
                                  export_rows, save_figure, sha, write_json)
    from v12_workflow_costs import WORKLOADS, LABELS, inspect_suite, inspect_dedicated_benchmark, export_costs

CONTRACT = dict(figure='1', backend='python', final_size_mm=[180, 225], minimum_font_pt=8,
    question='How do explicit probability conventions connect analysis and measured computational cost?',
    panels=dict(a='Methods workflow; no quantitative claim.',
        b='Two sibling axes: all matched fresh-process external/Rust wall-time and peak-RSS ratios; all planned failures retained.',
        c='Actual formal API workflow wall time; complete fresh attempts only; failures and partial resumes remain in source tables.'),
    selection='All formal external benchmark cells. Main c uses the separately locked 16/64/256-state six-workflow cost benchmark. Actual calibration/case workflow timings and P0 are supplementary only.',
    uncertainty='Raw timing attempts; repeated attempts are not independent datasets. No pooled confidence interval across tree sizes, hosts or workers.',
    scope='Dedicated 16/64/256-state full-workflow protocol implemented; formal data and rendering pending. Calibration/case costs do not substitute for this grid. Peak RSS is the whole fresh process measured by wait4, including interpreter/runtime, not kernel allocation. No standard-DEC comparison for Lagrange-ng: two-area v=2 stratum only. No optimizer speedup from fixed-likelihood timings.',
    legends=dict(b='one bounded panel header shared by both axes: software color; host marker; unavailable/failed status markers; tips and repeats retained by horizontal offsets and source rows', c='locked range-state color and reserved-CPU marker; wholly missing planned strata remain in the legend; counts show unique tasks observed/planned'),
    missingness='N=0 remains NA; failed calls and partial checkpoint resumes cannot become full-cost points. Returned bootstrap/tree calls with retained replicate/tree failures are marked. Profile box censoring, disconnected support and missing Hessian SE remain distinct source-table diagnostics, not numerical failures.',
    raster_dpi=600,
    outputs=['figure1.pdf', 'figure1.svg', 'figure1.png', 'figure1.tiff', 'source_data', 'source_inventory', 'caption', 'QA'])

BENCHMARK_COLORS={'original_biogeobears':'#257C87','revbayes':'#A56529','lagrange_ng':'#767676'}
BENCHMARK_NAMES={'original_biogeobears':'BioGeoBEARS','revbayes':'RevBayes','lagrange_ng':'Lagrange-ng (v=2)'}


def workflow_legend_spec(costs):
    """Keep encodings fixed by the locked design, never by surviving points."""
    states=costs['settings']['states'];workers=costs['settings']['workers']
    palette=['#257C87','#A56529','#767676','#7C568D'];symbols=['o','s','^','D']
    if not states or not workers or len(set(states))!=len(states) or len(set(workers))!=len(workers):
        raise ValueError('workflow legend requires unique planned states and workers')
    if len(states)>len(palette) or len(workers)>len(symbols):raise ValueError('new timing strata require a prospective legend revision')
    if any(p['states'] not in states or p['workers'] not in workers for p in costs['points']):raise ValueError('unplanned cost point in workflow legend')
    counts=costs['stratum_counts']
    def label(field,value,text):
        selected=[r for r in counts if r[field]==value]
        planned=sum(r['planned_tasks'] for r in selected);observed=sum(r['unique_tasks_observed'] for r in selected)
        if planned<=0:raise ValueError('planned legend stratum lacks its denominator')
        return text+(f' (0/{planned})' if observed==0 else '')
    return dict(state_colors=dict(zip(states,palette)),cpu_markers=dict(zip(workers,symbols)),
        state_labels={s:label('states',s,f'{s} states') for s in states},
        cpu_labels={w:label('workers',w,f'{w} CPUs reserved') for w in workers})


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
    ax.text(0,1.11,'a)',fontweight='bold',fontsize=10)


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
                timeout=left.get('timed_out',False) or right.get('timed_out',False)
                passed=comparison['alignment_passed'] and left['returncode']==right['returncode']==0 and not timeout
                for execution in (left,right):
                    for key in ('wall_seconds','peak_rss_bytes'):
                        value=execution.get(key)
                        if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or value<0):
                            raise ValueError(f'{key} must be a finite nonnegative number or unavailable')
                    if passed and (execution.get('wall_seconds') is None or execution['wall_seconds']<=0):
                        raise ValueError('aligned process time must be finite and strictly positive')
                ratio=left['wall_seconds']/right['wall_seconds'] if passed else None
                if passed and (not np.isfinite(ratio) or ratio<=0): raise ValueError('timing ratio must be finite and strictly positive')
                # An absent/zero RSS means no usable measurement, never zero memory.
                rss_left,rss_right=left.get('peak_rss_bytes'),right.get('peak_rss_bytes')
                rss_available=all(value is not None and value>0 for value in (rss_left,rss_right))
                memory_ratio=rss_left/rss_right if passed and rss_available else None
                if memory_ratio is not None and (not np.isfinite(memory_ratio) or memory_ratio<=0):
                    raise ValueError('memory ratio must be finite and strictly positive')
                rows.append(dict(host=host,tips=record['tips'],states=record['states'],replicate=record['replicate'],
                    software=software,aligned=passed,time_ratio=ratio,memory_ratio=memory_ratio,
                    external_seconds=left.get('wall_seconds'),rust_seconds=right.get('wall_seconds'),
                    external_peak_rss_bytes=rss_left,rust_peak_rss_bytes=rss_right,
                    external_peak_rss_mib=rss_left/(1024**2) if rss_left is not None and rss_left>0 else None,
                    rust_peak_rss_mib=rss_right/(1024**2) if rss_right is not None and rss_right>0 else None,
                    memory_available=rss_available,external_returncode=left['returncode'],rust_returncode=right['returncode'],
                    external_timed_out=left.get('timed_out',False),rust_timed_out=right.get('timed_out',False),
                    maximum_lnL_error=comparison['maximum_lnL_error'],
                    status='aligned' if passed else 'timeout' if timeout else 'failed_or_unsupported',
                    input_tree_sha256=record['input_sha256']['tree.newick']))
        if len(observed)!=len(set(observed)) or set(observed)!=expected:
            raise ValueError('benchmark grid is incomplete or duplicated; retain failed cells as records')
    return rows, grid


def benchmark_counts(rows):
    """Denominators count every planned comparison, not only successful pairs."""
    return dict(planned=len(rows),time_available=sum(r['time_ratio'] is not None for r in rows),
        memory_ratio_available=sum(r['memory_ratio'] is not None for r in rows),
        timeouts=sum(r['status']=='timeout' for r in rows),
        other_failed=sum(r['status']=='failed_or_unsupported' for r in rows),
        missing_rss_pairs=sum(not r['memory_available'] for r in rows))


def draw_benchmark_panel(fig,rows,grid,hosts):
    """Matched process cost ratios; missing values use a non-numeric status band."""
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    # The header/footer and both equal plot rectangles belong to one panel b.
    # The panel ownership is conceptual; journal production prohibits a visible
    # box around the complete graph. Shared header and equal axes decode it.
    header=fig.add_axes([.075,.705,.87,.068]);header.set_axis_off()
    axes=[fig.add_axes([.115,.482,.345,.195]),fig.add_axes([.615,.482,.345,.195])]
    states,tips=grid[1],grid[0];markers=['o','s']
    counts=benchmark_counts(rows)
    for ax,metric,label in zip(axes,['time_ratio','memory_ratio'],
            ['External / Rust\nprocess wall time','External / Rust\npeak process RSS']):
        for row in rows:
            host=hosts.index(row['host']);software=list(BENCHMARK_COLORS).index(row['software'])
            x=states.index(row['states'])+(software-1)*.20+(host-.5)*.07
            x+=(tips.index(row['tips'])-(len(tips)-1)/2)*.018
            x+=(row['replicate']-(grid[2]-1)/2)*(.012/max(1,grid[2]-1))
            color=BENCHMARK_COLORS[row['software']]
            if row[metric] is not None:
                ax.scatter(x,row[metric],s=13,marker=markers[host],color=color,alpha=.65,linewidths=0)
            else:
                marker='+' if row['status']=='timeout' else 'x' if not row['aligned'] else '|'
                ax.scatter(x,-.14,s=18,marker=marker,color=color,linewidths=.8,
                    transform=ax.get_xaxis_transform(),clip_on=False)
        ax.axhline(1,color='#333333',linewidth=.7,linestyle=':')
        ax.set_yscale('log');ax.set_xlim(-.55,len(states)-.45)
        ax.set_xticks(range(len(states)),list(map(str,states)));ax.set_ylabel(label)
        ax.set_xlabel('Range states',labelpad=20)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda value,_:format(value,'g')))
        ax.set_title('Time' if metric=='time_ratio' else 'Peak RSS',loc='left',fontsize=8,pad=5)
    handles=[Line2D([],[],color=color,marker='o',linestyle='',label=BENCHMARK_NAMES[name],markersize=4)
             for name,color in BENCHMARK_COLORS.items()]
    handles += [Line2D([],[],color='black',marker=markers[i],linestyle='',label=host,markersize=4) for i,host in enumerate(hosts)]
    for number,marker,label in [(counts['timeouts'],'+','Timeout'),(counts['other_failed'],'x','Failed / unaligned'),
                                 (counts['missing_rss_pairs'],'|','RSS unavailable')]:
        if number:handles.append(Line2D([],[],color='black',marker=marker,linestyle='',label=label,markersize=4))
    header.legend(handles=handles,loc='center',frameon=False,ncol=3,columnspacing=.8,handletextpad=.3)
    fig.text(.025,.795,'b)',fontweight='bold',fontsize=10)
    fig.text(.055,.414,f"Time {counts['time_available']}/{counts['planned']}; RSS ratios {counts['memory_ratio_available']}/{counts['planned']}",fontsize=8)
    fig.text(.055,.399,f"Timeout {counts['timeouts']}; other failed {counts['other_failed']}; RSS unavailable {counts['missing_rss_pairs']}",fontsize=8)
    return counts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--benchmark',action='append',default=[],metavar='HOST=RESULTS_JSON')
    p.add_argument('--p0',action='append',default=[],metavar='HOST=RESULTS_JSON')
    group=p.add_mutually_exclusive_group()
    group.add_argument('--workflow-input',type=Path,help='Supplementary actual calibration/case cost inventory')
    group.add_argument('--workflow-benchmark-input',type=Path,help='Dedicated formal 16/64/256-state cost inventory for main panel c')
    p.add_argument('--render',action='store_true');p.add_argument('--workflow-only',action='store_true')
    p.add_argument('--audit-scripts',type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    write_json(args.output/'figure_contract.json',CONTRACT)
    write_json(args.output/'workflow_input_template.json',dict(schema='rustbears.workflow_costs.inputs.v1',formal=False,
        budget=dict(path='LOCKED_BUDGET.json',sha256='PENDING'),suite='FORMAL_SUITE_DIRECTORY'))
    write_json(args.output/'workflow_benchmark_input_template.json',dict(schema='rustbears.workflow_benchmark.inputs.v1',formal=False,
        budget=dict(path='LOCKED_DEDICATED_COST_BUDGET.json',sha256='PENDING'),benchmark='DEDICATED_BENCHMARK_DIRECTORY'))
    costs=inspect_dedicated_benchmark(args.workflow_benchmark_input) if args.workflow_benchmark_input else inspect_suite(args.workflow_input) if args.workflow_input else None
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
    if costs is None or args.workflow_benchmark_input is None:raise ValueError('Main Figure 1c requires a formal --workflow-benchmark-input; calibration costs and P0 cannot substitute')
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
    copy_evidence([path for _,path in benchmarks+p0]+[args.workflow_benchmark_input]+archive,args.output)
    export_rows(args.output/'benchmark_source.csv',rows)
    if fitrows:export_rows(args.output/'supplementary_p0_source.csv',fitrows)
    fig_width_mm,fig_height_mm=CONTRACT['final_size_mm']
    fig=plt.figure(figsize=(fig_width_mm/25.4,fig_height_mm/25.4))
    axa=fig.add_axes([.025,.825,.95,.135]);workflow(axa)
    benchmark_summary=draw_benchmark_panel(fig,rows,grid,[host for host,_ in benchmarks])
    write_json(args.output/'benchmark_counts.json',benchmark_summary)
    write_json(args.output/'legend_coverage.json',dict(b=dict(axes=['b1','b2'],
        owner='bounded panel b',position='shared internal header',
        mappings=['software color','host marker','non-numeric status markers']),
        c=dict(owner='panel c',position='local header',mappings=['state color','reserved-CPU marker'])))
    axc=fig.add_axes([.31,.075,.63,.24])
    legend=workflow_legend_spec(costs)
    state_colors=legend['state_colors'];cpu_markers=legend['cpu_markers'];labels=[]
    write_json(args.output/'workflow_legend.json',legend)
    order=costs.get('workload_order',WORKLOADS);cost_labels=costs.get('workload_labels',LABELS)
    for i,(category,label) in enumerate(zip(order,cost_labels)):
        chosen=[r for r in costs['points'] if r['workload']==category];counts=costs['category_counts'][category]
        labels.append(f'{label}\n{counts["unique_tasks_observed"]}/{counts["planned_tasks"]} tasks')
        for row,offset in zip(chosen,np.linspace(-.15,.15,len(chosen))):
            axc.scatter(row['wall_seconds'],i+offset,s=13,color=state_colors[row['states']],marker=cpu_markers[row['workers']],alpha=.55,linewidths=0)
            if row.get('n_failed',0) or row.get('n_trees_incomplete',0):axc.scatter(row['wall_seconds'],i+offset,marker='x',s=24,color='black')
        if not chosen:axc.text(.02,i,'NA',transform=axc.get_yaxis_transform(),va='center',fontsize=8)
    axc.set_xscale('log');axc.set_ylim(5.5,-.5);axc.set_yticks(range(6),labels)
    axc.xaxis.set_major_formatter(plt.FuncFormatter(lambda value,_:format(value,'g')))
    axc.set_xlabel('API workflow wall time (s)')
    handles=[Line2D([],[],color=color,marker='o',linestyle='',label=legend['state_labels'][state],markersize=4) for state,color in state_colors.items()]
    handles += [Line2D([],[],color='black',marker=marker,linestyle='',label=legend['cpu_labels'][cpu],markersize=4) for cpu,marker in cpu_markers.items()]
    if any(r.get('n_failed',0) or r.get('n_trees_incomplete',0) for r in costs['points']):handles.append(Line2D([],[],color='black',marker='x',linestyle='',label='Bootstrap/tree failures',markersize=4))
    if handles:axc.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.22,1.03),frameon=False,ncol=3,columnspacing=.8,handletextpad=.3)
    fig.text(.025,.355,'c)',fontweight='bold',fontsize=10)
    audit_alignment(fig,args.output,'figure1',args.audit_scripts)
    save_figure(fig,args.output,'figure1');plt.close(fig)
    (args.output/'caption.md').write_text('Figure 1. Analysis flow and measured computational workloads. (a) The probability convention travels with each analysis. (b) Left, time; right, peak RSS. Each point is a matched external/Rust ratio for whole fresh-process wall time or peak resident memory (RSS), respectively. Both measurements include interpreter/runtime startup, imports, construction and one fixed likelihood evaluation; neither establishes optimizer speedup or kernel-allocation efficiency. Peak RSS is measured separately for each process using wait4. All registered tips, states, hosts and repeats are retained; horizontal offsets separate software, host, tip count and repeat. Lagrange-ng is restricted to the two-area v=2 stratum. The shared header owns both axes. A missing or zero RSS yields an unavailable ratio, never a zero-memory value. Timeout, failed/unaligned and unavailable marks appear in a non-numeric strip below each axis; they do not encode a numeric ratio. Counts retain every planned comparison, and RSS-unavailable counts may overlap failures. Raw seconds, bytes, return codes and timeout flags remain in the source table. (c) All fresh complete API attempts for the predeclared workflows. Profiles include d, e and j together; the dedicated LRT includes the registered null replicates without surrogate power. Six-model costs include every initial fit and averaging within the same phase. Branch-scaled tree copies define a computational workload, not a posterior sample; map and tree budgets are displayed on the axes. Crosses mark retained bootstrap-replicate or tree failures. A returned profile API is distinct from a numerical certification failure: box-censored endpoints, disconnected support, and unavailable Hessian standard errors are preserved separately in the source tables, without reclassifying them as failed calls. Colors identify all locked range-state levels and symbols all reserved-CPU levels, including wholly missing levels labelled 0/N; they do not encode measured CPU utilization. Per-workflow/state/worker observed and planned counts remain in the source table. Counts are unique tasks with a full-cost observation divided by planned tasks; attempts can repeat after interruption. Timing includes internal checkpoint I/O and excludes later table exports. Partial resumes, raised calls, missing observations and every unplotted phase are retained, without imputation. Points span the locked state/worker grid and repetitions under tip-survival conditioning; no pooled mean, confidence interval or suite elapsed time is inferred. Partial or timed-out phases have no full-cost point. Actual calibration/case costs and P0 fits are supplementary evidence.\n')
    binding={'sha256':sha(__file__),'helper_sha256':sha(Path(__file__).with_name('v12_figure_utils.py')),'workflow_reader_sha256':sha(Path(__file__).with_name('v12_workflow_costs.py'))}
    if costs.get('scope')=='core_v1':
        with (args.output/'caption.md').open('a') as stream:stream.write('\nCore scope amendment: '+costs['repetition_interpretation']+'\n')
        binding.update(scope='core_v1',collection_reader_sha256=sha(Path(__file__).with_name('v12_core_reporting.py')))
    write_json(args.output/'script_binding.json',binding)


if __name__=='__main__':main()
