#!/usr/bin/env python3
"""Plan Figure 1; render only actual formal benchmark and measured P0 inputs.

The workflow-only option draws a methods diagram without inventing data panels.
Matched process ratios include startup and are never called kernel speedups.
"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from v12_figure_utils import (audit_alignment, configure_plotting, copy_evidence,
                              export_rows, save_figure, sha, write_json)

CONTRACT = dict(figure='1', backend='python', final_size_mm=[180, 158], minimum_font_pt=8,
    question='How do explicit probability conventions connect analysis and measured computational cost?',
    panels=dict(a='Methods workflow; no quantitative claim.',
        b='All matched fresh-process external/Rust time ratios; startup included; all planned failures retained.',
        c='All single-worker P0 conditioned fits; exploratory four-tip workload, not larger-tree extrapolation.'),
    selection='All formal benchmark cells. All single-worker conditioned P0 replicates; legacy and parallel P0 groups retained in source data and supplementary P0 figure.',
    uncertainty='Raw independent timing repetitions; no pooled inferential confidence interval across tree sizes or hosts.',
    scope='No standard-DEC comparison for Lagrange-ng: two-area v=2 stratum only. No optimizer speedup from fixed-likelihood timings.',
    legends=dict(b='software color; host marker; tips recorded in source and point horizontal offset', c='host color/marker'),
    outputs=['pdf', 'svg', 'png', 'source_data', 'source_inventory', 'caption', 'QA'])


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
                if passed and (not np.isfinite(ratio) or ratio<=0): raise ValueError('invalid timing ratio')
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
    p.add_argument('--render',action='store_true');p.add_argument('--workflow-only',action='store_true')
    p.add_argument('--audit-scripts',type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    write_json(args.output/'figure_contract.json',CONTRACT)
    if not args.render and not args.workflow_only:return
    if args.audit_scripts is None:p.error('rendering requires --audit-scripts for the final alignment check')
    configure_plotting()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if args.workflow_only:
        fig=plt.figure(figsize=(180/25.4,55/25.4));ax=fig.add_axes([.025,.12,.95,.68]);workflow(ax)
        audit_alignment(fig,args.output,'workflow',args.audit_scripts)
        save_figure(fig,args.output,'workflow');plt.close(fig)
        write_json(args.output/'script_binding.json',{'sha256':sha(__file__),'helper_sha256':sha(Path(__file__).with_name('v12_figure_utils.py'))})
        (args.output/'caption.md').write_text('Figure 1a. An explicit probability convention and root specification connect fitting with parameter, model/tree, and mapping summaries. Bootstrap repeats simulation and fitting; external tree samples retain their provenance. Each workflow retains failed tasks and its numerical diagnostics. This methods diagram contains no measured performance or calibration result.\n')
        return
    parse=lambda items:[(item.split('=',1)[0],Path(item.split('=',1)[1])) for item in items]
    benchmarks=parse(args.benchmark);p0=parse(args.p0)
    if len(benchmarks)!=2 or len({x[0] for x in benchmarks})!=2 or {x[0] for x in benchmarks}!={x[0] for x in p0}:
        raise ValueError('provide two distinct matched hosts for both --benchmark and --p0')
    rows,grid=benchmark_rows(benchmarks);fitrows=[]
    for host,path in p0:
        data=json.loads(path.read_text())
        if data.get('status')!='complete' or data.get('smoke') or not data.get('p0_repetition_requirement_met'):
            raise ValueError('P0 inputs must be completed measurements')
        for group in data['groups']:
            for row in group['rows']:
                fitrows.append(dict(host=host,states=group['states'],workers=group['workers'],contract=group['contract'],
                    replicate=row['replicate'],seconds=row['elapsed_sec'],converged=row['converged'],
                    plotted=group['workers']==1 and group['contract']=='tip_survival'))
    copy_evidence([path for _,path in benchmarks+p0],args.output)
    export_rows(args.output/'benchmark_source.csv',rows);export_rows(args.output/'p0_source.csv',fitrows)
    fig=plt.figure(figsize=(180/25.4,158/25.4))
    axa=fig.add_axes([.025,.70,.95,.235]);workflow(axa)
    axb=fig.add_axes([.10,.13,.37,.43]);axc=fig.add_axes([.61,.13,.37,.43])
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
    axb.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.12,1.015),frameon=False,ncol=2,columnspacing=.8,handletextpad=.3)
    failed=sum(not row['aligned'] for row in rows)
    axb.text(0,-.24,f'{failed}/{len(rows)} comparisons failed alignment',transform=axb.transAxes,fontsize=8)
    hostcolors=['#257C87','#A56529'];fitstates=sorted({r['states'] for r in fitrows if r['plotted']})
    for i,(host,_) in enumerate(p0):
        for state in fitstates:
            chosen=[r for r in fitrows if r['host']==host and r['states']==state and r['plotted']]
            values=np.array([r['seconds'] for r in chosen]);x=fitstates.index(state)+(i-.5)*.17
            axc.scatter(x+np.linspace(-.035,.035,len(values)),values,s=11,color=hostcolors[i],marker=markers[i],alpha=.6,linewidths=0)
            for row,value,offset in zip(chosen,values,np.linspace(-.035,.035,len(values))):
                if not row['converged']:axc.scatter(x+offset,value,marker='x',s=24,color='black')
    axc.set_yscale('log');axc.yaxis.set_major_formatter(plt.FuncFormatter(lambda value,_:format(value,'g')))
    axc.set_xticks(range(len(fitstates)),list(map(str,fitstates)));axc.set_xlabel('Range states');axc.set_ylabel('Conditioned fit time (s)')
    axc.legend(handles=[Line2D([],[],color=hostcolors[i],marker=markers[i],linestyle='',label=host) for i,(host,_) in enumerate(p0)],frameon=False,loc='lower left',bbox_to_anchor=(0,1.015))
    for letter,ax in [('b',axb),('c',axc)]:ax.text(-.19,1.28,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
    audit_alignment(fig,args.output,'figure1',args.audit_scripts)
    save_figure(fig,args.output,'figure1');plt.close(fig)
    (args.output/'caption.md').write_text('Figure 1. Analysis flow and distinct computational workloads. (a) The probability convention travels with each analysis. (b) Each point is one matched fresh-process timing ratio, including interpreter startup, imports, model construction and one fixed-parameter likelihood. The reference and Rust Q, C, root weights and states are aligned in untimed checks, and each timed likelihood agrees within 1e-8 before a ratio is shown. Point offsets retain all tree sizes and repetitions; no pooled confidence interval is implied. Lagrange-ng is only the two-area v=2 stratum, not standard DEC. Memory, failure and unsupported results remain in the full source table. (c) Every single-worker P0 DEC+J fit to the fixed four-tip input, with exact tip-survival conditioning; crosses identify nonconvergence. This measured pilot workload does not predict optimization cost for larger trees or complete inferential procedures. Host load and executable hashes are archived with the input manifests.\n')
    write_json(args.output/'script_binding.json',{'sha256':sha(__file__),'helper_sha256':sha(Path(__file__).with_name('v12_figure_utils.py'))})


if __name__=='__main__':main()
