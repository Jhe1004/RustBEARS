#!/usr/bin/env python3
"""Plot completed P0 measurements only; this is supplementary pilot evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mac',type=Path,required=True);p.add_argument('--amd',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--audit-scripts',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(args.audit_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    contract=dict(question='Does the marginal cost of exact survival conditioning predict full optimization cost?',
        conclusion='Prepared likelihood cost and complete fit cost differ, especially at 128 states; the optimization workload must be timed directly.',
        archetype='quantitative grid',backend='python',figure_scope='supplementary P0 pilot, not the full formal performance benchmark',
        panels=dict(a='Prepared-evaluation cost provides the marginal arithmetic comparison.',b='Complete optimization tests whether that marginal comparison predicts realized cost.'),
        final_size_mm=[180,88],minimum_font_pt=8,formats=['pdf','svg','png'],
        glyph_and_alignment_requirement='8 pt floor; final plot rectangles within 1.5 pt; PDF text/collision audits and visual inspection',
        legend_coverage=dict(a='local four-entry contract/host legend',b='same local entries/order/style in reserved upper-left area'),
        inputs={key:dict(path=str(path.resolve()),sha256=sha(path)) for key,path in [('Mac',args.mac),('AMD',args.amd)]})
    (args.output/'figure_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    records=[];groups={}
    for host,path in [('Mac',args.mac),('AMD',args.amd)]:
        result=json.loads(path.read_text())
        if result['status']!='complete' or result['smoke'] or not result['p0_repetition_requirement_met']:
            raise ValueError('requires completed real P0 data, not smoke or planned values')
        for group in result['groups']:
            plotted=group['workers']==1
            if plotted:groups[host,group['contract'],group['states']]=group
            for row in group['rows']:
                records.append(dict(host=host,states=group['states'],contract=group['contract'],workers=group['workers'],
                    replicate=row['replicate'],seed=row['seed'],seconds=row['elapsed_sec'],n_evals=row['n_evals'],
                    prepared_evaluation_ms=group['prepared_evaluation']['seconds_per_evaluation']*1000,
                    converged=row['converged'],lnL=row['lnL'],plotted=plotted))
    with (args.output/'source_data.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
        'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7})
    fig,axes=plt.subplots(1,2,figsize=(7.0866141732,3.4645669291),gridspec_kw=dict(left=.10,right=.98,bottom=.19,top=.87,wspace=.36))
    states=(16,57,128);colors={'tip_survival':'#257C87','legacy_unconditioned':'#767676'}
    labels={'tip_survival':'Conditioned','legacy_unconditioned':'Legacy score'}
    offsets={('Mac','tip_survival'):-.12,('Mac','legacy_unconditioned'):-.04,
             ('AMD','tip_survival'):.04,('AMD','legacy_unconditioned'):.12}
    handles=[]
    for host,marker,style in [('Mac','o','-'),('AMD','s','--')]:
        for likelihood in ('tip_survival','legacy_unconditioned'):
            selected=[groups[host,likelihood,state] for state in states]
            x=np.arange(3)+offsets[host,likelihood];color=colors[likelihood]
            prepared=[g['prepared_evaluation']['seconds_per_evaluation']*1000 for g in selected]
            values=[np.array([row['elapsed_sec'] for row in g['rows']]) for g in selected]
            if not np.isfinite(prepared).all() or np.any(np.asarray(prepared)<=0) or any(not np.isfinite(v).all() or np.any(v<=0) for v in values):
                raise ValueError('log-scale timings must be finite and strictly positive')
            med=np.array([np.median(v) for v in values]);low=np.array([np.quantile(v,.1) for v in values]);high=np.array([np.quantile(v,.9) for v in values])
            axes[0].plot(x,prepared,color=color,marker=marker,linestyle=style,markersize=3.5,linewidth=1)
            axes[1].errorbar(x,med,yerr=[med-low,high-med],color=color,marker=marker,linestyle=style,
                markersize=3.5,linewidth=1,elinewidth=.7,capsize=2)
            # Every single-worker fit is visible; small deterministic offsets
            # spread points horizontally without changing measured values.
            for position,v in zip(x,values):
                axes[1].scatter(position+np.linspace(-.02,.02,len(v)),v,s=3,color=color,alpha=.4,linewidths=0,zorder=1)
            handles.append(Line2D([],[],color=color,marker=marker,linestyle=style,markersize=3.5,
                linewidth=1,label=f'{host}: {labels[likelihood]}'))
    for letter,ax,title,ylabel in zip('ab',axes,['Prepared evaluation','Complete optimization'],['Evaluation time (ms)','Fit time (s)']):
        ax.set_yscale('log');ax.set_xlim(-.3,2.3);ax.set_xticks(range(3),list(map(str,states)));ax.set_xlabel('Range states')
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda value,_position:format(value,'g')))
        ax.set_ylabel(ylabel);ax.set_title(title,fontsize=9,pad=10)
        ax.annotate(letter,xy=(0,1),xycoords='axes fraction',xytext=(-24,10),textcoords='offset points',fontsize=10,fontweight='bold')
        ax.legend(handles=handles,loc='upper left',frameon=False,handlelength=1.8,borderaxespad=.4,labelspacing=.35)
        ax.tick_params(width=.7,length=3)
    axes[0].set_ylim(.01,50);axes[1].set_ylim(.07,80)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,json_out=str(args.output/'p0_costs.alignment.json'),
        overlay_svg=str(args.output/'p0_costs.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
    fig.savefig(args.output/'p0_costs.pdf')
    fig.savefig(args.output/'p0_costs.svg')
    fig.savefig(args.output/'p0_costs.png',dpi=600)
    plt.close(fig)
    caption=('Supplementary pilot figure. Prepared evaluation and complete optimization have different cost profiles. '
        '(a) Mean time per evaluation from one timed batch of 100 calls after 10 warm-up calls; no independent timing-batch interval is available. '
        '(b) Median and 10th–90th percentiles of 20 complete DEC+J fits to the same fixed four-tip dataset; small points show every fit. '
        'Symbols distinguish Mac (circles) and AMD (squares); color distinguishes exact whole-tree tip-survival conditioning from the legacy unconditioned score. '
        'All shown groups use one outer worker, one native thread and one BLAS thread, with maxiter=500. '
        'The AMD machine was shared and had a load average near 241 on 256 logical CPUs; these data do not establish an isolated hardware ranking. '
        'The complete source table retains 480 fits, including 240 fits with two/four outer workers omitted from this contract-cost view. '
        'All 240 displayed fits converged. Prepared costs must not be extrapolated to full optimization, larger trees, or formal inferential coverage.\n')
    (args.output/'caption.md').write_text(caption)
    evidence=dict(source_sha256=sha(__file__),input_sha256={key:sha(path) for key,path in [('Mac',args.mac),('AMD',args.amd)]},
        input_fits=len(records),displayed_fits=sum(row['plotted'] for row in records),
        omitted_from_plot_rule='workers>1: retained in source table but not comparable to the serial legacy-score groups',
        outputs={f'p0_costs.{ext}':sha(args.output/f'p0_costs.{ext}') for ext in ('pdf','svg','png')},
        status='rendered_awaiting_external_text_collision_and_visual_QA')
    (args.output/'provenance.json').write_text(json.dumps(evidence,indent=2)+'\n')


if __name__=='__main__':main()
