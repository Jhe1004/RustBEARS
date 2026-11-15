#!/usr/bin/env python3
"""Plan Figure 3; render the locked, fully processed D2 empirical analysis.

This script cannot promote smoke output to a formal empirical figure. The
summary tree and every planned model/contract remain visible in source tables.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from v12_figure_utils import (audit_alignment, configure_plotting, copy_evidence,
                              save_figure, sha, write_json)

CONTRACT=dict(figure='3',backend='python',size_mm=[180,245],minimum_font_pt=8,
    question='How do parameter uncertainty, a fixed external tree sample and mapping precision qualify one empirical reconstruction?',
    panels=dict(a='All 57 tips and conditioned model-averaged reference-node probabilities.',
        b='j profiles for all three +J models under both contracts; censoring and noninferential legacy support labelled.',
        c='All three nested bootstrap comparisons under both contracts; failure bounds and Monte Carlo intervals distinguished.',
        d='All internal-node nonempty-range probabilities: reference vs strict monophyletic external-tree averages; missing cases counted.',
        e='d/e/j event-count quantiles and numerical mean intervals, by contract; all eight event types retained in source tables.'),
    predeclared_selection='D2 Trimezieae, 57 tips, E/N areas; no choice based on empirical estimates. d/e/j selected for their dispersal/extirpation/founder roles.',
    inference='External deposited trees are a fixed sample. AICc weights are not model posterior probabilities. Legacy profiles are support intervals, not calibrated confidence intervals.',
    QA='Editable PDF/SVG, 600 dpi PNG, exact sources/hashes, alignment, 8 pt PDF audit, collision audit and final-size visual inspection.')


def load_case(case, tree_path, budget_path):
    manifest=json.loads((case/'manifest.json').read_text());result=json.loads((case/'results.json').read_text())
    budget=json.loads(budget_path.read_text())
    if manifest.get('formal') is not True or result.get('formal') is not True or manifest.get('mode')!='formal':
        raise ValueError('Figure 3 rejects smoke/pilot results')
    if result.get('status') not in ('complete','incomplete_profiles','complete_with_tree_failures'):
        raise ValueError('case workflow has unprocessed tasks')
    if manifest['case_id']!='D2_Trimezieae' or manifest['tip_count']!=57 or manifest['areas']!=['E','N']:
        raise ValueError('inputs differ from the prospectively selected D2 case')
    lock=manifest.get('locked_budget') or {}
    if budget.get('status')!='locked' or lock.get('sha256')!=sha(budget_path) or lock.get('payload')!=budget:
        raise ValueError('case is not bound to this locked budget')
    if sha(tree_path)!=manifest['input_sha256']['tree']:raise ValueError('reference tree hash mismatch')
    jobs=[j for j in budget['jobs'] if j.get('analysis')=='d2' and j.get('case_id')==manifest['case_id']]
    if len(jobs)!=1:raise ValueError('budget must contain exactly one D2 case')
    for name in ['n_boot','bsm_maps','max_trees','bsm_maps_per_tree','seed']:
        if jobs[0][name]!=manifest[name]:raise ValueError('budget setting mismatch: '+name)
    if budget['source_sha256']['scripts/run_uncertainty_case_study.py']!=manifest['workflow_script_sha256']:
        raise ValueError('case script differs from the frozen source')
    for name in ['tip_survival','legacy_unconditioned']:
        contract=result['contracts'][name]
        if len(contract['profiles'])!=6 or len(contract['bootstrap'])!=3 or contract['multitree']['n_planned']!=manifest['max_trees']:
            raise ValueError('case is missing planned tasks')
        if contract['multitree']['n_complete'] and contract['multitree']['joint_mapping'].get('tree_target')!='fixed_sample':
            raise ValueError('D2 deposited trees must retain the fixed-sample target')
    # The same read-only terminal verifier used by the execution suite checks
    # every actual replicate/tree/map identity and its retained failure status.
    # A complete-looking summary file is never sufficient for a paper figure.
    from v12_case_evidence import verify_case
    verify_case(case,jobs[0],budget,sha(budget_path))
    return manifest,result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--case',type=Path)
    parser.add_argument('--tree',type=Path);parser.add_argument('--budget-json',type=Path)
    parser.add_argument('--render',action='store_true');parser.add_argument('--audit-scripts',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    write_json(args.output/'figure_contract.json',CONTRACT)
    if not args.render:return
    if any(x is None for x in [args.case,args.tree,args.budget_json,args.audit_scripts]):
        parser.error('--render needs --case, --tree, --budget-json and --audit-scripts')
    manifest,result=load_case(args.case,args.tree,args.budget_json)
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    from rustbears.io.newick import read_newick
    tree=read_newick(args.tree)
    contracts=['tip_survival','legacy_unconditioned'];colors=['#257C87','#767676']
    paths=[args.case/'manifest.json',args.case/'results.json',args.tree,args.budget_json]
    profiles=[];events=[];boot=[]
    for contract in contracts:
        directory=args.case/contract
        for filename in ['models.csv','multitree_probability_bounds.csv','ancestral_arrays.npz','fitted_parameters.json','profile_status.json']:
            paths.append(directory/filename)
        fitted=json.loads((directory/'fitted_parameters.json').read_text())
        statuses=json.loads((directory/'profile_status.json').read_text())
        if statuses!=result['contracts'][contract]['profiles']:
            raise ValueError('profile status ledger and case report differ')
        for index in range(1,7):
            path=directory/f'profile_{index}.csv';status=statuses[index-1]
            if status['status']=='complete':
                paths.append(path);table=pd.read_csv(path)
            elif status['status']=='numerically_uncertified':
                # An explicitly processed failure remains an unavailable row;
                # it cannot disappear simply because no interval CSV exists.
                table=pd.DataFrame([dict(param=param,mle=fitted[index-1]['params'][param],lower=np.nan,upper=np.nan,
                    status_lower='numerically_uncertified',status_upper='numerically_uncertified',
                    interval_kind='unavailable',error=status.get('error'))
                    for param in manifest['model_parameter_bounds'][manifest['models'][index-1]]])
            else:raise ValueError('unprocessed profile status')
            table['model']=manifest['models'][index-1];table['contract']=contract
            profiles.append(table)
        event_path=directory/'joint_event_precision.csv'
        complete=result['contracts'][contract]['multitree']['n_complete']
        if complete:
            event=pd.read_csv(event_path)
            for filename in ['joint_event_precision.csv','joint_route_comparisons.csv','joint_bsm_summary.json']:
                paths.append(directory/filename)
        else:
            event=pd.DataFrame([dict(event_type=kind,mean=np.nan,rpi_025=np.nan,rpi_975=np.nan,
                mean_mc_ci_lower=np.nan,mean_mc_ci_upper=np.nan,assessed=False,assessment='no_complete_trees')
                for kind in ['d','e','a','j','y','s','v','t']])
        event['contract']=contract;events.append(event)
        for row in result['contracts'][contract]['bootstrap']:boot.append(dict(contract=contract,**row))
    profiles=pd.concat(profiles,ignore_index=True);events=pd.concat(events,ignore_index=True)
    profiles.to_csv(args.output/'profile_source.csv',index=False);events.to_csv(args.output/'event_source.csv',index=False)
    pd.DataFrame(boot).to_csv(args.output/'bootstrap_source.csv',index=False)
    copy_evidence(paths,args.output)
    with np.load(args.case/'tip_survival/ancestral_arrays.npz',allow_pickle=False) as archive:prob=archive['average_top']
    if prob.shape!=(len(tree.labels),4) or not np.isfinite(prob).all() or np.max(abs(prob.sum(axis=1)-1))>1e-8:
        raise ValueError('invalid reference ancestral probabilities')
    bounds=pd.read_csv(args.case/'tip_survival/multitree_probability_bounds.csv')
    if set(bounds['node'])!=set(range(1,len(tree.labels)+1)) or set(bounds['range'])!={'_','E','N','EN'}:
        raise ValueError('missing node/range source rows')
    if bounds.duplicated(['node','range']).any() or len(bounds)!=4*len(tree.labels):
        raise ValueError('duplicated or incomplete node/range source rows')
    if set(bounds['match_mode'])!={'monophyletic'}:raise ValueError('requires predeclared strict monophyletic matching')
    bounds['reference_probability']=[prob[int(n)-1,['_','E','N','EN'].index(r)] for n,r in zip(bounds.node,bounds.range)]
    bounds.to_csv(args.output/'ancestral_source.csv',index=False)
    configure_plotting()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.offsetbox import AnnotationBbox,DrawingArea
    from matplotlib.patches import Wedge
    fig=plt.figure(figsize=(180/25.4,245/25.4))
    a=fig.add_axes([.025,.07,.59,.875])
    b,c,d,e=[fig.add_axes([.765,y,.215,.155]) for y in [.775,.55,.325,.10]]
    range_colors=['#FFFFFF','#257C87','#D49A3B','#79669F']
    ys={};leaf_order=[]
    def locate(node):
        if not tree.children[node]:
            ys[node]=len(leaf_order);leaf_order.append(node)
        else:
            for child in tree.children[node]:locate(child)
            ys[node]=float(np.mean([ys[child] for child in tree.children[node]]))
    locate(tree.root)
    xs=tree.ages[tree.root]-tree.ages;depth=float(max(xs));width=depth if depth>0 else 1
    for node,children in enumerate(tree.children):
        if children:
            a.plot([xs[node],xs[node]],[min(ys[ch] for ch in children),max(ys[ch] for ch in children)],color='#626262',linewidth=.5,zorder=1)
            for ch in children:a.plot([xs[node],xs[ch]],[ys[ch],ys[ch]],color='#626262',linewidth=.5,zorder=1)
            disk=DrawingArea(4.6,4.6,0,0);angle=0
            for value,color in zip(prob[node],range_colors):
                if value>0:disk.add_artist(Wedge((2.3,2.3),2.3,angle,angle+360*value,facecolor=color,edgecolor='none'))
                angle+=360*value
            a.add_artist(AnnotationBbox(disk,(xs[node],ys[node]),frameon=False,pad=0))
        else:
            a.text(width*1.035,ys[node],tree.labels[node].replace('_',' '),va='center',ha='left',fontsize=8,fontstyle='italic')
            a.scatter(width*1.012,ys[node],color=range_colors[int(np.argmax(prob[node]))],marker='s',s=10,linewidths=0)
    a.set_xlim(-.025*width,3.35*width);a.set_ylim(tree.tip_count-.25,-1.25)
    a.set_axis_off();a.text(0,1.02,'a  Reference tree',transform=a.transAxes,fontsize=10,fontweight='bold')
    scale=width/4;a.plot([0,scale],[tree.tip_count+.1]*2,color='black',lw=.8,clip_on=False)
    a.text(scale/2,tree.tip_count+.9,f'{scale:.3g} tree-time units',ha='center',va='top',fontsize=8)
    a.legend(handles=[Line2D([],[],marker='o',color=range_colors[i],linestyle='',label=name) for i,name in enumerate(['_','E','N','EN']) if i>0],
        loc='lower left',bbox_to_anchor=(0,-.062),frameon=False,ncol=3,handletextpad=.3,columnspacing=.6)
    plus=['DEC+J','DIVALIKE+J','BAYAREALIKE+J'];short=['DEC+J','DIVA+J','BAY+J']
    missing_profiles=0
    for ci,contract in enumerate(contracts):
        for mi,model in enumerate(plus):
            chosen=profiles[(profiles.contract==contract)&(profiles.model==model)&(profiles.param=='j')]
            if len(chosen)!=1:raise ValueError('missing j profile record')
            row=chosen.iloc[0];y=mi+(ci-.5)*.22
            if not np.isfinite([row.lower,row.upper,row.mle]).all():
                missing_profiles+=1;continue
            b.plot([row.lower,row.upper],[y,y],color=colors[ci],linewidth=1)
            b.scatter(row.mle,y,s=14,marker=['o','s'][ci],color=colors[ci])
            if 'censor' in str(row.status_upper):b.scatter(row.upper,y,marker='>',s=20,color=colors[ci])
            if 'censor' in str(row.status_lower):b.scatter(row.lower,y,marker='<',s=20,color=colors[ci])
        for mi,row in enumerate(result['contracts'][contract]['bootstrap']):
            y=plus.index(row['alt'])+(ci-.5)*.22
            if row['mc_ci95'][0] is not None:c.plot(row['mc_ci95'],[y,y],color=colors[ci],linewidth=.7)
            c.plot(row['failure_bounds'],[y,y],color=colors[ci],linewidth=3)
            c.scatter(row['pval'],y,color=colors[ci],s=14,marker=['o','s'][ci],zorder=3)
    for ax in [b,c]:ax.set_yticks(range(3),short);ax.set_ylim(2.5,-.5)
    b.set_xlabel('j profile endpoints');b.set_xlim(left=-.002)
    b.text(0,-.38,f'{missing_profiles}/6 unavailable',transform=b.transAxes,fontsize=8)
    c.set_xlabel('Bootstrap p');c.set_xlim(-.04,1.04);c.set_xticks([0,.5,1])
    c.axvline(.05,color='black',linestyle=':',linewidth=.6)
    handle=[Line2D([],[],color=colors[i],marker=['o','s'][i],label=name,linewidth=1) for i,name in enumerate(['Conditioned','Legacy score'])]
    b.legend(handles=handle,frameon=False,loc='lower left',bbox_to_anchor=(-.2,1.03),handlelength=1.4,labelspacing=.25)
    internal=bounds[(bounds.node>tree.tip_count)&(bounds['range']!='_')]
    unavailable=~np.isfinite(internal.complete_case_probability)
    for ri,name in enumerate(['E','N','EN'],1):
        selected=internal[(internal['range']==name)&~unavailable]
        for _,row in selected.iterrows():
            d.plot([row.reference_probability]*2,[row.failure_lower,row.failure_upper],color=range_colors[ri],linewidth=.4,alpha=.35)
            d.scatter(row.reference_probability,row.complete_case_probability,s=9,
                facecolor=range_colors[ri] if row.clade_support_all_trees>=.95 else 'white',edgecolor=range_colors[ri],linewidth=.45)
    d.plot([0,1],[0,1],color='#666666',linewidth=.5,linestyle=':');d.set_xlim(-.04,1.04);d.set_ylim(-.04,1.04)
    d.set_xticks([0,.5,1]);d.set_yticks([0,.5,1]);d.set_xlabel('Reference probability');d.set_ylabel('Tree-average probability')
    d.text(0,-.39,f'{int(unavailable.sum())}/{len(internal)} unavailable',transform=d.transAxes,fontsize=8)
    for ci,contract in enumerate(contracts):
        for ei,event_type in enumerate(['d','e','j']):
            selected=events[(events.contract==contract)&(events.event_type==event_type)]
            if len(selected)!=1:raise ValueError('missing predeclared event type')
            row=selected.iloc[0];y=ei+(ci-.5)*.22
            if not np.isfinite(row['mean']):
                e.text(.02,y,'NA',transform=e.get_yaxis_transform(),color=colors[ci],fontsize=8)
                continue
            e.plot([row.rpi_025,row.rpi_975],[y,y],color=colors[ci],linewidth=.7)
            if np.isfinite([row.mean_mc_ci_lower,row.mean_mc_ci_upper]).all():
                e.plot([row.mean_mc_ci_lower,row.mean_mc_ci_upper],[y,y],color=colors[ci],linewidth=3)
            e.scatter(row['mean'],y,s=14,marker=['o','s'][ci],facecolor=colors[ci] if row.assessed else 'white',edgecolor=colors[ci],linewidth=.7)
    e.set_yticks(range(3),['Dispersal','Extirpation','Founder']);e.set_ylim(2.5,-.5);e.set_xlabel('Event count')
    for letter,ax,title in [('b',b,'Profile'),('c',c,'Nested comparison'),('d',d,'Across trees'),('e',e,'History vs precision')]:
        ax.text(-.55,1.08,letter,transform=ax.transAxes,fontsize=10,fontweight='bold')
        if letter!='b':ax.set_title(title,fontsize=9,pad=10)
    audit_alignment(fig,args.output,'figure3',args.audit_scripts)
    save_figure(fig,args.output,'figure3');plt.close(fig)
    (args.output/'caption.md').write_text('Figure 3. Parameter, tree-sample and mapping uncertainty in the locked Trimezieae example. (a) All 57 labelled tips and reference-node probability pies, averaged across six models using within-contract AICc weights. These weights are not posterior model probabilities. (b) All three +J profiles: conditioned asymptotic confidence intervals versus legacy-score support intervals; arrows denote computational-box censoring. (c) All nested comparisons; thin lines show Monte Carlo binomial intervals and thick segments show bounds from failed replicates. The dotted line marks 0.05. (d) Every nonempty-range probability for every internal reference clade, compared with strict monophyletic averages over the deposited tree sample. Open symbols have clade support below 0.95; vertical lines retain failure-mass bounds. Missing comparisons are counted, never filled with zeros. Colors follow panel a. (e) Dispersal, extirpation and founder-event counts; thin lines are empirical 2.5–97.5% history quantiles, thick lines are numerical intervals for the mean, and open means are not assessed for the precision target. The complete source table retains all eight event types and route comparisons. Deposited trees are treated as a fixed sample; this analysis does not establish population convergence or a biological causal claim.\n')
    write_json(args.output/'script_binding.json',dict(sha256=sha(__file__),helper_sha256=sha(Path(__file__).with_name('v12_figure_utils.py'))))


if __name__=='__main__':main()
