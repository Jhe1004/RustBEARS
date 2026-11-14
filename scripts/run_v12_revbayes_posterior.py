#!/usr/bin/env python3
"""Run two real RevBayes molecular chains and audit their sampled posterior.

Requires ArviZ 0.22.0 (a validation dependency, not a RustBEARS dependency).
The model uses JC69, a fixed known molecular clock, and a constant-population
coalescent with population ~ Exponential(rate=0.1). This prior differs from
the root-rescaled simulation distribution; it is explicitly part of the
external phylogenetic analysis, never a claim of prior-matched calibration.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears.uncertainty.averaging import node_clades
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
from scripts.validate_v12_v3_external_posterior import _load_dataset,newick,sha,write_json


def read_trace(path):
    return pd.read_csv(path,sep='\t',comment='#')


def diagnostics(trace_arrays,clade_arrays):
    import arviz as az
    rows=[];constants=[];constant_numeric=[]
    for kind,arrays in [('numeric',trace_arrays),('clade_indicator',clade_arrays)]:
        for key,values in arrays.items():
            values=np.asarray(values,dtype=float)
            if not np.isfinite(values).all():raise ValueError(f'nonfinite trace: {key}')
            if np.ptp(values)==0:
                constants.append(dict(kind=kind,name=key,value=float(values.flat[0])))
                if kind=='numeric':constant_numeric.append(key)
                continue
            rhat=float(az.rhat(values,method='rank'))
            bulk=float(az.ess(values,method='bulk'));tail=float(az.ess(values,method='tail'))
            if not np.isfinite([rhat,bulk,tail]).all():
                rows.append(dict(kind=kind,name=key,rhat=None,ess_bulk=None,ess_tail=None));continue
            rows.append(dict(kind=kind,name=key,rhat=rhat,ess_bulk=bulk,ess_tail=tail))
    finite=all(r['rhat'] is not None for r in rows)
    min_ess=min(min(r['ess_bulk'],r['ess_tail']) for r in rows) if finite and rows else 0.
    max_rhat=max(r['rhat'] for r in rows) if finite and rows else None
    return dict(chains=2,arviz_version=az.__version__,method='rank-normalized split R-hat; bulk and tail ESS',
        min_ess=min_ess,max_rhat=max_rhat,metrics=rows,constant_traces=constants,
        constant_numeric_traces=constant_numeric,
        diagnostic_pass=finite and bool(rows) and not constant_numeric and min_ess>=200 and max_rhat<=1.01,
        limitations='Numeric and observed clade-indicator diagnostics cannot prove exploration of unvisited topology modes. Constant observed clades carry no mixing evidence.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',type=Path,required=True)
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--generations',type=int,default=40000)
    parser.add_argument('--warmup',type=int,default=5000)
    parser.add_argument('--sample-every',type=int,default=20)
    parser.add_argument('--burnin-fraction',type=float,default=.25)
    parser.add_argument('--seed',type=int,default=20261115)
    args=parser.parse_args()
    if min(args.generations,args.warmup,args.sample_every)<1 or not 0<=args.burnin_fraction<1:
        raise ValueError('invalid MCMC lengths or burnin')
    folder,simulation=_load_dataset(args.dataset.resolve())
    executable=args.executable.resolve()
    if not executable.is_file():raise FileNotFoundError(executable)
    template=ROOT/'scripts/v3_jc69_coalescent.Rev.template'
    record_path=folder/'external_run.json'
    if record_path.exists():raise FileExistsError('External run already recorded; use a new prepared dataset directory.')
    tasks=[TaskIdentity('v3_molecular_mcmc','JC69',str(simulation['dataset']),outer_rep=i,role='chain') for i in (1,2)]
    seeds=derive_task_seeds(args.seed,tasks)
    # RevBayes accepts a 32-bit seed. Record and check this second truncation.
    seeds32=[seeds[t]%(2**31-1)+1 for t in tasks]
    if seeds32[0]==seeds32[1]:raise ValueError('RevBayes chain-seed collision')
    record=dict(status='running',engine='RevBayes',executable_sha256=sha(executable),
        template_sha256=sha(template),runner_sha256=sha(__file__),started=time.time(),
        sequence_sha256=sha(folder/'alignment.nex'),commands=[],outputs={},
        engine_config=dict(engine='RevBayes',branch_length_unit='time',clock_rate=simulation['clock_rate'],
            fixed_clock=True,geography_used=False,true_tree_used=False,model='JC69; constant coalescent; population Exp(rate 0.1)',
            posterior_already_burned_in=True,chains=2,generations=args.generations,warmup=args.warmup,
            sample_every=args.sample_every,burnin_fraction_per_chain=args.burnin_fraction))
    write_json(record_path,record)
    numeric=[];trees_by_chain=[];raw_outputs=[]
    for index,seed in enumerate(seeds32,1):
        prefix=folder/f'chain_{index}'
        tokens={'@SEED32@':str(seed),'@ALIGNMENT@':str(folder/'alignment.nex'),'@PREFIX@':str(prefix),
            '@CLOCK_RATE@':format(simulation['clock_rate'],'.17g'),'@GENERATIONS@':str(args.generations),
            '@WARMUP@':str(args.warmup),'@SAMPLE_EVERY@':str(args.sample_every)}
        script=template.read_text()
        for key,value in tokens.items():script=script.replace(key,value)
        script_path=folder/f'chain_{index}.Rev';script_path.write_text(script)
        command=[str(executable),'-o','outputPrecision=17',str(script_path)]
        record['commands'].append(dict(argv=command,seed32=seed,seed64=f'{seeds[tasks[index-1]]:016x}',input_sha256=sha(script_path)))
        write_json(record_path,record)
        with (folder/f'chain_{index}.stdout.log').open('w') as out,(folder/f'chain_{index}.stderr.log').open('w') as err:
            result=subprocess.run(command,cwd=folder,stdout=out,stderr=err,check=False)
        if result.returncode:
            record.update(status='failed',returncode=result.returncode,finished=time.time());write_json(record_path,record)
            raise RuntimeError(f'RevBayes chain {index} failed')
        log=read_trace(prefix.with_suffix('.log'))
        tree_path=prefix.with_suffix('.trees.jsonl')
        records=[json.loads(line) for line in tree_path.read_text().splitlines() if line.strip()]
        tree_trace=pd.DataFrame([row for row in records if 'Iteration' in row])
        if len(log)!=len(tree_trace) or not np.array_equal(log.iloc[:,0],tree_trace['Iteration']):
            raise ValueError('tree and parameter sample indices differ')
        start=int(len(log)*args.burnin_fraction)
        numeric.append(log.iloc[start:].drop(columns=log.columns[0]).reset_index(drop=True))
        trees=[b.read_newick(value) for value in tree_trace['phylogeny'].iloc[start:].astype(str)]
        trees_by_chain.append(trees)
        raw_outputs.extend([prefix.with_suffix('.log'),tree_path,script_path])
    if len(trees_by_chain[0])!=len(trees_by_chain[1]):raise ValueError('chain lengths differ')
    clades_by_chain=[[set(frozenset(c) for c in node_clades(tree) if 1<len(c)<tree.tip_count)
                     for tree in trees] for trees in trees_by_chain]
    all_clades=set().union(*(c for chain in clades_by_chain for c in chain))
    clade_arrays={','.join(sorted(c)):np.array([[float(c in sample) for sample in chain] for chain in clades_by_chain]) for c in all_clades}
    columns=set(numeric[0].columns)&set(numeric[1].columns)
    arrays={key:np.array([frame[key].to_numpy(dtype=float) for frame in numeric]) for key in sorted(columns)}
    diag=diagnostics(arrays,clade_arrays)
    diag.update(post_burnin_draws_per_chain=len(trees_by_chain[0]),burnin_fraction_per_chain=args.burnin_fraction,
        tree_serialization='RevBayes JSON tree serialization (15 decimal places in v1.4.2)',
        max_tip_height_roundoff=max(float(np.ptp(tree.ages[:tree.tip_count])) for chain in trees_by_chain for tree in chain),
        raw_files={str(p.name):sha(p) for p in raw_outputs})
    write_json(folder/'chain_diagnostics.json',diag)
    # Choose the sampled tree that maximizes the product of clade frequencies.
    # Its original branch lengths are retained; this is not a mean-age MCC.
    pooled=[tree for pair in zip(*trees_by_chain) for tree in pair]
    pooled_clades=[clades for pair in zip(*clades_by_chain) for clades in pair]
    counts=Counter(c for clades in pooled_clades for c in clades)
    scores=[sum(np.log(counts[c]/len(pooled)) for c in clades) for clades in pooled_clades]
    chosen=int(np.argmax(scores))
    (folder/'posterior_pooled.newick').write_text('\n'.join(newick(tree) for tree in pooled)+'\n')
    (folder/'mcc_sampled.newick').write_text(newick(pooled[chosen])+'\n')
    write_json(folder/'mcc_definition.json',dict(method='sampled maximum clade credibility tree; original sampled node ages',
        sample_index=chosen,log_clade_credibility=float(scores[chosen]),n_samples=len(pooled),
        pooling='interleaved chains after per-chain burnin; no further burnin required'))
    for key,name in [('posterior_trees','posterior_pooled.newick'),('mcc_tree','mcc_sampled.newick'),('diagnostics_file','chain_diagnostics.json')]:
        path=folder/name;record['outputs'][key]=dict(path=str(path),sha256=sha(path))
    record.update(status='completed',diagnostic_pass=diag['diagnostic_pass'],finished=time.time())
    write_json(record_path,record)
    print(json.dumps(dict(status='actual_two_chain_run_completed',diagnostic_pass=diag['diagnostic_pass'],
        min_ess=diag['min_ess'],max_rhat=diag['max_rhat'],posterior_samples=len(pooled))))


if __name__=='__main__':main()
