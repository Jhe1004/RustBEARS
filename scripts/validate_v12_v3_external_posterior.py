#!/usr/bin/env python3
"""V3 sequence simulation and an audited external Bayesian posterior interface.

prepare generates true coalescent trees, JC69 alignments and native geographic
histories. run requires a real executable plus a separately validated Bayesian
input template. collect accepts only successful recorded external runs. No
posterior-tree perturbation, bootstrap tree, or true-tree copy is substituted.

Engine JSON fields: executable, engine, argv (argument list), input_template,
main_input, posterior_trees, mcc_tree, branch_length_unit ('time' or
'substitutions_per_site'), clock_rate, fixed_clock, optional postprocess_argv,
diagnostics_file. Templates/argv use @ALIGNMENT@, @FASTA@, @PREFIX@, @INPUT@,
@SEED32@, @THREADS@, @CLOCK_RATE@, @OUTPUT_DIR@ tokens. Posterior/MCC paths are
relative to the dataset directory. Diagnostics are external, explicitly
labelled as such, and must state chains>=2, min_ess>=200 and max_rhat<=1.01.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears.uncertainty.averaging import node_clades
from rustbears.uncertainty.multitree import read_tree_sample,fit_over_trees
from rustbears.uncertainty.rng import TaskIdentity,derive_task_seeds
from rustbears.uncertainty.checkpoint import canonical_manifest,task_fingerprint,model_manifest


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n")


def newick(tree):
    def visit(node):
        children=tree.children[node]
        label="("+",".join(visit(c) for c in children)+")" if children else "'"+tree.labels[node].replace("'","''")+"'"
        return label+(":"+format(float(tree.lengths[node]),".17g") if node!=tree.root else "")
    return visit(tree.root)+";"


def coalescent_tree(tips,rng,root_age=10):
    lineages=[(f"T{i+1}",0.) for i in range(tips)];age=0.
    while len(lineages)>1:
        k=len(lineages);age+=rng.exponential(2/(k*(k-1)))
        indices=sorted(rng.choice(k,2,replace=False),reverse=True)
        left,right=[lineages.pop(int(i)) for i in indices]
        lineages.append((f"({left[0]}:{age-left[1]:.17g},{right[0]}:{age-right[1]:.17g})",age))
    tree=b.read_newick(lineages[0][0]+";")
    scale=root_age/age;tree.lengths*=scale;tree.ages*=scale
    return tree


def jc69_sequences(tree,sites,clock_rate,rng):
    """Exact JC69 transition draws, with branch lengths in simulation time."""
    sequences=np.zeros((len(tree.labels),sites),dtype=np.uint8)
    sequences[tree.root]=rng.integers(0,4,size=sites)
    stack=[tree.root]
    while stack:
        parent=stack.pop()
        for child in tree.children[parent]:
            value=sequences[parent].copy()
            probability=.75*(-np.expm1(-4*clock_rate*tree.lengths[child]/3))
            selected=rng.random(sites)<probability
            value[selected]=(value[selected]+rng.integers(1,4,size=int(selected.sum())))%4
            sequences[child]=value;stack.append(child)
    alphabet=np.array(list("ACGT"))
    return {name:"".join(alphabet[sequences[i]]) for i,name in enumerate(tree.tip_labels)}


def prepare(args):
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    if min(args.replicates,args.tips,args.areas,args.sites)<1 or args.tips<3 or args.areas>8 or args.clock_rate<=0:
        raise ValueError("invalid simulation dimensions or clock rate")
    budget_record=None
    if getattr(args,"mode","pilot")=="formal":
        if args.budget_json is None:raise ValueError("formal V3 preparation requires a locked budget")
        from scripts.v12_formal_protocol import load_locked_budget,select_budget_job
        budget,budget_record=load_locked_budget(args.budget_json,source_script=__file__,execute=True)
        job=select_budget_job(budget,analysis="v3")
        if any(job[key]!=getattr(args,key) for key in ("replicates","tips","areas","sites","clock_rate","seed")):
            raise ValueError("V3 simulation dimensions differ from the locked budget")
    elif getattr(args,"budget_json",None) is not None:
        raise ValueError("a formal budget cannot relabel a pilot")
    tasks=[TaskIdentity("v3_simulation","DEC+J",str(i),role=role) for i in range(1,args.replicates+1) for role in ("tree","sequence","geography")]
    seeds=derive_task_seeds(args.seed,tasks)
    rows=[]
    for i in range(1,args.replicates+1):
        folder=out/f"dataset_{i:03d}";folder.mkdir(exist_ok=True)
        if (folder/"simulation.json").exists():raise FileExistsError(f"refusing to overwrite {folder}")
        seed=lambda role:seeds[TaskIdentity("v3_simulation","DEC+J",str(i),role=role)]
        tree=coalescent_tree(args.tips,np.random.default_rng(seed("tree")))
        (folder/"true_tree.newick").write_text(newick(tree)+"\n")
        sequences=jc69_sequences(tree,args.sites,args.clock_rate,np.random.default_rng(seed("sequence")))
        (folder/"alignment.fasta").write_text("".join(f">{name}\n{seq}\n" for name,seq in sequences.items()))
        matrix="\n".join(f"{name} {seq}" for name,seq in sequences.items())
        (folder/"alignment.nex").write_text(f"#NEXUS\nBegin data;\nDimensions ntax={args.tips} nchar={args.sites};\nFormat datatype=dna missing=? gap=-;\nMatrix\n{matrix}\n;\nEnd;\n")
        areas=tuple(f"A{x+1}" for x in range(args.areas))
        dummy=b.TipRanges(areas,{name:(0,) for name in tree.tip_labels})
        model=b.BioGeoBEARS(tree,dummy,"DEC+J",likelihood_conditioning="tip_survival",min_branchlength=0)
        params=model.model.resolve(dict(d=.1,e=.05,j=.02))
        with b.PreparedSimulation(model,params) as sampler:
            batch=sampler.draw(1,seed=seed("geography"))
        geography=batch.to_tip_ranges(0)
        (folder/"geography.data").write_text(f"{args.tips} {args.areas} ({' '.join(areas)})\n"+"\n".join(
            name+" "+"".join("1" if a in geography.ranges[name] else "0" for a in range(args.areas)) for name in tree.tip_labels)+"\n")
        truth=[dict(node=node+1,taxa=sorted(clade),state=int(batch.node_states[0,node]),
                    range=[int(x) for x in model.states[int(batch.node_states[0,node])]]) for node,clade in enumerate(node_clades(tree))]
        write_json(folder/"true_states.json",truth)
        manifest=dict(dataset=i,sequence_model="JC69",sites=args.sites,clock_rate=args.clock_rate,locked_budget=budget_record,
                      branch_length_unit="simulation_time",tree_generation="Kingman coalescent conditioned on n tips; rescaled root age 10",
                      seeds={role:f"{seed(role):016x}" for role in ("tree","sequence","geography")},
                      model=model_manifest(model),generating_params=params,
                      files={name:sha(folder/name) for name in ("true_tree.newick","alignment.fasta","alignment.nex","geography.data","true_states.json")})
        result=dict(status="simulated_waiting_for_external_bayesian_inference",posterior_inference_completed=False,
                    task_fingerprint=task_fingerprint(manifest),manifest=canonical_manifest(manifest),
                    dataset=i,clock_rate=args.clock_rate,files=manifest["files"],tips=args.tips,areas=args.areas,
                    mode=getattr(args,"mode","pilot"),locked_budget=budget_record)
        write_json(folder/"simulation.json",result);rows.append(dict(dataset=i,path=str(folder),status=result["status"]))
    write_json(out/"study_plan.json",dict(protocol="V3-sequence-posterior-v1",formal_dimensions=dict(replicates=100,tips=50,areas=4,posterior_trees=200),
        requested_dimensions=dict(replicates=args.replicates,tips=args.tips,areas=args.areas,sites=args.sites),
        full_v3_completed=False,mode=getattr(args,"mode","pilot"),locked_budget=budget_record,datasets=rows,
        evaluation_rule="Score only internal clades present in both the true tree and MCC; report all excluded true/MCC clades. Bin paired scores by external clade support. No perturbed-tree substitute is permitted."))
    print(json.dumps(dict(status="prepared_external_inference_pending",datasets=len(rows),output=str(out))))


def _load_dataset(folder):
    folder=Path(folder);manifest=json.loads((folder/"simulation.json").read_text())
    for name,digest in manifest["files"].items():
        if sha(folder/name)!=digest:raise ValueError(f"simulated input changed: {name}")
    return folder,manifest


def run_engine(args):
    folder,simulation=_load_dataset(args.dataset)
    config_path=Path(args.engine_config).resolve();config=json.loads(config_path.read_text())
    required={"engine","executable","argv","input_template","main_input","posterior_trees","mcc_tree","branch_length_unit","clock_rate","fixed_clock","diagnostics_file"}
    if required-set(config):raise ValueError(f"engine configuration missing {sorted(required-set(config))}")
    executable=Path(config["executable"]).expanduser().resolve()
    if not executable.is_file():raise FileNotFoundError("a real Bayesian executable is required")
    for name in ("main_input","posterior_trees","mcc_tree","diagnostics_file"):
        value=Path(config[name])
        if value.is_absolute() or ".." in value.parts:
            raise ValueError(f"{name} must be a relative path confined to the dataset directory")
    if not isinstance(config["argv"],list) or not all(isinstance(x,str) for x in config["argv"]):raise ValueError("argv must be an argument list; shell commands are not accepted")
    if not config["fixed_clock"] or float(config["clock_rate"])!=simulation["clock_rate"]:
        raise ValueError("this V3 protocol requires the known generating clock rate fixed in the Bayesian model")
    tokens={"@ALIGNMENT@":str((folder/"alignment.nex").resolve()),"@FASTA@":str((folder/"alignment.fasta").resolve()),
            "@OUTPUT_DIR@":str(folder.resolve()),"@PREFIX@":str((folder/"posterior").resolve()),
            "@INPUT@":str((folder/config["main_input"]).resolve()),"@THREADS@":str(args.threads),
            "@SEED32@":str(args.seed),"@CLOCK_RATE@":str(config["clock_rate"])}
    def expand(text):
        for key,value in tokens.items():text=text.replace(key,value)
        return text
    template=Path(config["input_template"]).expanduser().resolve()
    input_path=folder/config["main_input"]
    input_path.write_text(expand(template.read_text()))
    commands=[[str(executable),*[expand(x) for x in config["argv"]]]]
    if config.get("postprocess_argv"):
        commands.append([expand(x) for x in config["postprocess_argv"]])
    record=dict(status="running",engine=config["engine"],engine_config=config,engine_config_sha256=sha(config_path),
                executable_sha256=sha(executable),template_sha256=sha(template),input_sha256=sha(input_path),
                sequence_sha256=sha(folder/"alignment.nex"),commands=commands,seed=args.seed,started=time.time(),outputs={})
    write_json(folder/"external_run.json",record)
    for index,command in enumerate(commands):
        if any("@" in token for token in command):raise ValueError("unexpanded engine argument token")
        with (folder/f"external_{index}.stdout.log").open("w") as out,(folder/f"external_{index}.stderr.log").open("w") as err:
            completed=subprocess.run(command,cwd=folder,stdout=out,stderr=err,check=False)
        if completed.returncode:
            record.update(status="failed",returncode=completed.returncode,finished=time.time());write_json(folder/"external_run.json",record)
            raise RuntimeError(f"external Bayesian workflow returned {completed.returncode}")
    for name in ("posterior_trees","mcc_tree","diagnostics_file"):
        path=folder/config[name]
        if not path.is_file():raise FileNotFoundError(f"external workflow did not produce {name}: {path}")
        record["outputs"][name]=dict(path=str(path.resolve()),sha256=sha(path))
    record.update(status="completed",finished=time.time());write_json(folder/"external_run.json",record)
    print(json.dumps(dict(status="external_workflow_completed_not_yet_imported",dataset=str(folder))))


def collect(args):
    folder,simulation=_load_dataset(args.dataset);run=json.loads((folder/"external_run.json").read_text())
    if run["status"]!="completed":raise ValueError("a completed recorded Bayesian run is required")
    if run.get('sequence_sha256')!=sha(folder/'alignment.nex'):
        raise ValueError('External posterior was not inferred from this dataset alignment.')
    config=run['engine_config']
    if config.get('fixed_clock') is not True or config.get('clock_rate')!=simulation['clock_rate']:
        raise ValueError('External posterior clock does not match the fixed generating clock.')
    for item in run["outputs"].values():
        if sha(item["path"])!=item["sha256"]:raise ValueError("external posterior output changed after execution")
    diagnostics=json.loads(Path(run["outputs"]["diagnostics_file"]["path"]).read_text())
    if diagnostics.get('diagnostic_pass') is False:
        raise ValueError('External diagnostic workflow explicitly failed a convergence check.')
    ess=diagnostics.get('min_ess');rhat=diagnostics.get('max_rhat')
    if (diagnostics.get("chains",0)<2 or ess is None or rhat is None
            or not np.isfinite([ess,rhat]).all() or ess<200 or rhat>1.01):
        raise ValueError("external chain diagnostics have not met the preregistered convergence thresholds")
    if run['engine_config'].get('posterior_already_burned_in') and args.burnin_fraction!=0:
        raise ValueError('This posterior already removed per-chain burnin; use --burnin-fraction 0.')
    trees=read_tree_sample(Path(run["outputs"]["posterior_trees"]["path"]))
    burn=int(np.floor(len(trees)*args.burnin_fraction));remaining=trees[burn:]
    if len(remaining)<args.n_trees:raise ValueError("insufficient actual post-burnin posterior trees")
    indices=np.linspace(0,len(remaining)-1,args.n_trees,dtype=int);selected=[remaining[int(i)] for i in indices]
    mcc=b.read_newick(Path(run["outputs"]["mcc_tree"]["path"]))
    if config["branch_length_unit"]=="substitutions_per_site":scale=1/simulation["clock_rate"]
    elif config["branch_length_unit"]=="time":scale=1.
    else:raise ValueError("posterior branch_length_unit must be time or substitutions_per_site")
    truth=b.read_newick(folder/"true_tree.newick")
    for tree in selected+[mcc]:
        if set(tree.tip_labels)!=set(truth.tip_labels):raise ValueError("posterior taxa differ from the simulated alignment")
        tree.lengths*=scale;tree.ages*=scale;tree.root_edge*=scale
    (folder/"posterior_selected.newick").write_text("\n".join(newick(tree) for tree in selected)+"\n")
    (folder/"mcc_time_units.newick").write_text(newick(mcc)+"\n")
    result=dict(status="actual_external_posterior_imported",source_run_sha256=sha(folder/"external_run.json"),
        source_tree_indices=(indices+burn).tolist(),n_trees=len(selected),burnin_fraction=args.burnin_fraction,
        unit_conversion_multiplier=scale,convergence_evidence=diagnostics,
        convergence_evidence_scope="Reported by the configured external diagnostic workflow; RustBEARS does not silently infer convergence from tree count.",
        files={name:sha(folder/name) for name in ("posterior_selected.newick","mcc_time_units.newick")})
    write_json(folder/"posterior_import.json",result);print(json.dumps(dict(status=result["status"],n_trees=len(selected))))


def evaluate(args):
    folder,_simulation=_load_dataset(args.dataset);imported=json.loads((folder/"posterior_import.json").read_text())
    if imported["status"]!="actual_external_posterior_imported":raise ValueError("genuine imported posterior is required")
    for name,digest in imported["files"].items():
        if sha(folder/name)!=digest:raise ValueError("imported posterior changed")
    geography=b.read_geog(folder/"geography.data");mcc=b.read_newick(folder/"mcc_time_units.newick")
    single=b.BioGeoBEARS(mcc,geography,"DEC+J",likelihood_conditioning="tip_survival").fit(num_threads=args.threads)
    if not single.converged:raise RuntimeError("MCC fit failed")
    ensemble=fit_over_trees(folder/"posterior_selected.newick",geography,"DEC+J",summary_tree=mcc,
        likelihood_conditioning="tip_survival",criterion='AIC',init_strategy="full_cold_start",num_threads=args.threads,
        checkpoint_file=str(folder/"biogeography.jsonl"),resume=args.resume)
    true_states=json.loads((folder/"true_states.json").read_text());lookup={frozenset(r["taxa"]):r for r in true_states if len(r["taxa"])>1}
    paired=[];excluded=[]
    for node,clade in enumerate(node_clades(mcc)):
        if len(clade)<2:continue
        if clade not in lookup:
            excluded.append(sorted(clade));continue
        true_state=lookup[clade]["state"];target=np.zeros(len(single.model.states));target[true_state]=1
        support=float(np.clip(ensemble.clade_support_all_trees[node],0.,1.));row=dict(node=node+1,taxa=sorted(clade),clade_support=support,
            support_bin="<0.7" if support<.7 else "0.7-0.95" if support<=.95 else ">0.95",truth_state=int(true_state),
            complete_trees=int(ensemble.clade_complete_trees_count[node]),
            multitree_failure_bounds=[bound[node].tolist() for bound in ensemble.prob_failure_bounds_top])
        for name,probability in (("mcc",single.evaluation.ancestral_probs_branch_top[node]),("multitree",ensemble.complete_case_probs_top[node])):
            if not np.isfinite(probability).all():raise RuntimeError("undefined complete-case probability in paired score")
            p=float(probability[true_state]);row[name]=dict(brier=float(np.sum((probability-target)**2)),
                log_score=float(-np.log(p)) if p>0 else None,zero_probability=(p==0),probabilities=probability.tolist())
        paired.append(row)
    write_json(folder/"v3_scores.json",dict(status="evaluated_actual_external_posterior",posterior_import_sha256=sha(folder/"posterior_import.json"),
        scope='One dataset; formal V3 completion requires every preregistered dataset and aggregate calibration.',
        dataset_tips=mcc.tip_count,areas=len(geography.areas),posterior_trees=imported['n_trees'],
        weighting_criterion='AIC with one candidate model; no model selection or AICc correction',
        n_trees_planned=ensemble.n_trees_planned,n_trees_completed=ensemble.n_trees_completed,n_trees_incomplete=ensemble.n_trees_incomplete,
        paired_clades=paired,excluded_mcc_clades_not_in_truth=excluded,
        true_internal_clades_not_in_mcc=[sorted(c) for c in lookup if c not in set(node_clades(mcc))],
        scoring_boundary="Strict shared true/MCC internal clades only; exclusions are reported rather than assigning a false ancestral ground truth."))
    print(json.dumps(dict(status="evaluated",paired_clades=len(paired),excluded_mcc_clades=len(excluded))))


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="action",required=True)
    p=sub.add_parser("prepare");p.add_argument("--output-dir",required=True);p.add_argument("--replicates",type=int,default=100)
    p.add_argument("--tips",type=int,default=50);p.add_argument("--areas",type=int,default=4);p.add_argument("--sites",type=int,default=2000)
    p.add_argument("--clock-rate",type=float,default=.01);p.add_argument("--seed",type=int,default=20261115)
    p.add_argument("--mode",choices=("pilot","formal"),default="pilot");p.add_argument("--budget-json",type=Path)
    p=sub.add_parser("run");p.add_argument("--dataset",required=True);p.add_argument("--engine-config",required=True);p.add_argument("--threads",type=int,default=1);p.add_argument("--seed",type=int,default=20261115)
    p=sub.add_parser("collect");p.add_argument("--dataset",required=True);p.add_argument("--burnin-fraction",type=float,default=.25);p.add_argument("--n-trees",type=int,default=200)
    p=sub.add_parser("evaluate");p.add_argument("--dataset",required=True);p.add_argument("--threads",type=int,default=1);p.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    if args.action=="collect" and (not 0<=args.burnin_fraction<1 or args.n_trees<1):parser.error("invalid posterior selection settings")
    {"prepare":prepare,"run":run_engine,"collect":collect,"evaluate":evaluate}[args.action](args)


if __name__=="__main__":main()
