"""Conditional fits on an ordered external tree sample, with missing-fit bounds."""
import copy
from dataclasses import dataclass, replace
from pathlib import Path
import re
import time
import warnings

import numpy as np
import pandas as pd

from ..io import PhyloTree, read_newick, read_geog
from ..model import BioGeoBEARS, BioGeoBEARSModelObject, FitResult, model_template
from .averaging import node_clades, topology_hash, model_average, allocate_maps, sample_strata
from .checkpoint import CheckpointStore, model_manifest, frame_arrays, restore_frame
from .diagnostics import _make_bsm
from .parallel import resolve_parallel, parallel_scope
from .rng import TaskIdentity, derive_task_seeds, task_fingerprint as rng_fingerprint


class TreeFailureError(RuntimeError):
    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def _statements(text):
    """Split Newick/Nexus statements while respecting quotes and comments."""
    start = 0
    quote = False
    comment = 0
    i = 0
    while i < len(text):
        c = text[i]
        if c == "'" and not comment:
            if quote and i+1 < len(text) and text[i+1] == "'":
                i += 2
                continue
            quote = not quote
        elif not quote:
            if c == "[":
                comment += 1
            elif c == "]":
                comment -= 1
            elif c == ";" and not comment:
                yield text[start:i+1].strip()
                start = i+1
        i += 1
    if quote or comment or text[start:].strip():
        raise ValueError("unterminated tree statement, comment or quote")


def read_tree_sample(source):
    """Read plain Newick collections or Nexus trees with a Translate table."""
    if isinstance(source, PhyloTree):
        return [copy.deepcopy(source)]
    if isinstance(source, (str, Path)):
        raw = str(source)
        is_path = isinstance(source, Path) or not any(c in raw for c in "();\n")
        if not is_path and "\n" not in raw and len(raw)<4096:
            try:
                is_path = Path(raw).is_file()
            except OSError:
                is_path = False
        text = Path(source).read_text() if is_path else raw
        statements = list(_statements(text))
        translate = {}
        trees = []
        nexus = text.lstrip().upper().startswith("#NEXUS")
        for statement in statements:
            match = re.search(r"\btranslate\b(.*);\s*$", statement, re.I | re.S)
            if nexus and match:
                # Taxon names follow the standard quoted/unquoted Nexus tokens.
                token = r"('(?:[^']|'')*'|[^\s,]+)"
                for entry in re.finditer(r"(?:^|,)\s*(\S+)\s+"+token+r"\s*(?=,|$)", match[1].strip()):
                    label = entry[2]
                    translate[entry[1]] = label[1:-1].replace("''", "'") if label.startswith("'") else label
                if not translate:
                    raise ValueError("invalid Nexus Translate block")
            elif nexus:
                match = re.search(r"\b(?:tree|utree)\s+[^=]+=(.*);\s*$", statement, re.I | re.S)
                if not match:
                    continue
                tree = read_newick(match[1].strip()+";")
                if translate:
                    if not set(tree.tip_labels).issubset(translate):
                        raise ValueError("Nexus Translate lacks a tip label")
                    tree.labels = tuple(translate[label] for label in tree.tip_labels) + tree.labels[tree.tip_count:]
                    if len(set(tree.tip_labels)) != tree.tip_count:
                        raise ValueError("duplicate translated taxa")
                trees.append(tree)
            else:
                trees.append(read_newick(statement))
        return trees
    return [copy.deepcopy(read_newick(t)) for t in source]


def _mapping(tree, summary, match_mode):
    clades, wanted = node_clades(tree), node_clades(summary)
    lookup = {c: i for i, c in enumerate(clades)}
    indices, mono, excess = [], [], []
    for clade in wanted:
        present = clade in lookup
        mono.append(present)
        if present:
            node = lookup[clade]
        elif match_mode == "mrca":
            node = min((i for i,c in enumerate(clades) if clade.issubset(c)), key=lambda i: len(clades[i]))
        else:
            node = -1
        indices.append(node)
        excess.append(len(clades[node])-len(clade) if node >= 0 else 0)
    return np.asarray(indices), np.asarray(mono), np.asarray(excess)


@dataclass
class MultiTreeFitResult:
    summary_tree: object
    n_trees_planned: int
    n_trees_completed: int
    n_trees_incomplete: int
    models: tuple
    likelihood_conditioning: str
    precision_mode: str
    match_mode: str
    cold_start_audit_max_delta: float | None
    clade_support_all_trees: np.ndarray
    clade_complete_trees_count: np.ndarray
    mean_mrca_excess_tips: np.ndarray
    complete_case_probs_top: np.ndarray
    complete_case_probs_bottom: np.ndarray
    prob_failure_bounds_top: tuple
    prob_failure_bounds_bottom: tuple
    cross_tree_prob_sd_top: np.ndarray
    per_tree_df: pd.DataFrame
    model_winner_freq: pd.Series
    param_tree_sensitivity_df: pd.DataFrame
    joint_bsm: object
    elapsed_sec: float
    ranges_list: tuple
    states_list: tuple
    tree_order: tuple
    task_fingerprint: str | None
    fit_method: str
    interpretation: str = "Conditional ML/Akaike summaries over external trees; complete-case means and deterministic missing-fit bounds"

    def to_dataframe(self):
        rows = []
        low, high = self.prob_failure_bounds_top
        for node in range(len(self.summary_tree.labels)):
            for state, label in enumerate(self.ranges_list):
                rows.append(dict(node=node+1, range=label,
                    complete_case_probability=self.complete_case_probs_top[node,state],
                    failure_lower=low[node,state], failure_upper=high[node,state],
                    clade_support_all_trees=self.clade_support_all_trees[node],
                    clade_complete_trees_count=self.clade_complete_trees_count[node],
                    mean_mrca_excess_tips=self.mean_mrca_excess_tips[node], match_mode=self.match_mode))
        return pd.DataFrame(rows)

    def plot(self, filepath=None, min_clade_support=.5):
        import matplotlib.pyplot as plt
        if not 0 <= min_clade_support <= 1:
            raise ValueError("min_clade_support must be in [0,1]")
        data = self.complete_case_probs_top.copy()
        data[self.clade_support_all_trees < min_clade_support] = np.nan
        fig, (ax, support) = plt.subplots(1,2, gridspec_kw={"width_ratios":[4,1]}, sharey=True)
        image = ax.imshow(data, aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(len(self.ranges_list)), self.ranges_list, rotation=90)
        ax.set(title="Complete-case conditional probabilities", ylabel="Summary-tree node")
        support.barh(np.arange(len(data)),self.clade_support_all_trees)
        support.set(xlim=(0,1),title="Clade support")
        fig.colorbar(image, ax=ax)
        if filepath:
            fig.savefig(filepath,bbox_inches="tight")
        return fig


def fit_over_trees(trees, geog, models="DEC", summary_tree=None,
                  likelihood_conditioning="legacy_unconditioned", root_prior=None,
                  match_mode="monophyletic", init_strategy="warm_hybrid", cold_start_audit_fraction=.05,
                  model_averaging_within_tree=True, max_tree_failure_rate=.02,
                  on_excess_failure="abort", bsm_maps_per_tree=0, max_trees=None,
                  num_threads=8, seed=12345, checkpoint_file=None, *, resume=False,
                  criterion="AICc", tree_weights=None, parallel=None, fit_options=None,
                  cold_start_action="warn", **model_options):
    """Fit each complete candidate set; preserve original tree order and failures.

    Warm starts use a fixed reference-tree fit, making results independent of
    worker completion order. The most topologically distant trees also receive
    a cold fit. Set cold_start_action='refit_all' to rerun all trees cold when
    that audit finds an improvement greater than 1e-3.
    """
    started = time.monotonic()
    if resume and checkpoint_file is None:
        raise ValueError("resume requires checkpoint_file")
    sample = read_tree_sample(trees)
    if max_trees is not None:
        if isinstance(max_trees,bool) or not isinstance(max_trees,(int,np.integer)) or max_trees < 1:
            raise ValueError("max_trees must be a positive integer")
        sample = sample[:max_trees]
    if not sample:
        raise ValueError("tree sample is empty")
    geography = copy.deepcopy(read_geog(geog))
    taxa = set(geography.ranges)
    if any(set(t.tip_labels) != taxa for t in sample):
        raise ValueError("every tree must have exactly the geography taxon set")
    summary = copy.deepcopy(read_newick(summary_tree)) if summary_tree is not None else copy.deepcopy(sample[0])
    if set(summary.tip_labels) != taxa:
        raise ValueError("summary tree has a different taxon set")
    if match_mode not in ("monophyletic", "mrca") or init_strategy not in ("warm_hybrid", "full_cold_start", "fixed_ref_params"):
        raise ValueError("unknown clade matching or initialization strategy")
    if on_excess_failure not in ("abort", "conservative_bounds") or cold_start_action not in ("warn", "refit_all"):
        raise ValueError("unknown failure or cold-start action")
    if not 0 <= max_tree_failure_rate <= 1 or not 0 <= cold_start_audit_fraction <= 1:
        raise ValueError("invalid failure rate or audit fraction")
    if isinstance(bsm_maps_per_tree,bool) or not isinstance(bsm_maps_per_tree,(int,np.integer)) or bsm_maps_per_tree < 0:
        raise ValueError("bsm_maps_per_tree must be a nonnegative integer")
    definitions = [models] if isinstance(models,(str,BioGeoBEARSModelObject)) else list(models)
    definitions = tuple(model_template(x) if isinstance(x,str) else x for x in definitions)
    if not definitions or any(not isinstance(x,BioGeoBEARSModelObject) for x in definitions):
        raise ValueError("models must be a nonempty set of model definitions")
    k, g = len(sample), len(definitions)
    weights = np.ones(k)/k if tree_weights is None else np.asarray(tree_weights,float)
    if weights.shape != (k,) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("tree_weights must be finite, positive and match the tree sample")
    weights = weights/weights.sum()
    config = resolve_parallel(parallel,num_threads)
    options = dict(fit_options or {})
    if {"seed","num_threads","initial"}.intersection(options):
        raise ValueError("fit_options cannot override managed seed, threads or initial values")
    model_options = dict(model_options)
    model_options.setdefault("allow_null_tips",likelihood_conditioning=="unconditional_complete")

    def build(tree, definition):
        return BioGeoBEARS(copy.deepcopy(tree), copy.deepcopy(geography), definition,
                          likelihood_conditioning=likelihood_conditioning, root_prior=root_prior, **copy.deepcopy(model_options))

    templates = [build(sample[0],definition) for definition in definitions]
    # Fail unsupported models before fitting or opening a checkpoint.
    for model in templates:
        model._require_tier1a()
    if any(m.states != templates[0].states for m in templates):
        raise ValueError("candidate models have different state spaces")
    if criterion not in ("AIC", "AICc"):
        raise ValueError("criterion must be AIC or AICc")
    if criterion == "AICc" and any(sum(m.model.params[name].min < m.model.params[name].max
                                      for name in m.model.free)+1 >= m.tree.tip_count for m in templates):
        raise ValueError("AICc undefined for this taxon/parameter count; choose criterion='AIC' explicitly")
    if bsm_maps_per_tree and bsm_maps_per_tree < 2*(g if model_averaging_within_tree else 1):
        raise ValueError("BSM budget requires two maps per active model")
    manifest = dict(kind="fit_over_trees", trees=sample, summary_tree=summary, geography=geography,
                    models=[model_manifest(m) for m in templates], weights=weights, match_mode=match_mode,
                    init_strategy=init_strategy, cold_start_audit_fraction=cold_start_audit_fraction,
                    cold_start_action=cold_start_action, criterion=criterion,
                    model_averaging_within_tree=model_averaging_within_tree,
                    max_tree_failure_rate=max_tree_failure_rate, on_excess_failure=on_excess_failure,
                    bsm_maps_per_tree=bsm_maps_per_tree,max_trees=max_trees,seed=seed,fit_options=options)
    store = CheckpointStore(checkpoint_file,manifest,resume) if checkpoint_file else None
    # Derive every optimization key before dispatch; no scheduler-dependent seed state.
    tasks = [TaskIdentity("multitree_fit",str(j+1),str(i+1),role=role)
             for i in range(k) for j in range(g) for role in ("fit","retry","audit")]
    tasks += [TaskIdentity("multitree_fit",str(j+1),"reference",role="fit") for j in range(g)]
    if bsm_maps_per_tree:
        tasks += [TaskIdentity("multitree_bsm",str(j+1),str(i+1),role="draw") for i in range(k) for j in range(g)]
    seeds = derive_task_seeds(seed,tasks)
    def key(i,j,role="fit"):
        # NumPy/scipy optimization generators accept 64-bit integer seeds.
        return seeds[TaskIdentity("multitree_fit",str(j+1),str(i+1),role=role)]
    mappings = [_mapping(t,summary,match_mode) for t in sample]
    ref_clades = set(node_clades(summary))
    distances = [len(set(node_clades(t)).symmetric_difference(ref_clades)) for t in sample]
    audit_count = min(k,max(5,int(np.ceil(cold_start_audit_fraction*k)))) if init_strategy=="warm_hybrid" else 0
    audit_ids = set(sorted(range(k),key=lambda i:(-distances[i],i))[:audit_count])
    refs = None
    completed = {r["id"]: (r,store.load_arrays(r["id"])) for r in store.records} if store else {}
    if len(completed) < k and init_strategy != "full_cold_start":
        refs = []
        for j, definition in enumerate(definitions):
            try:
                fit = build(summary,definition).fit(num_threads=config.inner_threads,
                        seed=seeds[TaskIdentity("multitree_fit",str(j+1),"reference",role="fit")],**options)
                if not fit.converged or not np.isfinite(fit.lnL):
                    fit = None
            except (ValueError,RuntimeError,FloatingPointError):
                fit = None
            if fit is None:
                if init_strategy=="fixed_ref_params":
                    raise RuntimeError("fixed_ref_params requires converged reference-tree fits")
                warnings.warn(f"reference fit for {definition.name} failed; using cold starts for this model",RuntimeWarning)
            refs.append(fit)

    def worker(i):
        fits, rows, arrays, deltas = [], [], {}, []
        for j, definition in enumerate(definitions):
            model = build(sample[i],definition)
            fit, retried, error = None, False, None
            try:
                if init_strategy=="fixed_ref_params":
                    evaluation = model.evaluate(refs[j].params,num_threads=config.inner_threads)
                    nparams = refs[j].n_params
                    aic = 2*nparams-2*evaluation.lnL
                    aicc = aic+2*nparams*(nparams+1)/(model.tree.tip_count-nparams-1) if model.tree.tip_count>nparams+1 else np.inf
                    fit = FitResult(evaluation.lnL,dict(refs[j].params),nparams,aic,aicc,1,True,evaluation,model,"fixed reference parameters")
                else:
                    initial = refs[j].params if refs is not None and refs[j] is not None else None
                    try:
                        fit = model.fit(initial=initial,num_threads=config.inner_threads,seed=key(i,j),**options)
                    except (ValueError,RuntimeError,FloatingPointError):
                        fit = None
                    if fit is None or not fit.converged:
                        retried = True
                        retry = model.fit(initial=fit.params if fit is not None else initial,
                                          num_threads=config.inner_threads,seed=key(i,j,"retry"),**options)
                        if fit is None or retry.converged and (not fit.converged or retry.lnL >= fit.lnL):
                            fit = retry
                    if i in audit_ids:
                        cold = build(sample[i],definition).fit(num_threads=config.inner_threads,seed=key(i,j,"audit"),**options)
                        if cold.converged:
                            delta = cold.lnL-fit.lnL
                            deltas.append(float(delta))
                            if not fit.converged or delta>0:
                                fit=cold
                        else:
                            raise RuntimeError("cold-start audit did not converge")
                if not fit.converged or not np.isfinite(fit.lnL):
                    raise RuntimeError("optimizer did not return a finite converged fit")
            except (ValueError,RuntimeError,FloatingPointError) as exc:
                error = f"{type(exc).__name__}: {exc}"
            successful = error is None
            fits.append(fit if successful else None)
            row = dict(tree_id=i+1,tree_topology_hash=topology_hash(sample[i]),model_id=j+1,
                       model=definition.name,fit_status="converged" if successful else "failed",
                       retried=retried,error=error,params=fit.params if fit is not None else None,
                       lnL=fit.lnL if fit is not None else np.nan,
                       aic=fit.aic if fit is not None else np.nan,aicc=fit.aicc if fit is not None else np.nan,
                       weight=np.nan,map_id=None,stratum_weight=np.nan,
                       precision_mode=model.precision,likelihood_conditioning=likelihood_conditioning,
                       fit_method=init_strategy,
                       rng_task=TaskIdentity("multitree_fit",str(j+1),str(i+1),role="fit").to_dict(),
                       rng_seed_hex=f"{key(i,j):016x}",
                       rng_task_sha256=rng_fingerprint(seed,TaskIdentity("multitree_fit",str(j+1),str(i+1),role="fit")))
            if fit is not None:
                arrays[f"model{j+1}_top"] = fit.evaluation.ancestral_probs_branch_top
                arrays[f"model{j+1}_bottom"] = fit.evaluation.ancestral_probs_branch_bottom
            rows.append(row)
        record = dict(id=i+1,models=rows,tree_status="complete" if all(f is not None for f in fits) else "incomplete_model_set",
                      audit_delta=max(deltas) if deltas else None)
        if record["tree_status"]=="complete":
            averaged = model_average(fits,criterion)
            local_weights = averaged.weights_df.weight.to_numpy()
            if not model_averaging_within_tree:
                winner = int(np.argmax(local_weights))
                local_weights = np.eye(1,g,winner).ravel()
            for row,w in zip(rows,local_weights):
                row["weight"] = float(w)
                row["stratum_weight"] = float(weights[i]*w)
            arrays["top"] = sum(w*f.evaluation.ancestral_probs_branch_top for w,f in zip(local_weights,fits))
            arrays["bottom"] = sum(w*f.evaluation.ancestral_probs_branch_bottom for w,f in zip(local_weights,fits))
            if bsm_maps_per_tree:
                maps = sample_strata(fits,local_weights,allocate_maps(local_weights,bsm_maps_per_tree),
                                     seed=seed,tree_id=i+1,purpose="multitree_bsm",num_threads=config.inner_threads,task_seeds=seeds)
                for name, table in (("clado",maps.clado_events_df),("ana",maps.ana_events_df),("maps",maps.map_metadata)):
                    metadata, data = frame_arrays(table,name)
                    record[name] = metadata
                    arrays.update(data)
        for row in rows:
            row["tree_status"] = record["tree_status"]
        return record,arrays

    pending = [i for i in range(k) if i+1 not in completed]
    with parallel_scope(config) as dispatch:
        for record, arrays in dispatch(worker,pending):
            completed[record["id"]] = (record,arrays)
            if store:
                store.append(record,arrays)
    audits = [r["audit_delta"] for r,_ in completed.values() if r["audit_delta"] is not None]
    audit_max = max(audits) if audits else None
    if audit_max is not None and audit_max>1e-3:
        warnings.warn(f"cold-start audit improved a warm fit by {audit_max:.6g}; full cold refitting is recommended",RuntimeWarning)
        if cold_start_action=="refit_all":
            # A separate checkpoint avoids silently mixing two fitting protocols.
            return fit_over_trees(sample,geography,definitions,summary,likelihood_conditioning,root_prior,
                match_mode,"full_cold_start",cold_start_audit_fraction,model_averaging_within_tree,
                max_tree_failure_rate,on_excess_failure,bsm_maps_per_tree,None,num_threads,seed,
                None if checkpoint_file is None else str(checkpoint_file)+".cold",resume=resume,
                criterion=criterion,tree_weights=weights,parallel=config,fit_options=options,**model_options)
    ok = np.array([completed[i+1][0]["tree_status"]=="complete" for i in range(k)])
    nnode,nstate = len(summary.labels),len(templates[0].states)
    top,bottom = np.full((k,nnode,nstate),np.nan),np.full((k,nnode,nstate),np.nan)
    mono = np.array([mapping[1] for mapping in mappings])
    for i,(indices,_mono,_excess) in enumerate(mappings):
        if ok[i]:
            usable=indices>=0
            top[i,usable]=completed[i+1][1]["top"][indices[usable]]
            bottom[i,usable]=completed[i+1][1]["bottom"][indices[usable]]
    eligible = mono if match_mode=="monophyletic" else np.ones_like(mono)
    success = eligible & ok[:,None]
    denominator = weights @ eligible
    complete_mass = weights @ success

    def aggregate(values):
        numerator = np.sum(np.where(success[:,:,None],np.nan_to_num(values),0)*weights[:,None,None],axis=0)
        means=np.divide(numerator,complete_mass[:,None],out=np.full((nnode,nstate),np.nan),where=complete_mass[:,None]>0)
        lower=np.divide(numerator,denominator[:,None],out=np.full_like(means,np.nan),where=denominator[:,None]>0)
        upper=np.divide(numerator+(denominator-complete_mass)[:,None],denominator[:,None],out=np.full_like(means,np.nan),where=denominator[:,None]>0)
        # Root-branch bottom is undefined when no root edge exists, not a zero probability vector.
        undefined = np.all(~np.isfinite(values),axis=(0,2)) & (complete_mass>0)
        means[undefined]=lower[undefined]=upper[undefined]=np.nan
        return means,(lower,upper)

    mean_top,bounds_top=aggregate(top)
    mean_bottom,bounds_bottom=aggregate(bottom)
    # An MRCA can be the root in only some trees. Its incoming branch bottom
    # then has no definition; never turn those NaNs into zero probabilities.
    undefined_bottom=np.any(np.array([mapping[0]==tree.root for tree,mapping in zip(sample,mappings)]) & eligible,axis=0)
    mean_bottom[undefined_bottom]=np.nan
    bounds_bottom[0][undefined_bottom]=bounds_bottom[1][undefined_bottom]=np.nan
    sd=np.sqrt(np.divide(np.sum(np.where(success[:,:,None],np.nan_to_num((top-mean_top)**2),0)*weights[:,None,None],axis=0),
                          complete_mass[:,None],out=np.full_like(mean_top,np.nan),where=complete_mass[:,None]>0))
    rows=[row for i in range(k) for row in completed[i+1][0]["models"]]
    frame=pd.DataFrame(rows,columns=["tree_id","tree_topology_hash","model_id","model","fit_status",
        "retried","error","params","lnL","aic","aicc","weight","precision_mode",
        "likelihood_conditioning","fit_method","tree_status","map_id","stratum_weight",
        "rng_task","rng_seed_hex","rng_task_sha256"])
    winners={definition.name:(0.0 if ok.any() else np.nan) for definition in definitions}
    parameter_rows=[]
    if ok.any():
        for i in np.flatnonzero(ok):
            local=completed[i+1][0]["models"]
            winners[local[int(np.argmax([r["weight"] for r in local]))]["model"]]+=weights[i]/weights[ok].sum()
        for j,definition in enumerate(definitions):
            for name in definition.params:
                values=np.array([completed[i+1][0]["models"][j]["params"][name] for i in np.flatnonzero(ok)])
                from .diagnostics import _weighted_quantile
                quantiles=[_weighted_quantile(values,weights[ok],q) for q in (.025,.5,.975)]
                parameter_rows.append(dict(model_id=j+1,model=definition.name,parameter=name,q025=quantiles[0],median=quantiles[1],q975=quantiles[2],
                    interpretation="complete-case cross-tree conditional point-estimate range; not a confidence interval"))
    joint=None
    if bsm_maps_per_tree and ok.any():
        clados,anas,metas,offset=[],[],[],0
        for i in np.flatnonzero(ok):
            record,arrays=completed[i+1]
            for name, target in (("clado",clados),("ana",anas),("maps",metas)):
                table=restore_frame(record[name],arrays)
                table["map"]+=offset
                table["stratum_weight"]*=weights[i]/weights[ok].sum()
                target.append(table)
            offset+=bsm_maps_per_tree
        joint=_make_bsm(pd.concat(clados,ignore_index=True),pd.concat(anas,ignore_index=True),templates[0].areas,offset,pd.concat(metas,ignore_index=True))
        joint.tree_order=tuple(int(i)+1 for i in np.flatnonzero(ok))
        joint.n_trees_missing=int((~ok).sum())
        joint.summary_scope="complete_case" if (~ok).any() else "all_planned_trees"
    result=MultiTreeFitResult(summary,k,int(ok.sum()),int((~ok).sum()),tuple(d.name for d in definitions),
        likelihood_conditioning,templates[0].precision,match_mode,audit_max,weights @ mono,success.sum(axis=0),
        weights @ np.array([m[2] for m in mappings]),mean_top,mean_bottom,bounds_top,bounds_bottom,sd,
        frame,pd.Series(winners,name="complete_case_winner_frequency"),pd.DataFrame(parameter_rows),joint,time.monotonic()-started,
        templates[0].ranges,templates[0].states,tuple(range(1,k+1)),store.task_fingerprint if store else None,init_strategy)
    if init_strategy=="fixed_ref_params":
        result.interpretation="Conditional fixed-reference-parameter summaries; parameters were not reoptimized per tree"
    if (~ok).mean()>max_tree_failure_rate:
        message=f"{int((~ok).sum())}/{k} trees have incomplete model sets; means are complete-case summaries"
        if on_excess_failure=="abort":
            raise TreeFailureError(message,result)
        warnings.warn(message+"; returning deterministic missing-fit bounds",RuntimeWarning)
    return result
