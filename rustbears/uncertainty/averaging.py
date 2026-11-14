"""Conditional ML summaries under frequency-based Akaike model weights."""
from dataclasses import dataclass
import hashlib
import json

import numpy as np
import pandas as pd

from .guards import validate_compatible_fits


def node_clades(tree):
    result = {}
    stack = [(tree.root, False)]
    while stack:
        node, visited = stack.pop()
        children = tree.children[node]
        if not children:
            result[node] = frozenset((tree.labels[node],))
        elif visited:
            result[node] = frozenset().union(*(result[c] for c in children))
        else:
            stack.append((node, True))
            stack.extend((c, False) for c in children)
    return tuple(result[i] for i in range(len(tree.labels)))


def topology_hash(tree):
    clades = sorted(tuple(sorted(c)) for c in node_clades(tree))
    return hashlib.sha256(json.dumps(clades, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def align_node_array(array, source_tree, target_tree):
    if array is None:
        raise ValueError("model averaging requires fitted ancestral probabilities")
    source = {c: i for i, c in enumerate(node_clades(source_tree))}
    target = node_clades(target_tree)
    if set(source) != set(target):
        raise ValueError("node arrays cannot be aligned across different topologies")
    return np.asarray(array)[[source[c] for c in target]].copy()


def allocate_maps(weights, total):
    if isinstance(total, bool) or not isinstance(total, (int, np.integer)):
        raise ValueError("num_maps must be an integer")
    active = np.flatnonzero(np.asarray(weights)>0)
    if total < 2*len(active):
        raise ValueError("num_maps must provide at least two maps per active model")
    counts = np.zeros(len(weights), dtype=int)
    extra = (total-2*len(active))*np.asarray(weights)[active]
    counts[active] = 2+np.floor(extra).astype(int)
    remainder = int(total-counts.sum())
    ranked = np.argsort(-(extra-np.floor(extra)), kind="stable")
    counts[active[ranked[:remainder]]] += 1
    return counts


def sample_strata(fits, weights, counts, *, seed, tree_id=1, purpose="model_average_bsm", num_threads=1, task_seeds=None):
    from .._ffi import PreparedSampler
    from ..bsm import _frame, _CLADO, _ANA
    from .simulate import sampler_buffers
    from .rng import TaskIdentity, derive_task_seeds, task_fingerprint
    from .diagnostics import _make_bsm
    tasks = [TaskIdentity(purpose=purpose, model_id=str(i+1), tree_id=str(tree_id), role="draw")
             for i, n in enumerate(counts) if n]
    seeds = derive_task_seeds(seed, tasks) if task_seeds is None else task_seeds
    clado, ana, metadata, offset = [], [], [], 0
    for task, index in zip(tasks, np.flatnonzero(counts)):
        fit, n, weight = fits[index], int(counts[index]), float(weights[index])
        m = fit.model
        m._require_tier1a()
        m._validate_contract_params(fit.params, force=True)
        if m.likelihood_conditioning=="legacy_unconditioned" and np.any(m.tips[:,0]>0):
            from .contracts import UnsupportedModelContractError
            raise UnsupportedModelContractError("prepared observed BSM with empty tips requires an explicit A/B probability contract; legacy empty-clade scores must not be completed silently")
        with PreparedSampler(sampler_buffers(m, fit.params), m.root_sampling_probs, 0) as sampler:
            raw = sampler.draw(1, n, seeds[task], num_threads, bsm=True)
        c, a = _frame(raw[0], _CLADO, m), _frame(raw[1], _ANA, m)
        provenance = dict(tree_id=tree_id, tree_topology_hash=topology_hash(m.tree), model_id=int(index)+1,
                          stratum_weight=weight, precision_mode=fit.precision_mode,
                          likelihood_conditioning=fit.likelihood_conditioning,
                          rng_seed_hex=f"{seeds[task]:016x}", rng_task_sha256=task_fingerprint(seed,task))
        for table in (c, a):
            table["map_id"] = table["map"].astype(int)
            table["map"] += offset
            for key, value in provenance.items():
                table[key] = value
        clado.append(c)
        ana.append(a)
        metadata.extend(dict(map=offset+j, map_id=j, **provenance) for j in range(1, n+1))
        offset += n
    result = _make_bsm(pd.concat(clado, ignore_index=True), pd.concat(ana, ignore_index=True),
                       fits[0].model.areas, offset, pd.DataFrame(metadata))
    return result


@dataclass
class ModelAverageResult:
    fits: tuple
    criterion: str
    likelihood_conditioning: str
    precision_mode: str
    weights_df: pd.DataFrame
    ancestral_probs_branch_top: np.ndarray
    ancestral_probs_branch_bottom: np.ndarray
    node_sensitivity_df: pd.DataFrame
    states_list: tuple
    ranges_list: tuple
    tree: object
    interpretation: str = "Akaike-weighted conditional ML summary; not a Bayesian model posterior or parameter integral"

    def to_dataframe(self):
        frame = pd.DataFrame(self.ancestral_probs_branch_top, columns=self.ranges_list)
        frame.insert(0, "node", np.arange(len(frame))+1)
        return frame

    def run_bsm(self, num_maps=500, num_threads=1, seed=12345):
        validate_compatible_fits(self.fits)
        weights = self.weights_df.weight.to_numpy()
        return sample_strata(self.fits, weights, allocate_maps(weights, num_maps), seed=seed, num_threads=num_threads)

    def plot(self, filepath=None):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        image = ax.imshow(self.ancestral_probs_branch_top, aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(len(self.ranges_list)), self.ranges_list, rotation=90)
        ax.set(xlabel="Range", ylabel="Reference-tree node", title="Akaike-weighted ancestral probabilities")
        fig.colorbar(image, ax=ax)
        if filepath:
            fig.savefig(filepath, bbox_inches="tight")
        return fig


def model_average(fits, criterion="AICc"):
    fits = validate_compatible_fits(fits)
    if criterion not in ("AIC", "AICc"):
        raise ValueError("criterion must be AIC or AICc")
    scores = np.array([f.aic if criterion == "AIC" else f.aicc for f in fits])
    if not np.isfinite(scores).all():
        raise ValueError("all candidate scores must be finite; undefined AICc requires an explicit choice of AIC")
    weights = np.exp(-.5*(scores-scores.min()))
    weights /= weights.sum()
    tree = fits[0].model.tree
    top = np.array([align_node_array(f.evaluation.ancestral_probs_branch_top, f.model.tree, tree) for f in fits])
    bottom = np.array([align_node_array(f.evaluation.ancestral_probs_branch_bottom, f.model.tree, tree) for f in fits])
    if not np.isfinite(top).all() or not np.allclose(top.sum(axis=2), 1, atol=1e-10):
        raise ValueError("invalid ancestral probability arrays")
    top_mean, bottom_mean = np.einsum("g,gns->ns", weights, top), np.einsum("g,gns->ns", weights, bottom)
    sensitivity = np.sqrt(np.einsum("g,gns->ns", weights, (top-top_mean)**2))
    table = pd.DataFrame([dict(model_id=i+1, model=f.model.model.name, lnL=f.lnL,
                              n_params=f.n_params, aic=f.aic, aicc=f.aicc, weight=float(w),
                              likelihood_conditioning=f.likelihood_conditioning, precision_mode=f.precision_mode)
                          for i, (f, w) in enumerate(zip(fits, weights))])
    nodes = pd.DataFrame([dict(node=i+1, range=state, probability=top_mean[i,j], model_sd=sensitivity[i,j],
                               min_model_probability=top[:,i,j].min(), max_model_probability=top[:,i,j].max())
                          for i in range(len(tree.labels)) for j, state in enumerate(fits[0].model.ranges)])
    return ModelAverageResult(fits, criterion, fits[0].likelihood_conditioning, fits[0].precision_mode,
                              table, top_mean, bottom_mean, nodes, fits[0].model.states, fits[0].model.ranges, tree)
