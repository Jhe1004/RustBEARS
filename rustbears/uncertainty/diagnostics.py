"""Mapping precision, realization spread and exploratory paired route contrasts.

Fixed-sample intervals are nominal pointwise normal approximations; posterior
tree-population intervals use Student-t critical values (df = K - 1 for
independent trees, or df = ESS - 1 when MCMC effective_n_trees >= 50, otherwise
explicitly unavailable NaN endpoints). Data-selected route ranking and batch
precision stopping do not supply simultaneous or anytime coverage. Unobserved
and low-count quantities are explicitly not assessed.
"""
from dataclasses import dataclass
import math

import numpy as np
import pandas as pd
from scipy.special import betaincinv
from scipy.stats import norm


MIN_MCMC_EFFECTIVE_TREES = 50.0


def _t_two_sided_crit(alpha, df):
    """Full float64 two-sided Student-t critical multiplier via inverse incomplete beta."""
    v = float(betaincinv(0.5, 0.5 * float(df), 1.0 - float(alpha)))
    return math.sqrt(float(df) * v / (1.0 - v))


def rank_status(lower, upper, tie_margin):
    if not np.isfinite(tie_margin) or tie_margin <= 0:
        raise ValueError("tie_margin must be positive and finite")
    if not np.isfinite([lower, upper]).all() or lower > upper:
        return "unresolved"
    if lower >= -tie_margin and upper <= tie_margin:
        return "near_tie"
    if lower > tie_margin or upper < -tie_margin:
        return "resolved"
    return "unresolved"


def _ips_variance(values):
    """Geyer initial-positive, monotone pair sequence for one scalar series."""
    x = np.asarray(values, float)
    n = len(x)
    centered = x - x.mean()
    gamma = np.array([np.dot(centered[:n-k], centered[k:])/n for k in range(n)])
    if gamma[0] <= 0:
        return 0.0, float("nan")
    pairs = []
    for k in range(0, n-1, 2):
        pair = gamma[k] + gamma[k+1]
        if pair <= 0:
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    tau = max(1.0, (-gamma[0] + 2*sum(pairs))/gamma[0])
    return float(x.var(ddof=1)*tau/n), float(n/tau)


def _weighted_quantile(values, weights, q):
    order = np.argsort(values, kind="stable")
    x, w = np.asarray(values)[order], np.asarray(weights)[order]
    cumulative = np.cumsum(w)/np.sum(w)
    return float(x[min(np.searchsorted(cumulative, q, side="left"), len(x)-1)])


def mapping_metadata(result):
    """Every map, including zero-event maps, must have a provenance row."""
    metadata = getattr(result, "map_metadata", None)
    if metadata is None:
        metadata = pd.DataFrame({"map": np.arange(1, result.num_maps+1), "tree_id": 1,
                                 "model_id": 1, "stratum_weight": 1.0})
    metadata = metadata.copy()
    required = {"map", "tree_id", "model_id", "stratum_weight"}
    if not required.issubset(metadata):
        raise ValueError("mapping metadata lacks map/tree/model/weight columns")
    if len(metadata) != result.num_maps or metadata["map"].duplicated().any():
        raise ValueError("mapping metadata must include each map exactly once")
    if set(metadata["map"]) != set(range(1, result.num_maps+1)):
        raise ValueError("global map ids must cover 1..num_maps")
    metadata = metadata.sort_values("map", kind="stable").reset_index(drop=True)
    if not np.isfinite(metadata.stratum_weight).all() or (metadata.stratum_weight <= 0).any():
        raise ValueError("stratum weights must be positive and finite")
    for _, group in metadata.groupby(["tree_id", "model_id"], sort=False):
        if group.stratum_weight.nunique() != 1:
            raise ValueError("a stratum must have one fixed weight")
    weights = metadata.groupby(["tree_id", "model_id"], sort=False).stratum_weight.first()
    if not np.isclose(weights.sum(), 1, rtol=0, atol=1e-12):
        raise ValueError("stratum weights must sum to one")
    return metadata


class _Design:
    def __init__(self, result, tree_target, tree_dependence):
        self.meta = mapping_metadata(result)
        self.groups = [(key, np.asarray(rows, int)) for key, rows in
                       self.meta.groupby(["tree_id", "model_id"], sort=False).indices.items()]
        self.weights = np.array([self.meta.iloc[i[0]].stratum_weight for _, i in self.groups])
        if any(len(i) < 2 for _, i in self.groups):
            raise ValueError("each active stratum requires at least two maps")
        # Preserve the supplied chain ordering, never topology-sort MCMC draws.
        self.tree_ids = tuple(dict.fromkeys(self.meta.tree_id))
        explicit_order = getattr(result, "tree_order", None)
        if explicit_order is not None:
            if set(explicit_order) != set(self.tree_ids) or len(explicit_order) != len(self.tree_ids):
                raise ValueError("tree_order must enumerate every tree exactly once")
            self.tree_ids = tuple(explicit_order)
        self.tree_weights = np.array([sum(w for ((t, _), _i), w in zip(self.groups, self.weights) if t == tree)
                                      for tree in self.tree_ids])
        self.target, self.dependence = tree_target, tree_dependence
        if tree_target not in ("fixed_sample", "posterior_population"):
            raise ValueError("unknown tree_target")
        if tree_target == "posterior_population":
            if getattr(result, "n_trees_missing", 0):
                raise ValueError("posterior_population precision is unavailable after incomplete-tree selection")
            if tree_dependence not in ("independent", "mcmc_chain"):
                raise ValueError("posterior_population requires explicit tree_dependence")
            if len(self.tree_ids) < (20 if tree_dependence == "mcmc_chain" else 2):
                raise ValueError("insufficient tree samples for the requested dependence estimator")
            if not np.allclose(self.tree_weights, 1/len(self.tree_ids), rtol=0, atol=1e-12):
                raise ValueError("posterior_population requires equal tree weights")
        elif tree_dependence is not None:
            raise ValueError("tree_dependence applies only to posterior_population")

    def summarize(self, x, alpha):
        x = np.asarray(x, float)
        means = np.array([x[i].mean() for _, i in self.groups])
        variance = np.array([x[i].var(ddof=1) if len(i) >= 2 else np.nan for _, i in self.groups])
        size = np.array([len(i) for _, i in self.groups])
        mean = float(self.weights @ means)
        fixed_var = float(np.sum(self.weights**2*variance/size))
        per_tree = np.array([sum(w*mu for ((t, _), _i), w, mu in zip(self.groups, self.weights, means) if t == tree)/tw
                             for tree, tw in zip(self.tree_ids, self.tree_weights)])
        ess = None
        target_var = fixed_var
        if self.target == "posterior_population":
            if self.dependence == "independent":
                target_var = float(per_tree.var(ddof=1)/len(per_tree))
                ess = float(len(per_tree))
            else:
                target_var, ess = _ips_variance(per_tree)
        se = float(np.sqrt(target_var))
        per_map_weight = np.zeros(len(x))
        for (_, i), w in zip(self.groups, self.weights):
            per_map_weight[i] = w/len(i)
        denom = 1-float(self.weights @ self.weights)
        between = max(0.0, float(np.sum(self.weights*(means-mean)**2 -
                      self.weights*(1-self.weights)*variance/size))/denom) if denom > 1e-15 else 0.0
        if len(self.tree_ids)>1:
            tree_noise = np.array([sum((w/tw)**2*v/n for ((t,_),_i),w,v,n in
                                      zip(self.groups,self.weights,variance,size) if t==tree)
                                   for tree,tw in zip(self.tree_ids,self.tree_weights)])
            denom = 1-float(self.tree_weights @ self.tree_weights)
            between = max(0.0, float(np.sum(self.tree_weights*(per_tree-mean)**2 -
                          self.tree_weights*(1-self.tree_weights)*tree_noise))/denom)
        if self.target == "posterior_population":
            if self.dependence == "independent":
                crit = _t_two_sided_crit(alpha, len(per_tree)-1)
                ci_lower, ci_upper = mean-crit*se, mean+crit*se
            elif ess is not None and np.isfinite(ess) and ess >= MIN_MCMC_EFFECTIVE_TREES:
                crit = _t_two_sided_crit(alpha, ess-1.0)
                ci_lower, ci_upper = mean-crit*se, mean+crit*se
            else:
                ci_lower = ci_upper = float("nan")
        else:
            z = norm.ppf(1-alpha/2)
            ci_lower, ci_upper = mean-z*se, mean+z*se
        return dict(mean=mean, se_target=se, mapping_mcse=math.sqrt(fixed_var),
                    rpi_025=_weighted_quantile(x, per_map_weight, .025),
                    rpi_975=_weighted_quantile(x, per_map_weight, .975),
                    mean_mc_ci_lower=ci_lower, mean_mc_ci_upper=ci_upper,
                    var_within=float(self.weights @ variance), var_between_trunc=between,
                    stratum_mean_q025=_weighted_quantile(means, self.weights, .025),
                    stratum_mean_q975=_weighted_quantile(means, self.weights, .975),
                    effective_n_trees=ess)


@dataclass
class BSMDiagnosticsResult:
    sampling_design: str
    stop_reason: str
    tree_target: str | None
    tree_dependence: str | None
    effective_n_trees: float | None
    num_maps_total: int
    precision_target_met: bool
    target_rmcse: float
    tie_margin: float
    max_observed_rmcse: float
    event_summary_df: pd.DataFrame
    route_summary_df: pd.DataFrame
    route_rank_comparisons_df: pd.DataFrame
    ett_df: pd.DataFrame
    interval_scope: str = "nominal pointwise; selected route comparisons are exploratory; not anytime-valid"

    def plot(self, filepath=None):
        import matplotlib.pyplot as plt
        table = self.event_summary_df
        fig, ax = plt.subplots()
        ax.errorbar(table.event_type, table["mean"], yerr=1.96*table.se_target, fmt="o")
        ax.set(ylabel="Conditional mean event count", title="Nominal Monte Carlo precision")
        if filepath:
            fig.savefig(filepath, bbox_inches="tight")
        return fig


def _counts(result, areas):
    n = result.num_maps
    event_names = tuple("deajysvt")
    events = np.zeros((n, len(event_names)))
    routes = np.zeros((n, len(areas), len(areas)))
    for table, route_types in ((result.clado_events_df, {"j"}), (result.ana_events_df, {"d", "a"})):
        for row in table.itertuples(index=False):
            mid = int(row.map)-1
            if not 0 <= mid < n or row.map != mid+1:
                raise ValueError("event map id outside the complete map index")
            if row.event_type not in event_names:
                raise ValueError("unknown event type")
            events[mid, event_names.index(row.event_type)] += 1
            if row.event_type in route_types and row.source_area in areas and row.destination_area in areas:
                routes[mid, areas.index(row.source_area), areas.index(row.destination_area)] += 1
    return events, routes


def bsm_diagnostics(bsm_result, target_rmcse=.05, tie_margin=.25,
                    tree_target="fixed_sample", tree_dependence=None,
                    top_k_routes=5, min_mean_count=1.0, alpha=.05):
    if not np.isfinite([target_rmcse, tie_margin, min_mean_count, alpha]).all() or target_rmcse <= 0 or tie_margin <= 0 or min_mean_count < 0 or not 0 < alpha < 1:
        raise ValueError("invalid precision, equivalence, activity or interval settings")
    if isinstance(top_k_routes, bool) or not isinstance(top_k_routes, int) or top_k_routes < 1:
        raise ValueError("top_k_routes must be a positive integer")
    if bsm_result.num_maps < 2:
        raise ValueError("at least two maps are required")
    design = _Design(bsm_result, tree_target, tree_dependence)
    areas = tuple(getattr(bsm_result, "areas", ()))
    if not areas:
        areas = tuple(bsm_result.summary_counts["all_dispersals_counts_fromto_means"].index)
    events, routes = _counts(bsm_result, areas)

    def row(x, **labels):
        out = design.summarize(x, alpha)
        active = out["mean"] >= min_mean_count and out["mean"] > 0
        out.update(labels, assessed=active, assessment="assessed" if active else "not_assessed",
                   rmcse=out["se_target"]/out["mean"] if out["mean"] > 0 else np.nan)
        if not active:
            out["mean_mc_ci_lower"] = out["mean_mc_ci_upper"] = np.nan
        return out

    event_table = pd.DataFrame([row(events[:, i], event_type=kind) for i, kind in enumerate("deajysvt")])
    route_table = pd.DataFrame([row(routes[:, i, j], source_area=a, destination_area=b)
                                for i, a in enumerate(areas) for j, b in enumerate(areas) if i != j])
    ranks = route_table[route_table.assessed].sort_values(["mean", "source_area", "destination_area"],
                                                         ascending=[False, True, True]).head(top_k_routes) if len(route_table) else route_table
    comparisons = []
    for a, b in zip(list(ranks.itertuples()), list(ranks.itertuples())[1:]):
        diff = routes[:, areas.index(a.source_area), areas.index(a.destination_area)] - routes[:, areas.index(b.source_area), areas.index(b.destination_area)]
        summary = design.summarize(diff, alpha)
        lower, upper = summary["mean_mc_ci_lower"], summary["mean_mc_ci_upper"]
        comparisons.append(dict(route_a=f"{a.source_area}->{a.destination_area}", route_b=f"{b.source_area}->{b.destination_area}",
                                mean_diff=summary["mean"], se_diff=summary["se_target"], ci_lower=lower, ci_upper=upper,
                                tie_margin=tie_margin, zero_excluded=bool(lower > 0 or upper < 0),
                                status=rank_status(lower, upper, tie_margin), effective_n_trees=summary["effective_n_trees"],
                                interval_scope="exploratory selected-pair nominal interval"))
    comparison_table = pd.DataFrame(comparisons, columns=["route_a", "route_b", "mean_diff", "se_diff", "ci_lower", "ci_upper", "tie_margin", "zero_excluded", "status", "effective_n_trees", "interval_scope"])
    all_rows = pd.concat([event_table, route_table], ignore_index=True)
    active_rows = all_rows.loc[all_rows.assessed]
    active = active_rows["rmcse"].to_numpy()
    ci_ok = bool(np.isfinite(active_rows[["mean_mc_ci_lower", "mean_mc_ci_upper"]].to_numpy(float)).all()) if len(active_rows) else False
    met = bool(len(active) and ci_ok and np.isfinite(active).all() and np.all(active <= target_rmcse))
    stop_reason = getattr(bsm_result, "stop_reason", "fixed_count")
    if stop_reason not in ("fixed_count", "target_precision_met", "max_maps_reached"):
        raise ValueError("unknown BSM stop_reason")
    # Event-through-time counts use the same estimator and zero-event map rows.
    tables = [x for x in (bsm_result.clado_events_df, bsm_result.ana_events_df) if len(x)]
    ett = []
    if tables:
        table = pd.concat(tables, ignore_index=True)
        maximum = float(table.time_bp.max())
        edges = np.linspace(0, max(maximum, 1e-12), 21)
        bins = np.clip(np.searchsorted(edges, table.time_bp, side="right")-1, 0, 19)
        counts = np.zeros((bsm_result.num_maps, 20))
        np.add.at(counts, (table["map"].to_numpy(int)-1, bins), 1)
        ett = [row(counts[:, i], time_lower=edges[i], time_upper=edges[i+1]) for i in range(20)]
    kind = "two_stage_trees" if len(design.tree_ids)>1 else "model_stratified" if len(design.groups)>1 else "single_tree_iid"
    return BSMDiagnosticsResult(kind, stop_reason, tree_target if len(design.tree_ids)>1 else None,
        tree_dependence, None, bsm_result.num_maps, met, target_rmcse, tie_margin,
        float(np.max(active)) if len(active) else np.nan, event_table, route_table, comparison_table, pd.DataFrame(ett))


def _make_bsm(clado, ana, areas, num_maps, metadata=None, stop_reason="fixed_count"):
    from ..bsm import BSMResult, count_ana_clado_events
    clado, ana = clado.copy(), ana.copy()
    for table in (clado, ana):
        if len(table) and (not np.isfinite(table["map"].to_numpy(float)).all() or
                           (table["map"].to_numpy(float) % 1 != 0).any()):
            raise ValueError("map ids must be integers")
        table["map"] = table["map"].astype(np.int64)
    result = BSMResult(clado, ana, count_ana_clado_events(clado, ana, areas, num_maps), num_maps)
    result.areas = tuple(areas)
    result.stop_reason = stop_reason
    if metadata is not None:
        result.map_metadata = metadata
        meta = mapping_metadata(result)
        per_map = np.zeros(num_maps)
        for _key, rows in meta.groupby(["tree_id", "model_id"], sort=False).indices.items():
            per_map[rows] = meta.iloc[rows[0]].stratum_weight/len(rows)
        counts, _ = _counts(result, tuple(areas))
        result.summary_counts["event_means"] = pd.Series(per_map @ counts, index=list("deajysvt"))
        result.summary_counts["event_sds"] = pd.Series(np.sqrt(per_map @ (counts-per_map @ counts)**2), index=list("deajysvt"))
        # Existing plotting/count APIs also report the specified mixture weights.
        for prefix, tables in (("ana", [(ana, {"d", "a"})]), ("clado", [(clado, {"j"})]),
                               ("all", [(ana, {"d", "a"}), (clado, {"j"})])):
            x = np.zeros((num_maps, len(areas), len(areas)))
            for table, kinds in tables:
                for row in table.itertuples(index=False):
                    if row.event_type in kinds and row.source_area in areas and row.destination_area in areas:
                        x[int(row.map)-1, areas.index(row.source_area), areas.index(row.destination_area)] += 1
            mean = np.einsum("m,mij->ij", per_map, x)
            sd = np.sqrt(np.einsum("m,mij->ij", per_map, (x-mean)**2))
            key = prefix+"_dispersals_counts_fromto_"
            result.summary_counts[key+"means"] = pd.DataFrame(mean, index=areas, columns=areas)
            result.summary_counts[key+"sds"] = pd.DataFrame(sd, index=areas, columns=areas)
    return result


def run_bsm_adaptive(fit, *, target_rmcse=.05, min_maps=200, batch_size=200,
                     max_maps=10000, num_threads=1, seed=12345, tie_margin=.25,
                     min_mean_count=1.0, checkpoint_file=None, resume=False):
    """Practical batch precision stopping; no anytime-valid coverage claim.

    Every batch persists both complete event tables, including provenance, so
    resuming reconstructs all diagnostics without drawing previous maps again.
    """
    from .._ffi import PreparedSampler
    from ..bsm import _frame, _CLADO, _ANA
    from .simulate import sampler_buffers
    from .checkpoint import CheckpointStore, model_manifest, frame_arrays, restore_frame
    if resume and checkpoint_file is None:
        raise ValueError("resume requires checkpoint_file")
    for value, name in ((min_maps, "min_maps"), (batch_size, "batch_size"), (max_maps, "max_maps")):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError(f"{name} must be an integer >= 2")
    if min_maps > max_maps:
        raise ValueError("min_maps exceeds max_maps")
    if not np.isfinite([target_rmcse, tie_margin, min_mean_count]).all() or target_rmcse <= 0 or tie_margin <= 0 or min_mean_count < 0:
        raise ValueError("invalid adaptive precision settings")
    m = fit.model
    m._require_tier1a()
    m._validate_contract_params(fit.params, force=True)
    if m.likelihood_conditioning=="legacy_unconditioned" and np.any(m.tips[:,0]>0):
        from .contracts import UnsupportedModelContractError
        raise UnsupportedModelContractError("prepared observed BSM with empty tips requires an explicit probability contract")
    manifest = dict(kind="adaptive_bsm", model=model_manifest(fit), target_rmcse=target_rmcse,
                    min_maps=min_maps, max_maps=max_maps, batch_size=batch_size,
                    seed=seed, tie_margin=tie_margin, min_mean_count=min_mean_count)
    store = CheckpointStore(checkpoint_file, manifest, resume) if checkpoint_file else None
    clado, ana, total = [], [], 0
    if store:
        for record in store.records:
            if record["start"] != total+1:
                raise ValueError("adaptive checkpoint batches are not contiguous")
            arrays = store.load_arrays(record["id"])
            clado.append(restore_frame(record["clado"], arrays))
            ana.append(restore_frame(record["ana"], arrays))
            total += record["count"]

    def combined(reason):
        from .averaging import topology_hash
        provenance = dict(tree_id=1, tree_topology_hash=topology_hash(m.tree), model_id=1,
                          stratum_weight=1.0, precision_mode=fit.precision_mode,
                          likelihood_conditioning=fit.likelihood_conditioning)
        metadata = pd.DataFrame([dict(map=i, map_id=i, **provenance) for i in range(1, total+1)])
        result = _make_bsm(pd.concat(clado, ignore_index=True), pd.concat(ana, ignore_index=True), m.areas, total,
                           metadata, stop_reason=reason)
        for table in (result.clado_events_df, result.ana_events_df):
            table["map_id"] = table["map"]
            for key, value in provenance.items():
                table[key] = value
        result.likelihood_conditioning, result.precision_mode = fit.likelihood_conditioning, fit.precision_mode
        result.diagnostics = bsm_diagnostics(result, target_rmcse=target_rmcse, tie_margin=tie_margin,
                                             min_mean_count=min_mean_count)
        return result

    if total >= min_maps:
        existing = combined("max_maps_reached" if total >= max_maps else "fixed_count")
        if existing.diagnostics.precision_target_met:
            existing.stop_reason = existing.diagnostics.stop_reason = "target_precision_met"
            return existing
        if total >= max_maps:
            return existing
    with PreparedSampler(sampler_buffers(m, fit.params), m.root_sampling_probs, 0) as sampler:
        while total < max_maps:
            n = min(batch_size, max_maps-total)
            raw = sampler.draw(total+1, n, seed, num_threads, bsm=True)
            c, a = _frame(raw[0], _CLADO, m), _frame(raw[1], _ANA, m)
            if store:
                cm, ca = frame_arrays(c, "clado")
                am, aa = frame_arrays(a, "ana")
                store.append(dict(id=len(store.records)+1, start=total+1, count=n, clado=cm, ana=am), {**ca, **aa})
            clado.append(c)
            ana.append(a)
            total += n
            if total >= min_maps:
                result = combined("max_maps_reached" if total == max_maps else "fixed_count")
                if result.diagnostics.precision_target_met:
                    result.stop_reason = result.diagnostics.stop_reason = "target_precision_met"
                    return result
    return result
