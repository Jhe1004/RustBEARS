#!/usr/bin/env python3
"""Final v3 supplementary experiments addressing second-round review audit items.

Covers four reproducible empirical/calibration/benchmark supplements:
  V3_1: (a) Complete S1 summary metrics (mean, SD, bias, RMSE, exact zero-LRT fraction);
        (b) Independent neutral birth-death + range-evolution simulations (N=200 trees per
            regime) where losing the last area causes true lineage extinction and extinct
            lineages are pruned from the reconstructed tree (addressing circularity B2),
            plus a DEC*-truth control (N=200);
        (c) Multi-stratum monotone island emergence simulation calibration (N=100 datasets)
            evaluating stratified profile_ci coverage and bootstrap_lrt (B=99) Type-I error.
  V3_2: Strictly fair equal-core (1-thread vs 1-thread) and equal-resource (8-thread vs
        8-process) runtime benchmarks across Original R BioGeoBEARS, RustBEARS R frontend
        (native_whole_tree & legacy_branch_sink), and RustBEARS Python frontend for both
        BSM and high-state (64/128/256-state, single-layer & 4-stratum) MLE optimization.
  V3_3: Empirical census of tip_survival across the 77-study archival benchmark corpus,
        quantifying exact applicability, frequency of box-upper-bound (e_hat -> e_max)
        censoring vs interior convergence, relationship with single-area endemicity, and
        impact on DEC vs DEC+J AICc selection.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.stats import beta

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import rustbears as b
from scripts.run_uncertainty_validation import balanced_tree


def cp95(k: int, n: int) -> list[float]:
    if n <= 0:
        return [float("nan"), float("nan")]
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


# ---------------------------------------------------------------------------
# Independent forward Gillespie birth-death + range evolution tree simulator
# ---------------------------------------------------------------------------
def _dec_cladogenesis_transitions(mask: int, n_areas: int, max_range: int, j_val: float = 0.0) -> list[tuple[int, int, float]]:
    """Standard BioGeoBEARS DEC/DEC+J cladogenetic daughter state weights for non-empty bitmask."""
    bits = [i for i in range(n_areas) if (mask >> i) & 1]
    sz = len(bits)
    assert sz >= 1
    y_wt = max(0.0, (3.0 - j_val) / 3.0) if j_val > 0 else 1.0
    s_wt = y_wt
    v_wt = y_wt
    out: list[tuple[int, int, float]] = []
    if sz == 1:
        # Narrow sympatry: A -> (A, A)
        if y_wt > 0:
            out.append((mask, mask, y_wt))
    else:
        # Subset sympatry: S -> (S, a) and (a, S) for each single area a in S
        if s_wt > 0:
            for a in bits:
                single = 1 << a
                out.append((mask, single, s_wt))
                out.append((single, mask, s_wt))
        # Narrow vicariance: partition S into single area a and S \ {a}
        if v_wt > 0:
            for a in bits:
                single = 1 << a
                comp = mask ^ single
                if sz == 2:
                    # Note: for size 2, iterating both elements generates (a,b) and (b,a); keep exact pair weights
                    out.append((single, comp, v_wt))
                else:
                    out.append((single, comp, v_wt))
                    out.append((comp, single, v_wt))
    if j_val > 0:
        other_bits = [i for i in range(n_areas) if not ((mask >> i) & 1)]
        for b_idx in other_bits:
            single_b = 1 << b_idx
            out.append((mask, single_b, j_val))
            out.append((single_b, mask, j_val))
    # Deduplicate identical ordered pairs by summing weights (matches BioGeoBEARS clado table)
    merged: dict[tuple[int, int], float] = {}
    for l_m, r_m, w in out:
        if bin(l_m).count("1") <= max_range and bin(r_m).count("1") <= max_range:
            merged[(l_m, r_m)] = merged.get((l_m, r_m), 0.0) + w
    total = sum(merged.values())
    return [(k[0], k[1], v / total) for k, v in merged.items()]


def simulate_bd_geog_reconstructed_tree(
    rng: np.random.Generator,
    n_tips: int = 30,
    n_areas: int = 4,
    max_range: int = 4,
    lam: float = 0.30,
    d_rate: float = 0.03,
    e_rate: float = 0.08,
    forbid_single_area_extinction: bool = False,
    max_attempts: int = 5000,
) -> tuple[str, dict[str, str]]:
    """Forward continuous-time birth-death + biogeographic evolution conditioned on n_tips extant species.

    When forbid_single_area_extinction=False (realistic local extinction), losing the final occupied area
    causes immediate species extinction; extinct side lineages are pruned to yield the reconstructed tree.
    When forbid_single_area_extinction=True (DEC* truth control), single-area lineages have zero extinction rate.
    Uses standard Stadler (2011) uniform-time stopping during the sojourn in state k = n_tips.
    """
    nonnull_masks = [
        m for m in range(1, 1 << n_areas) if 1 <= bin(m).count("1") <= max_range
    ]
    clado_cache = {
        m: _dec_cladogenesis_transitions(m, n_areas, max_range, 0.0)
        for m in nonnull_masks
    }

    for _attempt in range(max_attempts):
        root_mask = int(rng.choice(nonnull_masks))
        pairs = clado_cache[root_mask]
        probs = [p[2] for p in pairs]
        idx0 = int(rng.choice(len(pairs), p=probs))
        l_mask, r_mask, _ = pairs[idx0]

        # Node table: id -> dict(parent, children, birth_time, death_time, alive, mask)
        nodes = {
            0: {"children": [1, 2], "birth": 0.0, "end": 0.0, "alive": False, "mask": root_mask},
            1: {"children": [], "birth": 0.0, "end": None, "alive": True, "mask": l_mask},
            2: {"children": [], "birth": 0.0, "end": None, "alive": True, "mask": r_mask},
        }
        alive_ids = [1, 2]
        next_id = 3
        t_now = 0.0
        steps = 0

        while alive_ids and steps < 50000:
            steps += 1
            # Compute total event rate across all alive lineages
            lineage_events = []  # (kind, node_id, new_mask_or_None, rate)
            total_rate = 0.0
            for nid in alive_ids:
                m = nodes[nid]["mask"]
                sz = bin(m).count("1")
                # Speciation
                lineage_events.append(("spec", nid, 0, lam))
                total_rate += lam
                # Range expansion (dispersal to each unoccupied area if sz < max_range)
                if sz < max_range and d_rate > 0.0:
                    for a in range(n_areas):
                        if not ((m >> a) & 1):
                            lineage_events.append(("disp", nid, m | (1 << a), d_rate))
                            total_rate += d_rate
                # Local extinction from each occupied area
                if e_rate > 0.0 and (sz > 1 or not forbid_single_area_extinction):
                    for a in range(n_areas):
                        if (m >> a) & 1:
                            lineage_events.append(("ext", nid, m ^ (1 << a), e_rate))
                            total_rate += e_rate

            if total_rate <= 0.0:
                break
            dt = float(rng.exponential(1.0 / total_rate))

            # If currently at exact target n_tips, sample uniformly inside sojourn before next birth/death
            # if the next event changes lineage count (speciation or lineage extinction), or advance anagenetically
            u = float(rng.random() * total_rate)
            acc = 0.0
            chosen = lineage_events[-1]
            for ev in lineage_events:
                acc += ev[3]
                if u <= acc:
                    chosen = ev
                    break
            kind, nid, target_mask, _ = chosen

            if len(alive_ids) == n_tips and (kind == "spec" or (kind == "ext" and target_mask == 0)):
                # Stop at uniform point in [t_now, t_now + dt] so extant tips are observed at present
                t_stop = t_now + float(rng.random()) * dt
                for aid in alive_ids:
                    nodes[aid]["end"] = t_stop
                # Prune extinct lineages and format reconstructed Newick
                tip_map: dict[int, str] = {}
                for idx_tip, aid in enumerate(alive_ids, start=1):
                    tip_map[aid] = f"T{idx_tip}"

                def build_reconstructed(u_id: int) -> tuple[str | None, float]:
                    nd = nodes[u_id]
                    seg_len = float(nd["end"] - nd["birth"])
                    if u_id in tip_map:
                        return tip_map[u_id], seg_len
                    sub = []
                    for ch in nd["children"]:
                        s_newick, s_len = build_reconstructed(ch)
                        if s_newick is not None:
                            sub.append((s_newick, s_len))
                    if not sub:
                        return None, 0.0
                    if len(sub) == 1:
                        # Extinct sister pruned: collapse degree-2 node and add branch length
                        return sub[0][0], seg_len + sub[0][1]
                    l_str, l_len = sub[0]
                    r_str, r_len = sub[1]
                    return f"({l_str}:{max(l_len, 1e-5):.8f},{r_str}:{max(r_len, 1e-5):.8f})", seg_len

                root_res, _ = build_reconstructed(0)
                # Ensure root bifurcates among extant taxa (crown tree of n_tips)
                if root_res is not None and root_res.startswith("("):
                    geog_dict = {}
                    for aid, lbl in tip_map.items():
                        m = nodes[aid]["mask"]
                        geog_dict[lbl] = tuple(a for a in range(n_areas) if (m >> a) & 1)
                    # Require at least 2 areas represented across tips so dataset is non-degenerate
                    occupied_union = 0
                    for aid in alive_ids:
                        occupied_union |= nodes[aid]["mask"]
                    if bin(occupied_union).count("1") >= 2:
                        return root_res + ";", geog_dict
                break

            t_now += dt
            if kind == "disp":
                nodes[nid]["mask"] = target_mask
            elif kind == "ext":
                if target_mask == 0:
                    # True lineage extinction!
                    nodes[nid]["end"] = t_now
                    nodes[nid]["alive"] = False
                    alive_ids.remove(nid)
                else:
                    nodes[nid]["mask"] = target_mask
            elif kind == "spec":
                m = nodes[nid]["mask"]
                pairs = clado_cache[m]
                idx_c = int(rng.choice(len(pairs), p=[p[2] for p in pairs]))
                lm, rm, _ = pairs[idx_c]
                c1, c2 = next_id, next_id + 1
                next_id += 2
                nodes[nid]["end"] = t_now
                nodes[nid]["alive"] = False
                nodes[nid]["children"] = [c1, c2]
                nodes[c1] = {"children": [], "birth": t_now, "end": None, "alive": True, "mask": lm}
                nodes[c2] = {"children": [], "birth": t_now, "end": None, "alive": True, "mask": rm}
                alive_ids.remove(nid)
                alive_ids.extend([c1, c2])

    raise RuntimeError("Exceeded max_attempts in simulate_bd_geog_reconstructed_tree")


def _fit_six_worker(args):
    tree_str, areas, tip_dict, e_max = args
    geog = b.io.TipRanges(tuple(areas), {k: tuple(v) for k, v in tip_dict.items()})

    def make_tpl(fam):
        t = b.model_template(fam).with_parameters(
            d=dict(min=0.0, max=e_max, init=0.02),
            e=dict(min=0.0, max=e_max, init=0.02),
        )
        if "+J" in fam:
            t = t.with_parameters(j=dict(min=0.0, init=0.02))
        return t

    fc0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), likelihood_conditioning="legacy_unconditioned").fit()
    fc1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), likelihood_conditioning="legacy_unconditioned").fit()
    fs0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), include_null_range=False, likelihood_conditioning="legacy_unconditioned").fit()
    fs1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), include_null_range=False, likelihood_conditioning="legacy_unconditioned").fit()
    fa0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), likelihood_conditioning="tip_survival").fit()
    fa1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), likelihood_conditioning="tip_survival").fit()

    single_frac = float(np.mean([len(s) == 1 for s in tip_dict.values()]))
    return dict(
        single_area_tip_fraction=single_frac,
        d_C=float(fc0.params["d"]),
        e_C=float(fc0.params["e"]),
        d_DECstar=float(fs0.params["d"]),
        e_DECstar=float(fs0.params["e"]),
        d_A=float(fa0.params["d"]),
        e_A=float(fa0.params["e"]),
        d_stat_C=max(0.0, 2.0 * (fc1.lnL - fc0.lnL)),
        d_stat_DECstar=max(0.0, 2.0 * (fs1.lnL - fs0.lnL)),
        d_stat_A=max(0.0, 2.0 * (fa1.lnL - fa0.lnL)),
        aicc_fp_C=bool(fc1.aicc < fc0.aicc),
        aicc_fp_DECstar=bool(fs1.aicc < fs0.aicc),
        aicc_fp_A=bool(fa1.aicc < fa0.aicc),
    )


def _summarize_fit_rows(rows: list[dict], d_true: float, e_true: float, e_max: float = 5.0) -> dict:
    n_reps = len(rows)
    chi2_95 = 3.841458820694124
    chibar_95 = 2.705543454095404
    res = {
        "d_true": d_true,
        "e_true": e_true,
        "j_true": 0.0,
        "n_reps": n_reps,
        "single_area_tip_fraction_median": float(np.median([r["single_area_tip_fraction"] for r in rows])),
    }
    for tag in ("C", "DECstar", "A"):
        e_vals = np.array([r[f"e_{tag}"] for r in rows], dtype=float)
        d_vals = np.array([r[f"d_{tag}"] for r in rows], dtype=float)
        d_stats = np.array([r[f"d_stat_{tag}"] for r in rows], dtype=float)
        uncens = e_vals < (e_max - 1e-4)
        e_uncens = e_vals[uncens] if np.any(uncens) else e_vals
        aicc_hits = int(sum(r[f"aicc_fp_{tag}"] for r in rows))
        chi2_hits = int(np.sum(d_stats > chi2_95))
        chibar_hits = int(np.sum(d_stats > chibar_95))
        cens_hits = int(np.sum(~uncens))
        res[tag] = dict(
            e_median=float(np.median(e_vals)),
            e_iqr=[float(np.quantile(e_vals, 0.25)), float(np.quantile(e_vals, 0.75))],
            e_mean_all=float(np.mean(e_vals)),
            e_sd_all=float(np.std(e_vals, ddof=1)),
            e_rmse_all=float(np.sqrt(np.mean((e_vals - e_true) ** 2))),
            e_mean_uncensored=float(np.mean(e_uncens)),
            e_rmse_uncensored=float(np.sqrt(np.mean((e_uncens - e_true) ** 2))),
            e_bias_median=float(np.median(e_vals) - e_true),
            e_zero_fraction=float(np.mean(e_vals < 1e-5)),
            e_upper_censored_fraction=cens_hits / n_reps,
            d_median=float(np.median(d_vals)),
            d_rmse_all=float(np.sqrt(np.mean((d_vals - d_true) ** 2))),
            d_stat_exact_zero_fraction=float(np.mean(d_stats < 1e-6)),
            d_obs_median=float(np.median(d_stats)),
            d_obs_q90=float(np.quantile(d_stats, 0.90)),
            aicc_false_positive_rate=aicc_hits / n_reps,
            aicc_fp_cp95=cp95(aicc_hits, n_reps),
            chi2_df1_false_positive_rate=chi2_hits / n_reps,
            chi2_df1_fp_cp95=cp95(chi2_hits, n_reps),
            chibar_50_50_false_positive_rate=chibar_hits / n_reps,
            chibar_50_50_fp_cp95=cp95(chibar_hits, n_reps),
        )
    return res


def _strat_cal_worker(args):
    rep_idx, tree_str, areas, tip_dict, times, allowed_rows, d_true, e_true, n_boot = args
    geog = b.io.TipRanges(tuple(areas), {k: tuple(v) for k, v in tip_dict.items()})

    def make_tpl(fam):
        t = b.model_template(fam).with_parameters(
            d=dict(min=0.0, max=5.0, init=0.05),
            e=dict(min=0.0, max=5.0, init=0.05),
        )
        if "+J" in fam:
            t = t.with_parameters(j=dict(min=0.0, init=0.02))
        return t

    f0 = b.BioGeoBEARS(
        tree_str, geog, make_tpl("DEC"),
        timeperiods=times, areas_allowed=allowed_rows,
        likelihood_conditioning="tip_survival",
    ).fit(num_threads=1)
    f1 = b.BioGeoBEARS(
        tree_str, geog, make_tpl("DEC+J"),
        timeperiods=times, areas_allowed=allowed_rows,
        likelihood_conditioning="tip_survival",
    ).fit(num_threads=1)
    prof = f0.profile_ci(params=["d", "e"], n_grid=12, num_threads=1)
    df_p = prof.to_dataframe().set_index("param")

    def covers(row, truth):
        lo = float(row["lower"]) if pd.notna(row["lower"]) else 0.0
        is_cens = str(row["status_upper"]) == "censored_at_box_upper"
        if is_cens or pd.isna(row["upper"]):
            hi = float("inf")
        else:
            hi = float(row["upper"])
        return bool(lo <= truth <= hi), is_cens

    hit_d, cens_d = covers(df_p.loc["d"], d_true)
    hit_e, cens_e = covers(df_p.loc["e"], e_true)
    boot = b.bootstrap_lrt(f0, f1, n_boot=n_boot, num_threads=1, seed=20261200 + rep_idx)
    d_obs = max(0.0, 2.0 * (f1.lnL - f0.lnL))
    return dict(
        d_hat=float(f0.params["d"]),
        e_hat=float(f0.params["e"]),
        hit_d=hit_d,
        cens_d=cens_d,
        hit_e=hit_e,
        cens_e=cens_e,
        d_obs=d_obs,
        boot_pval=float(boot.pval_mc_conservative),
        aicc_fp=bool(f1.aicc < f0.aicc),
    )


import pandas as pd


def run_v3_1_simulations(workers: int = 32) -> dict:
    # 1a. Complete S1 metrics on fixed balanced tree (exact same seed as S1)
    tree_str = balanced_tree(30, height=10.0)
    geog_init = {f"T{i+1}": "1000" for i in range(30)}
    base_mod = b.BioGeoBEARS(tree_str, geog_init, "DEC", likelihood_conditioning="tip_survival")
    areas = list(base_mod.areas)
    states = list(base_mod.states)
    n_reps = 200
    s1_extended = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, d_true, e_true in [
            ("low_e002", 0.015, 0.02),
            ("moderate_e008", 0.015, 0.08),
            ("high_e015", 0.015, 0.15),
        ]:
            seed = 20261118 + int(round(e_true * 1000))
            sim = b.simulate_datasets(
                base_mod, params={"d": d_true, "e": e_true},
                n_reps=n_reps, seed=seed, num_threads=min(8, workers),
            )
            tasks = [
                (tree_str, areas, {lbl: states[int(st)] for lbl, st in zip(sim.tip_labels, sim.tip_states[r])}, 5.0)
                for r in range(n_reps)
            ]
            rows = list(pool.map(_fit_six_worker, tasks))
            s1_extended[name] = _summarize_fit_rows(rows, d_true, e_true, 5.0)

        # 1b. Independent neutral birth-death + biogeographic extinction pruned trees (N=200 per regime)
        bd_results = {}
        bd_regimes = [
            ("bd_pruned_moderate_e008", 0.30, 0.03, 0.08, False, 20262001),
            ("bd_pruned_high_e015", 0.35, 0.04, 0.15, False, 20262002),
            ("decstar_truth_control_e010", 0.30, 0.03, 0.10, True, 20262003),
        ]
        bd_areas = ["A", "B", "C", "D"]
        for rname, lam, d_true, e_true, forbid_single, base_seed in bd_regimes:
            tasks = []
            for r in range(200):
                rng = np.random.default_rng(base_seed + r)
                tr_s, g_d = simulate_bd_geog_reconstructed_tree(
                    rng, n_tips=30, n_areas=4, max_range=4,
                    lam=lam, d_rate=d_true, e_rate=e_true,
                    forbid_single_area_extinction=forbid_single,
                )
                tasks.append((tr_s, bd_areas, g_d, 5.0))
            rows = list(pool.map(_fit_six_worker, tasks))
            summary = _summarize_fit_rows(rows, d_true, e_true, 5.0)
            summary["speciation_lambda"] = lam
            summary["forbid_single_area_extinction_truth"] = forbid_single
            bd_results[rname] = summary

        # 1c. Stratified monotone island emergence calibration (N=100 datasets, 4 strata, B=99)
        strat_tree = balanced_tree(30, height=8.0)
        times = [1.0, 3.0, 6.0, 10.0]
        allowed_rows = [
            np.outer(v, v).tolist()
            for v in (
                [1, 1, 1, 1],  # 0-1 Ma: all 4 islands
                [1, 1, 1, 0],  # 1-3 Ma: A, B, C
                [1, 1, 0, 0],  # 3-6 Ma: A, B
                [1, 0, 0, 0],  # 6-10 Ma (root at 8 Ma): oldest island A only
            )
        ]
        d_st, e_st = 0.08, 0.06
        strat_mod = b.BioGeoBEARS(
            strat_tree, {f"T{i+1}": "1000" for i in range(30)}, "DEC",
            timeperiods=times, areas_allowed=allowed_rows,
            likelihood_conditioning="tip_survival",
        )
        sim_st = b.simulate_datasets(
            strat_mod, params={"d": d_st, "e": e_st},
            n_reps=100, seed=20263001, num_threads=min(8, workers),
        )
        st_areas = list(strat_mod.areas)
        st_states = list(strat_mod.states)
        st_tasks = [
            (
                r, strat_tree, st_areas,
                {lbl: st_states[int(st)] for lbl, st in zip(sim_st.tip_labels, sim_st.tip_states[r])},
                times, allowed_rows, d_st, e_st, 99,
            )
            for r in range(100)
        ]
        st_rows = list(pool.map(_strat_cal_worker, st_tasks))

    n_st = len(st_rows)
    hits_d = int(sum(r["hit_d"] for r in st_rows))
    hits_e = int(sum(r["hit_e"] for r in st_rows))
    boot_fp = int(sum(r["boot_pval"] <= 0.05 for r in st_rows))
    chi2_fp = int(sum(r["d_obs"] > 3.841458820694124 for r in st_rows))
    chibar_fp = int(sum(r["d_obs"] > 2.705543454095404 for r in st_rows))
    aicc_fp = int(sum(r["aicc_fp"] for r in st_rows))
    strat_cal = dict(
        n_reps=n_st,
        n_tips=30,
        n_areas=4,
        timeperiods=times,
        d_true=d_st,
        e_true=e_st,
        j_true=0.0,
        bootstrap_B=99,
        d_median=float(np.median([r["d_hat"] for r in st_rows])),
        e_median=float(np.median([r["e_hat"] for r in st_rows])),
        profile_ci_coverage_d=hits_d / n_st,
        profile_ci_coverage_d_cp95=cp95(hits_d, n_st),
        profile_ci_upper_censored_d_fraction=float(np.mean([r["cens_d"] for r in st_rows])),
        profile_ci_coverage_e=hits_e / n_st,
        profile_ci_coverage_e_cp95=cp95(hits_e, n_st),
        profile_ci_upper_censored_e_fraction=float(np.mean([r["cens_e"] for r in st_rows])),
        bootstrap_lrt_type1_error=boot_fp / n_st,
        bootstrap_lrt_type1_cp95=cp95(boot_fp, n_st),
        chi2_df1_type1_error=chi2_fp / n_st,
        chibar_50_50_type1_error=chibar_fp / n_st,
        aicc_false_positive_rate=aicc_fp / n_st,
    )

    return dict(
        s1_extended_metrics=s1_extended,
        independent_birth_death_pruned_simulations=bd_results,
        stratified_monotone_emergence_calibration=strat_cal,
    )


# ---------------------------------------------------------------------------
# V3_2: Fair equal-core (1 vs 1) & equal-resource (8 vs 8) runtime benchmarks
# ---------------------------------------------------------------------------
def _write_phylip_geog(path: Path, tip_dict: dict[str, str], areas: list[str]) -> None:
    lines = [f"{len(tip_dict)}\t{len(areas)}\t({' '.join(areas)})"]
    for k, v in tip_dict.items():
        lines.append(f"{k}\t{v}")
    path.write_text("\n".join(lines) + "\n")


def run_v3_2_fair_benchmarks(rscript: str, oracle_rlib: str, rust_rlib: str, deps_rlib: str) -> dict:
    load_before = list(os.getloadavg())
    rng = np.random.default_rng(20261118)

    # Prepare synthetic 64-state (6 areas), 128-state (7 areas), and 256-state (8 areas) datasets
    # both single-layer and 4-stratum!
    bgb_dir = REPO / "r_package/RustBEARS/inst/extdata/examples/Psychotria_M3strat/BGB"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_p = Path(tmp)
        datasets = {}
        # Psychotria single-layer & M3 stratified
        datasets["psychotria_16st_single"] = {
            "tree": str(bgb_dir / "Psychotria_5.2.newick"),
            "geog": str(bgb_dir / "Psychotria_geog.data"),
            "max_range": 4,
            "strat": False,
        }
        datasets["psychotria_16st_m3strat"] = {
            "tree": str(bgb_dir / "Psychotria_5.2.newick"),
            "geog": str(bgb_dir / "Psychotria_geog.data"),
            "max_range": 4,
            "strat": True,
            "timeperiods": str(bgb_dir / "timeperiods.txt"),
            "areas_allowed": str(bgb_dir / "areas_allowed.txt"),
            "dispersal_multipliers": str(bgb_dir / "dispersal_multipliers.txt"),
        }

        # Synthetic 64, 128, 256 state datasets (20 tips, height 6.0, realistic tip ranges)
        tr20 = balanced_tree(20, height=6.0)
        tr20_p = tmp_p / "tree20.newick"
        tr20_p.write_text(tr20 + "\n")
        tp4_p = tmp_p / "tp4.txt"
        tp4_p.write_text("1.5\n3.0\n4.5\n10.0\n")

        for n_a, label_st in [(6, "64st"), (7, "128st"), (8, "256st")]:
            areas = [chr(ord("A") + i) for i in range(n_a)]
            tip_d = {}
            for i in range(20):
                # 1 or 2 areas occupied so interior d, e are identifiable
                k_occ = 1 if (i % 3 != 0) else 2
                occ = rng.choice(n_a, size=k_occ, replace=False)
                s = "".join("1" if a in occ else "0" for a in range(n_a))
                tip_d[f"T{i+1}"] = s
            gp = tmp_p / f"geog_{label_st}.data"
            _write_phylip_geog(gp, tip_d, areas)
            # 4-stratum all-allowed + mild multiplier matrix so original R strat_fn runs cleanly
            aa_p = tmp_p / f"aa_{label_st}.txt"
            hdr = "\t".join(areas)
            aa_blocks = []
            dm_blocks = []
            for s_idx in range(4):
                aa_blocks.append(hdr + "\n" + "\n".join(["\t".join(["1"] * n_a)] * n_a) + "\n")
                mat_rows = []
                for r_i in range(n_a):
                    row_vals = ["1.0" if r_i == c_i else f"{0.8 - 0.1 * s_idx:.2f}" for c_i in range(n_a)]
                    mat_rows.append("\t".join(row_vals))
                dm_blocks.append(hdr + "\n" + "\n".join(mat_rows) + "\n")
            aa_p.write_text("\n".join(aa_blocks) + "\nEND\n")
            dm_p = tmp_p / f"dm_{label_st}.txt"
            dm_p.write_text("\n".join(dm_blocks) + "\nEND\n")

            datasets[f"syn_{label_st}_single"] = {
                "tree": str(tr20_p), "geog": str(gp), "max_range": n_a, "strat": False,
            }
            datasets[f"syn_{label_st}_4strat"] = {
                "tree": str(tr20_p), "geog": str(gp), "max_range": n_a, "strat": True,
                "timeperiods": str(tp4_p), "areas_allowed": str(aa_p), "dispersal_multipliers": str(dm_p),
            }

        # 1. Fair BSM benchmark (200 maps, 5 repeats) across Original R (1-core & 8-proc mclapply),
        #    RustBEARS R (native_whole_tree 1-thread & 8-thread; legacy_branch_sink 1-thread),
        #    and RustBEARS Python (1-thread & 8-thread).
        bsm_cases = ["psychotria_16st_single", "psychotria_16st_m3strat", "syn_64st_single", "syn_128st_single"]
        bsm_results = {}
        for cname in bsm_cases:
            cfg = datasets[cname]
            n_maps = 200
            reps = 5
            kw = {}
            if cfg["strat"]:
                kw = dict(
                    timeperiods=cfg["timeperiods"],
                    areas_allowed=cfg["areas_allowed"],
                    dispersal_multipliers=cfg["dispersal_multipliers"],
                )
            py_fit = b.BioGeoBEARS(
                cfg["tree"], cfg["geog"], "DEC+J",
                max_range_size=cfg["max_range"],
                likelihood_conditioning="legacy_unconditioned", **kw,
            ).fit(num_threads=4)
            py_fit.run_bsm(num_maps=10, num_threads=1, seed=1)  # warmup
            py_t1, py_t8 = [], []
            for rep in range(reps):
                t0 = time.perf_counter()
                py_fit.run_bsm(num_maps=n_maps, num_threads=1, seed=500 + rep)
                py_t1.append(time.perf_counter() - t0)
                t0 = time.perf_counter()
                py_fit.run_bsm(num_maps=n_maps, num_threads=8, seed=500 + rep)
                py_t8.append(time.perf_counter() - t0)

            # R RustBEARS (native_whole_tree 1 & 8 threads + legacy_branch_sink 1 thread)
            out_rust_r = tmp_p / f"r_rust_bsm_{cname}.json"
            strat_lines = ""
            if cfg["strat"]:
                strat_lines = f"""
x$timesfn <- {json.dumps(cfg['timeperiods'])}
x$areas_allowed_fn <- {json.dumps(cfg['areas_allowed'])}
x$dispersal_multipliers_fn <- {json.dumps(cfg['dispersal_multipliers'])}
x <- readfiles_BioGeoBEARS_run(x)
x <- section_the_tree(inputs=x, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
"""
            r_rust_script = f"""
.libPaths(c({json.dumps(rust_rlib)}, {json.dumps(deps_rlib)}, .libPaths()))
suppressPackageStartupMessages({{ library(ape); library(RustBEARS) }})
x <- define_BioGeoBEARS_run(trfn={json.dumps(cfg['tree'])}, geogfn={json.dumps(cfg['geog'])},
                            max_range_size={cfg['max_range']}, include_null_range=TRUE,
                            num_cores_to_use=1, print_optim=FALSE, printlevel=-1)
p <- x$BioGeoBEARS_model_object@params_table
p['j','type'] <- 'free'; p['j',c('init','est')] <- 0.02; p['j','max'] <- 2.99999
x$BioGeoBEARS_model_object@params_table <- p
x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
{strat_lines}
res_p <- x$BioGeoBEARS_model_object@params_table
res_p['d','est'] <- {float(py_fit.params['d'])}; res_p['d','init'] <- {float(py_fit.params['d'])}
res_p['e','est'] <- {float(py_fit.params['e'])}; res_p['e','init'] <- {float(py_fit.params['e'])}
res_p['j','est'] <- {float(py_fit.params['j'])}; res_p['j','init'] <- {float(py_fit.params['j'])}
x$BioGeoBEARS_model_object@params_table <- res_p
x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=FALSE)
x$return_condlikes_table <- TRUE; x$calc_TTL_loglike_from_condlikes_table <- TRUE; x$calc_ancprobs <- TRUE
invisible(capture.output(res <- bears_optim_run(x, skip_optim=TRUE, skip_optim_option='return_all')))
invisible(capture.output(runBSM(res, NULL, nummaps_goal=10, bsm_engine='native_whole_tree', num_threads=1L, save_after_every_try=FALSE, savedir=tempdir(), seedval=1L)))
t_nw1 <- numeric({reps}); t_nw8 <- numeric({reps}); t_leg1 <- numeric({reps})
for (r in seq_len({reps})) {{
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(runBSM(res, NULL, nummaps_goal={n_maps}, bsm_engine='native_whole_tree', num_threads=1L, save_after_every_try=FALSE, savedir=tempdir(), seedval=600L+r)))
  t_nw1[r] <- proc.time()[['elapsed']] - s0
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(runBSM(res, NULL, nummaps_goal={n_maps}, bsm_engine='native_whole_tree', num_threads=8L, save_after_every_try=FALSE, savedir=tempdir(), seedval=600L+r)))
  t_nw8[r] <- proc.time()[['elapsed']] - s0
}}
invisible(capture.output(inp <- get_inputs_for_stochastic_mapping(res)))
for (r in seq_len({reps})) {{
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(runBSM(res, stochastic_mapping_inputs_list=inp, nummaps_goal={n_maps}, bsm_engine='legacy_r_bridge', maxtries_per_branch=40000, save_after_every_try=FALSE, savedir=tempdir(), seedval=600L+r)))
  t_leg1[r] <- proc.time()[['elapsed']] - s0
}}
writeLines(sprintf('{{"nw1":[%s],"nw8":[%s],"leg1":[%s]}}',
  paste(sprintf('%.6f', t_nw1), collapse=','),
  paste(sprintf('%.6f', t_nw8), collapse=','),
  paste(sprintf('%.6f', t_leg1), collapse=',')), {json.dumps(str(out_rust_r))})
"""
            rs_file = tmp_p / f"run_rust_bsm_{cname}.R"
            rs_file.write_text(r_rust_script)
            r_env_rust = {**os.environ, "R_LIBS_USER": f"{rust_rlib}:{deps_rlib}"}
            subprocess.run([rscript, "--vanilla", str(rs_file)], check=True, env=r_env_rust)
            rr_d = json.loads(out_rust_r.read_text())

            # Original BioGeoBEARS R (1-core & 8-process parallel::mclapply splitting 200 maps = 8x25)
            out_ora_r = tmp_p / f"r_ora_bsm_{cname}.json"
            ora_reps = 3 if "128st" in cname else reps
            r_ora_script = f"""
.libPaths(c({json.dumps(oracle_rlib)}, {json.dumps(deps_rlib)}, .libPaths()))
suppressPackageStartupMessages({{ library(ape); library(rexpokit); library(cladoRcpp); library(BioGeoBEARS); library(parallel) }})
x <- define_BioGeoBEARS_run(trfn={json.dumps(cfg['tree'])}, geogfn={json.dumps(cfg['geog'])},
                            max_range_size={cfg['max_range']}, include_null_range=TRUE,
                            num_cores_to_use=1, print_optim=FALSE, printlevel=-1)
p <- x$BioGeoBEARS_model_object@params_table
p['j','type'] <- 'free'; p['j',c('init','est')] <- 0.02; p['j','max'] <- 2.99999
x$BioGeoBEARS_model_object@params_table <- p
x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
{strat_lines}
# Set fixed MLE params from Python fit so we time pure BSM mapping + input prep identically:
res_p <- x$BioGeoBEARS_model_object@params_table
res_p['d','est'] <- {float(py_fit.params['d'])}; res_p['d','init'] <- {float(py_fit.params['d'])}
res_p['e','est'] <- {float(py_fit.params['e'])}; res_p['e','init'] <- {float(py_fit.params['e'])}
res_p['j','est'] <- {float(py_fit.params['j'])}; res_p['j','init'] <- {float(py_fit.params['j'])}
x$BioGeoBEARS_model_object@params_table <- res_p
x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=FALSE)
x$return_condlikes_table <- TRUE; x$calc_TTL_loglike_from_condlikes_table <- TRUE; x$calc_ancprobs <- TRUE
invisible(capture.output(res <- bears_optim_run(x, skip_optim=TRUE, skip_optim_option='return_all')))
invisible(capture.output(inp <- get_inputs_for_stochastic_mapping(res)))
t_ora1 <- numeric({ora_reps}); t_ora8 <- numeric({ora_reps})
for (r in seq_len({ora_reps})) {{
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(runBSM(res, stochastic_mapping_inputs_list=inp, nummaps_goal={n_maps}, maxtries_per_branch=40000, save_after_every_try=FALSE, savedir=tempdir(), seedval=700L+r, wait_before_save=0.0, master_nodenum_toPrint=0)))
  t_ora1[r] <- proc.time()[['elapsed']] - s0
  s0 <- proc.time()[['elapsed']]
  chunks <- mclapply(1:8, function(k) {{
    td <- file.path(tempdir(), paste0('w_', r, '_', k)); dir.create(td, showWarnings=FALSE, recursive=TRUE)
    invisible(capture.output(runBSM(res, stochastic_mapping_inputs_list=inp, nummaps_goal=25, maxtries_per_branch=40000, save_after_every_try=FALSE, savedir=td, seedval=7000L+r*10L+k, wait_before_save=0.0, master_nodenum_toPrint=0)))
    TRUE
  }}, mc.cores=8L)
  t_ora8[r] <- proc.time()[['elapsed']] - s0
}}
writeLines(sprintf('{{"ora1":[%s],"ora8":[%s]}}',
  paste(sprintf('%.6f', t_ora1), collapse=','),
  paste(sprintf('%.6f', t_ora8), collapse=',')), {json.dumps(str(out_ora_r))})
"""
            ro_file = tmp_p / f"run_ora_bsm_{cname}.R"
            ro_file.write_text(r_ora_script)
            r_env_ora = {**os.environ, "R_LIBS_USER": f"{oracle_rlib}:{deps_rlib}"}
            subprocess.run([rscript, "--vanilla", str(ro_file)], check=True, env=r_env_ora)
            ora_d = json.loads(out_ora_r.read_text())

            ora1_med = float(np.median(ora_d["ora1"]))
            ora8_med = float(np.median(ora_d["ora8"]))
            py1_med = float(np.median(py_t1))
            py8_med = float(np.median(py_t8))
            rnw1_med = float(np.median(rr_d["nw1"]))
            rnw8_med = float(np.median(rr_d["nw8"]))
            rleg1_med = float(np.median(rr_d["leg1"]))

            bsm_results[cname] = dict(
                num_maps=n_maps,
                original_r_1core_median_sec=ora1_med,
                original_r_8proc_mclapply_median_sec=ora8_med,
                rustbears_r_legacy_1core_median_sec=rleg1_med,
                rustbears_r_native_1thread_median_sec=rnw1_med,
                rustbears_r_native_8thread_median_sec=rnw8_med,
                rustbears_py_native_1thread_median_sec=py1_med,
                rustbears_py_native_8thread_median_sec=py8_med,
                speedup_1v1_r_native_vs_original=ora1_med / rnw1_med,
                speedup_1v1_py_native_vs_original=ora1_med / py1_med,
                speedup_1v1_r_legacy_vs_original=ora1_med / rleg1_med,
                speedup_8v8_r_native_vs_original_mclapply=ora8_med / rnw8_med,
                speedup_8v8_py_native_vs_original_mclapply=ora8_med / py8_med,
                speedup_8v1_py_native_vs_original_1core=ora1_med / py8_med,
                speedup_8v1_r_native_vs_original_1core=ora1_med / rnw8_med,
            )

        # 2. Fair MLE optimization benchmark (DEC+J, 1-core vs 1-core and 8-core vs 8-core)
        mle_cases = [
            "syn_64st_single", "syn_128st_single", "syn_256st_single",
            "syn_64st_4strat", "syn_128st_4strat", "syn_256st_4strat",
        ]
        mle_results = {}
        for cname in mle_cases:
            cfg = datasets[cname]
            reps = 3
            kw = {}
            if cfg["strat"]:
                kw = dict(
                    timeperiods=cfg["timeperiods"],
                    areas_allowed=cfg["areas_allowed"],
                    dispersal_multipliers=cfg["dispersal_multipliers"],
                )
            py1_t, py8_t = [], []
            for _ in range(reps):
                t0 = time.perf_counter()
                f_py1 = b.BioGeoBEARS(
                    cfg["tree"], cfg["geog"], "DEC+J",
                    max_range_size=cfg["max_range"],
                    likelihood_conditioning="legacy_unconditioned", **kw,
                ).fit(num_threads=1)
                py1_t.append(time.perf_counter() - t0)
                t0 = time.perf_counter()
                f_py8 = b.BioGeoBEARS(
                    cfg["tree"], cfg["geog"], "DEC+J",
                    max_range_size=cfg["max_range"],
                    likelihood_conditioning="legacy_unconditioned", **kw,
                ).fit(num_threads=8)
                py8_t.append(time.perf_counter() - t0)

            strat_lines = ""
            if cfg["strat"]:
                strat_lines = f"""
x$timesfn <- {json.dumps(cfg['timeperiods'])}
x$areas_allowed_fn <- {json.dumps(cfg['areas_allowed'])}
x$dispersal_multipliers_fn <- {json.dumps(cfg['dispersal_multipliers'])}
x <- readfiles_BioGeoBEARS_run(x)
x <- section_the_tree(inputs=x, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
"""
            # RustBEARS R MLE (1 thread & 8 threads)
            out_rr_mle = tmp_p / f"rr_mle_{cname}.json"
            rr_mle_code = f"""
.libPaths(c({json.dumps(rust_rlib)}, {json.dumps(deps_rlib)}, .libPaths()))
suppressPackageStartupMessages({{ library(ape); library(RustBEARS) }})
run_one <- function(nc) {{
  x <- define_BioGeoBEARS_run(trfn={json.dumps(cfg['tree'])}, geogfn={json.dumps(cfg['geog'])},
                              max_range_size={cfg['max_range']}, include_null_range=TRUE,
                              num_cores_to_use=nc, print_optim=FALSE, printlevel=-1)
  p <- x$BioGeoBEARS_model_object@params_table
  p['j','type'] <- 'free'; p['j',c('init','est')] <- 0.02; p['j','max'] <- 2.99999
  x$BioGeoBEARS_model_object@params_table <- p
  x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
  {strat_lines}
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(res <- bears_optim_run(x)))
  c(proc.time()[['elapsed']] - s0, res$total_loglikelihood)
}}
t1 <- numeric({reps}); t8 <- numeric({reps}); lnl <- 0
for (r in seq_len({reps})) {{
  v1 <- run_one(1L); t1[r] <- v1[1]; lnl <- v1[2]
  v8 <- run_one(8L); t8[r] <- v8[1]
}}
writeLines(sprintf('{{"t1":[%s],"t8":[%s],"lnl":%.10f}}',
  paste(sprintf('%.6f', t1), collapse=','),
  paste(sprintf('%.6f', t8), collapse=','), lnl), {json.dumps(str(out_rr_mle))})
"""
            f_rr = tmp_p / f"run_rr_mle_{cname}.R"
            f_rr.write_text(rr_mle_code)
            r_env_rust = {**os.environ, "R_LIBS_USER": f"{rust_rlib}:{deps_rlib}"}
            subprocess.run([rscript, "--vanilla", str(f_rr)], check=True, env=r_env_rust)
            rr_m = json.loads(out_rr_mle.read_text())

            # Original BioGeoBEARS R MLE (1 core & 8 cores via optimx cluster)
            ora_mle_reps = 1 if ("256st" in cname) else 3
            out_ora_mle = tmp_p / f"ora_mle_{cname}.json"
            ora_mle_code = f"""
.libPaths(c({json.dumps(oracle_rlib)}, {json.dumps(deps_rlib)}, .libPaths()))
suppressPackageStartupMessages({{ library(ape); library(rexpokit); library(cladoRcpp); library(BioGeoBEARS); library(parallel) }})
run_one <- function(nc) {{
  x <- define_BioGeoBEARS_run(trfn={json.dumps(cfg['tree'])}, geogfn={json.dumps(cfg['geog'])},
                              max_range_size={cfg['max_range']}, include_null_range=TRUE,
                              num_cores_to_use=nc, print_optim=FALSE, printlevel=-1)
  p <- x$BioGeoBEARS_model_object@params_table
  p['j','type'] <- 'free'; p['j',c('init','est')] <- 0.02; p['j','max'] <- 2.99999
  x$BioGeoBEARS_model_object@params_table <- p
  x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
  {strat_lines}
  s0 <- proc.time()[['elapsed']]
  invisible(capture.output(res <- bears_optim_run(x)))
  c(proc.time()[['elapsed']] - s0, res$total_loglikelihood)
}}
t1 <- numeric({ora_mle_reps}); t8 <- numeric({ora_mle_reps}); lnl <- 0
for (r in seq_len({ora_mle_reps})) {{
  v1 <- run_one(1L); t1[r] <- v1[1]; lnl <- v1[2]
  v8 <- run_one(8L); t8[r] <- v8[1]
}}
writeLines(sprintf('{{"t1":[%s],"t8":[%s],"lnl":%.10f}}',
  paste(sprintf('%.6f', t1), collapse=','),
  paste(sprintf('%.6f', t8), collapse=','), lnl), {json.dumps(str(out_ora_mle))})
"""
            f_oo = tmp_p / f"run_ora_mle_{cname}.R"
            f_oo.write_text(ora_mle_code)
            r_env_ora = {**os.environ, "R_LIBS_USER": f"{oracle_rlib}:{deps_rlib}"}
            subprocess.run([rscript, "--vanilla", str(f_oo)], check=True, env=r_env_ora)
            oo_m = json.loads(out_ora_mle.read_text())

            ora1 = float(np.median(oo_m["t1"]))
            ora8 = float(np.median(oo_m["t8"]))
            rr1 = float(np.median(rr_m["t1"]))
            rr8 = float(np.median(rr_m["t8"]))
            py1 = float(np.median(py1_t))
            py8 = float(np.median(py8_t))

            mle_results[cname] = dict(
                original_r_1core_median_sec=ora1,
                original_r_8core_median_sec=ora8,
                rustbears_r_1thread_median_sec=rr1,
                rustbears_r_8thread_median_sec=rr8,
                rustbears_py_1thread_median_sec=py1,
                rustbears_py_8thread_median_sec=py8,
                lnL_original_r=float(oo_m["lnl"]),
                lnL_rustbears_r=float(rr_m["lnl"]),
                lnL_rustbears_py=float(f_py1.lnL),
                speedup_1v1_r_vs_original=ora1 / rr1,
                speedup_1v1_py_vs_original=ora1 / py1,
                speedup_8v8_r_vs_original=ora8 / rr8,
                speedup_8v8_py_vs_original=ora8 / py8,
            )

    return dict(
        host_loadavg_before=load_before,
        host_loadavg_after=list(os.getloadavg()),
        bsm_equal_core_and_resource=bsm_results,
        mle_equal_core_and_resource=mle_results,
    )


# ---------------------------------------------------------------------------
# V3_3: Empirical census across the 77-study archival benchmark corpus
# ---------------------------------------------------------------------------
def _corpus_worker(args):
    screen, m8_root_str, item = args
    m8_root = Path(m8_root_str)
    files = item["files"]
    tree_p = m8_root / files["tree"]
    geog_p = m8_root / files["geog"]
    max_r = item.get("max_range_size")
    kw = {}
    if max_r is not None:
        kw["max_range_size"] = int(max_r)
    for k_opt, k_arg in [
        ("timeperiods", "timeperiods"),
        ("areas_allowed", "areas_allowed"),
        ("dispersal_multipliers", "dispersal_multipliers"),
        ("distances", "distances"),
        ("areas_adjacency", "areas_adjacency"),
    ]:
        if k_opt in files and files[k_opt]:
            kw[k_arg] = m8_root / files[k_opt]

    try:
        mod_test = b.BioGeoBEARS(tree_p, geog_p, "DEC", likelihood_conditioning="tip_survival", **kw)
        n_states = len(mod_test.states)
        n_tips = len(mod_test.tree.tip_labels)
        n_areas = len(mod_test.areas)
        if n_states > 512:
            return dict(screen=screen, status="skipped_large_state_space", n_states=n_states, n_tips=n_tips, n_areas=n_areas)

        # Compute single-area endemic tip fraction
        single_frac = float(np.mean([len(v) == 1 for v in mod_test.geog.ranges.values()]))

        def make_tpl(fam):
            t = b.model_template(fam).with_parameters(
                d=dict(min=0.0, max=5.0, init=0.015),
                e=dict(min=0.0, max=5.0, init=0.015),
            )
            if "+J" in fam:
                t = t.with_parameters(j=dict(min=0.0, init=0.02))
            return t

        fc0 = b.BioGeoBEARS(tree_p, geog_p, make_tpl("DEC"), likelihood_conditioning="legacy_unconditioned", **kw).fit(num_threads=1)
        fc1 = b.BioGeoBEARS(tree_p, geog_p, make_tpl("DEC+J"), likelihood_conditioning="legacy_unconditioned", **kw).fit(num_threads=1)
        fa0 = b.BioGeoBEARS(tree_p, geog_p, make_tpl("DEC"), likelihood_conditioning="tip_survival", **kw).fit(num_threads=1)
        fa1 = b.BioGeoBEARS(tree_p, geog_p, make_tpl("DEC+J"), likelihood_conditioning="tip_survival", **kw).fit(num_threads=1)

        return dict(
            screen=screen,
            status="fitted",
            n_tips=n_tips,
            n_areas=n_areas,
            n_states=n_states,
            single_area_tip_fraction=single_frac,
            e_C_dec=float(fc0.params["e"]),
            e_C_decj=float(fc1.params["e"]),
            e_A_dec=float(fa0.params["e"]),
            e_A_decj=float(fa1.params["e"]),
            log_surv_A_dec=float(fa0.log_tip_survival_prob),
            log_surv_A_decj=float(fa1.log_tip_survival_prob),
            aicc_plusj_preferred_C=bool(fc1.aicc < fc0.aicc),
            aicc_plusj_preferred_A=bool(fa1.aicc < fa0.aicc),
            d_lrt_C=max(0.0, 2.0 * (fc1.lnL - fc0.lnL)),
            d_lrt_A=max(0.0, 2.0 * (fa1.lnL - fa0.lnL)),
        )
    except Exception as exc:
        return dict(
            screen=screen,
            status="rejected_or_unsupported",
            error_type=type(exc).__name__,
            error_msg=str(exc)[:240],
        )


def run_v3_3_empirical_corpus_census(m8_dir: Path, workers: int = 32) -> dict:
    inputs = json.loads((m8_dir / "source/empirical-inputs.json").read_text())
    by_screen: dict[str, list[dict]] = {}
    for x in inputs:
        by_screen.setdefault(x["screen"], []).append(x)

    rep_tasks = []
    invalid_screens = []
    for sc in sorted(by_screen):
        lst = [x for x in by_screen[sc] if x.get("valid")]
        if not lst:
            invalid_screens.append(sc)
            continue
        ft = [x for x in lst if x.get("track") == "FullTree_Direct"]
        pool_lst = ft if ft else lst
        b0 = [x for x in pool_lst if x.get("axis") == "baseline"]
        rep = (b0 if b0 else pool_lst)[0]
        rep_tasks.append((sc, str(m8_dir), rep))

    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_corpus_worker, rep_tasks))

    fitted = [r for r in rows if r["status"] == "fitted"]
    rejected = [r for r in rows if r["status"] == "rejected_or_unsupported"]
    skipped_large = [r for r in rows if r["status"] == "skipped_large_state_space"]

    n_fit = len(fitted)
    e_a_dec = np.array([r["e_A_dec"] for r in fitted], dtype=float)
    e_a_decj = np.array([r["e_A_decj"] for r in fitted], dtype=float)
    e_c_dec = np.array([r["e_C_dec"] for r in fitted], dtype=float)
    sf = np.array([r["single_area_tip_fraction"] for r in fitted], dtype=float)

    cens_dec_mask = e_a_dec >= 4.99
    cens_decj_mask = e_a_decj >= 4.99
    zero_a_dec = e_a_dec < 1e-5
    zero_c_dec = e_c_dec < 1e-5

    return dict(
        total_archival_studies=len(by_screen),
        invalid_upstream_inputs=len(invalid_screens),
        evaluated_valid_studies=len(rep_tasks),
        fitted_under_tip_survival=n_fit,
        skipped_over_512_states=len(skipped_large),
        explicitly_rejected_studies=len(rejected),
        rejection_reasons=[{"screen": r["screen"], "error_type": r["error_type"], "msg": r["error_msg"]} for r in rejected],
        dec_extinction_summary=dict(
            legacy_zero_fraction=float(np.mean(zero_c_dec)),
            legacy_zero_count=int(np.sum(zero_c_dec)),
            legacy_median_e=float(np.median(e_c_dec)),
            tip_survival_zero_fraction=float(np.mean(zero_a_dec)),
            tip_survival_zero_count=int(np.sum(zero_a_dec)),
            tip_survival_interior_fraction=float(np.mean((~zero_a_dec) & (~cens_dec_mask))),
            tip_survival_interior_count=int(np.sum((~zero_a_dec) & (~cens_dec_mask))),
            tip_survival_upper_censored_fraction=float(np.mean(cens_dec_mask)),
            tip_survival_upper_censored_count=int(np.sum(cens_dec_mask)),
            tip_survival_median_e_all=float(np.median(e_a_dec)),
            tip_survival_median_e_interior=float(np.median(e_a_dec[(~zero_a_dec) & (~cens_dec_mask)])) if np.any((~zero_a_dec) & (~cens_dec_mask)) else None,
            single_area_fraction_median_when_censored=float(np.median(sf[cens_dec_mask])) if np.any(cens_dec_mask) else None,
            single_area_fraction_median_when_uncensored=float(np.median(sf[~cens_dec_mask])) if np.any(~cens_dec_mask) else None,
        ),
        decj_extinction_summary=dict(
            tip_survival_upper_censored_fraction=float(np.mean(cens_decj_mask)),
            tip_survival_upper_censored_count=int(np.sum(cens_decj_mask)),
        ),
        plusj_model_selection_shift=dict(
            legacy_aicc_prefers_decj_count=int(sum(r["aicc_plusj_preferred_C"] for r in fitted)),
            legacy_aicc_prefers_decj_fraction=float(np.mean([r["aicc_plusj_preferred_C"] for r in fitted])),
            tip_survival_aicc_prefers_decj_count=int(sum(r["aicc_plusj_preferred_A"] for r in fitted)),
            tip_survival_aicc_prefers_decj_fraction=float(np.mean([r["aicc_plusj_preferred_A"] for r in fitted])),
            median_lrt_D_legacy=float(np.median([r["d_lrt_C"] for r in fitted])),
            median_lrt_D_tip_survival=float(np.median([r["d_lrt_A"] for r in fitted])),
        ),
        per_study_records=fitted,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--m8-dir", type=Path, required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--oracle-rlib", required=True)
    parser.add_argument("--rust-rlib", required=True)
    parser.add_argument("--deps-rlib", required=True)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output.parent / "v3_stage_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    def run_cached(name, fn):
        p = cache_dir / f"{name}.json"
        if p.exists():
            print(f"[CACHE HIT] {name}", flush=True)
            return json.loads(p.read_text())
        t_s = time.perf_counter()
        val = fn()
        print(f"[DONE] {name} in {time.perf_counter() - t_s:.1f}s", flush=True)
        p.write_text(json.dumps(val, indent=2) + "\n")
        return val

    t0 = time.perf_counter()
    print("[V3_2] Running fair 1v1 and 8v8 R/Python/Original benchmarks...", flush=True)
    v3_2 = run_cached("v3_2_benchmarks", lambda: run_v3_2_fair_benchmarks(
        args.rscript, args.oracle_rlib, args.rust_rlib, args.deps_rlib,
    ))
    print("[V3_3] Running 77-study archival empirical census under tip_survival...", flush=True)
    v3_3 = run_cached("v3_3_corpus", lambda: run_v3_3_empirical_corpus_census(
        args.m8_dir, workers=args.workers,
    ))
    print("[V3_1] Running S1 extended metrics, N=200 neutral BD pruned trees & N=100 stratified calibration...", flush=True)
    v3_1 = run_cached("v3_1_simulations", lambda: run_v3_1_simulations(workers=args.workers))

    payload = dict(
        schema="rustbears.reviewer_v3_supplements.v1",
        seed=20261118,
        total_elapsed_sec=time.perf_counter() - t0,
        V3_1_simulations_and_stratified_calibration=v3_1,
        V3_2_fair_equal_core_and_resource_benchmarks=v3_2,
        V3_3_archival_77study_tip_survival_census=v3_3,
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("ALL_V3_SUPPLEMENTS_DONE", flush=True)


if __name__ == "__main__":
    main()
