#!/usr/bin/env python3
"""Targeted supplementary experiments addressing independent outline review points.

Covers four reproducible empirical/calibration supplements:
  S1: Pure DEC (j=0) tip-survival simulations across extinction regimes (e=0.02, 0.08, 0.15)
      comparing Contract-C (legacy_unconditioned), DEC* (include_null_range=False), and
      Contract-A (tip_survival) on extinction bias and AICc / Chi2 / Chibar +J false positives.
  S2: Empirical supplements on Hawaii Psychotria (D1 + M3 stratified) and Trimezieae (D2):
      - DEC* / DEC*+J comparisons;
      - Log-scale e_max scan (1..200) on Psychotria D1 & M3 quantifying asymptotic convergence;
      - High-resolution Psychotria M3 parametric bootstrap (B=999) + 5,000-map model-averaged BSM.
  S3: Clean idle-machine repeated thread-scaling benchmarks (repeats=5, threads=1,2,4,8)
      plus Trimezieae (57 tips, 4 states) original R BioGeoBEARS vs RustBEARS timing & Peak RSS.
  S4: Long MCMC tree-chain calibration (K=800, m=20, persistence=0.8, 500 replicates) verifying
      nominal ~95% posterior-population coverage once effective_n_trees >= 50 is satisfied.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import beta

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import rustbears as b
from rustbears._ffi import PreparedSampler
from rustbears.uncertainty.diagnostics import _Design
from rustbears.uncertainty.rng import TaskIdentity, derive_task_seeds
from rustbears.uncertainty.simulate import sampler_buffers
from scripts.run_uncertainty_validation import balanced_tree
from scripts.validate_v12_v4_precision import make_fit, exact_event_counts


def cp95(k: int, n: int) -> list[float]:
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def _s1_worker(args):
    tree_str, areas, tip_dict, e_max = args
    geog = b.io.TipRanges(areas, tip_dict)

    def make_tpl(fam, min_zero=True):
        t = b.model_template(fam).with_parameters(
            d=dict(min=0.0 if min_zero else 1e-12, max=e_max, init=0.015),
            e=dict(min=0.0 if min_zero else 1e-12, max=e_max, init=0.015),
        )
        if "+J" in fam:
            t = t.with_parameters(j=dict(min=0.0 if min_zero else 1e-5, init=0.02))
        return t

    fc0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), likelihood_conditioning="legacy_unconditioned").fit()
    fc1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), likelihood_conditioning="legacy_unconditioned").fit()
    fs0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), include_null_range=False, likelihood_conditioning="legacy_unconditioned").fit()
    fs1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), include_null_range=False, likelihood_conditioning="legacy_unconditioned").fit()
    fa0 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC"), likelihood_conditioning="tip_survival").fit()
    fa1 = b.BioGeoBEARS(tree_str, geog, make_tpl("DEC+J"), likelihood_conditioning="tip_survival").fit()

    return dict(
        e_C=float(fc0.params["e"]),
        e_DECstar=float(fs0.params["e"]),
        e_A=float(fa0.params["e"]),
        d_stat_C=max(0.0, 2.0 * (fc1.lnL - fc0.lnL)),
        d_stat_DECstar=max(0.0, 2.0 * (fs1.lnL - fs0.lnL)),
        d_stat_A=max(0.0, 2.0 * (fa1.lnL - fa0.lnL)),
        aicc_fp_C=bool(fc1.aicc < fc0.aicc),
        aicc_fp_DECstar=bool(fs1.aicc < fs0.aicc),
        aicc_fp_A=bool(fa1.aicc < fa0.aicc),
    )


def run_s1_extinction_and_plusj_fp(workers: int = 32) -> dict:
    tree_str = balanced_tree(30, height=10.0)
    geog_init = {f"T{i+1}": "1000" for i in range(30)}
    base_mod = b.BioGeoBEARS(tree_str, geog_init, "DEC", likelihood_conditioning="tip_survival")
    areas = base_mod.areas
    states = base_mod.states
    n_reps = 200
    regimes = [
        ("low_e002", 0.015, 0.02),
        ("moderate_e008", 0.015, 0.08),
        ("high_e015", 0.015, 0.15),
    ]
    out = {}
    chi2_95 = 3.841458820694124
    chibar_95 = 2.705543454095404  # 50:50 mixture 0.5*chi2_0 + 0.5*chi2_1 at alpha=0.05

    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, d_true, e_true in regimes:
            seed = 20261118 + int(round(e_true * 1000))
            sim = b.simulate_datasets(
                base_mod, params={"d": d_true, "e": e_true},
                n_reps=n_reps, seed=seed, num_threads=min(8, workers),
            )
            tasks = []
            for r in range(n_reps):
                tip_dict = {lbl: states[int(st)] for lbl, st in zip(sim.tip_labels, sim.tip_states[r])}
                tasks.append((tree_str, areas, tip_dict, 5.0))
            rows = list(pool.map(_s1_worker, tasks))
            regime_res = {"d_true": d_true, "e_true": e_true, "j_true": 0.0, "n_reps": n_reps}
            for tag in ("C", "DECstar", "A"):
                e_vals = np.array([row[f"e_{tag}"] for row in rows])
                d_stats = np.array([row[f"d_stat_{tag}"] for row in rows])
                aicc_hits = int(sum(row[f"aicc_fp_{tag}"] for row in rows))
                chi2_hits = int(np.sum(d_stats > chi2_95))
                chibar_hits = int(np.sum(d_stats > chibar_95))
                cens_hits = int(np.sum(e_vals >= 5.0 - 1e-4))
                regime_res[tag] = dict(
                    e_median=float(np.median(e_vals)),
                    e_iqr=[float(np.quantile(e_vals, 0.25)), float(np.quantile(e_vals, 0.75))],
                    e_zero_fraction=float(np.mean(e_vals < 1e-6)),
                    e_upper_censored_fraction=cens_hits / n_reps,
                    aicc_false_positive_rate=aicc_hits / n_reps,
                    aicc_fp_cp95=cp95(aicc_hits, n_reps),
                    chi2_df1_false_positive_rate=chi2_hits / n_reps,
                    chi2_df1_fp_cp95=cp95(chi2_hits, n_reps),
                    chibar_50_50_false_positive_rate=chibar_hits / n_reps,
                    chibar_50_50_fp_cp95=cp95(chibar_hits, n_reps),
                    d_obs_median=float(np.median(d_stats)),
                    d_obs_q90=float(np.quantile(d_stats, 0.90)),
                )
            out[name] = regime_res
    return out


def run_s2_empirical_supplements(d2_dir: Path, workers: int = 32) -> dict:
    bgb_dir = REPO / "r_package/RustBEARS/inst/extdata/examples/Psychotria_M3strat/BGB"
    tree_p = bgb_dir / "Psychotria_5.2.newick"
    geog_p = bgb_dir / "Psychotria_geog.data"
    times_p = bgb_dir / "timeperiods.txt"
    allowed_p = bgb_dir / "areas_allowed.txt"
    manual_p = bgb_dir / "dispersal_multipliers.txt"
    strat_kw = dict(timeperiods=times_p, areas_allowed=allowed_p, dispersal_multipliers=manual_p)

    # 1. e_max log-scale convergence scan on Psychotria D1 and M3 under Contract-A:
    emax_scan = {"D1_single_layer": [], "M3_stratified": []}
    for label, kw in (("D1_single_layer", {}), ("M3_stratified", strat_kw)):
        for e_max in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 150.0, 200.0):
            try:
                fits = {}
                for fam in ("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"):
                    mt = b.model_template(fam).with_parameters(e=dict(max=e_max))
                    fits[fam] = b.BioGeoBEARS(tree_p, geog_p, mt, likelihood_conditioning="tip_survival", **kw).fit(num_threads=4)
                ma = b.model_average(list(fits.values()), criterion="AICc")
                bsm_dec = fits["DEC"].run_bsm(num_maps=500, num_threads=8, seed=20261118)
                no_j_w = float(ma.weights_df[ma.weights_df["model"].isin(["DEC", "DIVALIKE", "BAYAREALIKE"])]["weight"].sum())
                emax_scan[label].append(dict(
                    e_max=e_max,
                    status="ok",
                    lnL_DEC=float(fits["DEC"].lnL),
                    log_survival_DEC=float(fits["DEC"].evaluation.log_tip_survival_prob),
                    e_hat_DEC=float(fits["DEC"].params["e"]),
                    d_hat_DEC=float(fits["DEC"].params["d"]),
                    lnL_DECJ=float(fits["DEC+J"].lnL),
                    log_survival_DECJ=float(fits["DEC+J"].evaluation.log_tip_survival_prob),
                    e_hat_DECJ=float(fits["DEC+J"].params["e"]),
                    j_hat_DECJ=float(fits["DEC+J"].params["j"]),
                    D_obs_DEC_vs_DECJ=float(2.0 * (fits["DEC+J"].lnL - fits["DEC"].lnL)),
                    no_plusJ_total_AICc_weight=no_j_w,
                    bsm_DEC_mean_extinction_events=float(bsm_dec.summary_counts["event_means"]["e"]),
                ))
            except ValueError as exc:
                emax_scan[label].append(dict(
                    e_max=e_max,
                    status="ieee754_float64_survival_underflow",
                    error=str(exc),
                ))

    # 2. DEC* (include_null_range=False) comparisons on Psychotria D1, M3, and Trimezieae D2:
    decstar_comp = {}
    for label, tr, gg, kw in [
        ("Psychotria_D1", tree_p, geog_p, {}),
        ("Psychotria_M3", tree_p, geog_p, strat_kw),
        ("Trimezieae_D2_ref_tree", d2_dir / "reference_clean.tre", d2_dir / "geography.data", {}),
    ]:
        f_dec = b.BioGeoBEARS(tr, gg, "DEC", include_null_range=False, likelihood_conditioning="legacy_unconditioned", **kw).fit(num_threads=4)
        f_decj = b.BioGeoBEARS(tr, gg, "DEC+J", include_null_range=False, likelihood_conditioning="legacy_unconditioned", **kw).fit(num_threads=4)
        decstar_comp[label] = dict(
            DEC_star=dict(lnL=float(f_dec.lnL), aicc=float(f_dec.aicc), d=float(f_dec.params["d"]), e=float(f_dec.params["e"])),
            DEC_star_plusJ=dict(lnL=float(f_decj.lnL), aicc=float(f_decj.aicc), d=float(f_decj.params["d"]), e=float(f_decj.params["e"]), j=float(f_decj.params["j"])),
            D_obs=float(2.0 * (f_decj.lnL - f_dec.lnL)),
        )

    # 3. High-resolution Psychotria M3 (B=999 bootstrap + 5,000-map model-averaged BSM at standard e_max=5.0):
    m3_fits_a = {
        fam: b.BioGeoBEARS(tree_p, geog_p, fam, likelihood_conditioning="tip_survival", **strat_kw).fit(num_threads=4)
        for fam in ("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
    }
    boot_999 = b.bootstrap_lrt(
        m3_fits_a["DEC"], m3_fits_a["DEC+J"], n_boot=999,
        num_threads=workers, seed=20261118,
    )
    ma_m3 = b.model_average(list(m3_fits_a.values()), criterion="AICc")
    bsm_5000 = ma_m3.run_bsm(num_maps=5000, num_threads=min(16, workers), seed=20261118)
    diag_5000 = b.bsm_diagnostics(bsm_5000, target_rmcse=0.05)

    return dict(
        psychotria_emax_scan=emax_scan,
        dec_star_comparisons=decstar_comp,
        psychotria_m3_formal_b999_5000maps=dict(
            d_obs=float(boot_999.d_obs),
            pval_mc_conservative=float(boot_999.pval_mc_conservative),
            pval_mc_ci95=list(boot_999.pval_mc_ci95),
            n_converged=int(boot_999.n_converged),
            n_failed=int(boot_999.n_failed),
            top_routes_5000maps=diag_5000.route_summary_df.sort_values("mean", ascending=False)
                .head(6)[["source_area", "destination_area", "mean", "se_target", "rmcse"]]
                .to_dict(orient="records"),
        ),
    )


def run_s3_clean_scaling_and_d2_timing(d2_dir: Path, rscript: str, oracle_rlib: str) -> dict:
    load_before = os.getloadavg()

    # Build 6-area (64 states) and 8-area (256 states) clean benchmark fixtures:
    tree_20 = balanced_tree(20, height=5.0)
    geog_64 = {f"T{i+1}": "".join("1" if j == (i % 6) else "0" for j in range(6)) for i in range(20)}
    geog_256 = {f"T{i+1}": "".join("1" if j == (i % 8) else "0" for j in range(8)) for i in range(20)}
    fixed_dec = b.model_template("DEC").with_parameters(d=dict(type="fixed", init=0.06), e=dict(type="fixed", init=0.04))
    fit_64 = b.BioGeoBEARS(tree_20, geog_64, fixed_dec, likelihood_conditioning="tip_survival").fit()
    fit_256 = b.BioGeoBEARS(tree_20, geog_256, fixed_dec, likelihood_conditioning="tip_survival").fit()

    bgb_dir = REPO / "r_package/RustBEARS/inst/extdata/examples/Psychotria_M3strat/BGB"
    fit_m3 = b.BioGeoBEARS(
        bgb_dir / "Psychotria_5.2.newick", bgb_dir / "Psychotria_geog.data", "DEC+J",
        likelihood_conditioning="tip_survival",
        timeperiods=bgb_dir / "timeperiods.txt",
        areas_allowed=bgb_dir / "areas_allowed.txt",
        dispersal_multipliers=bgb_dir / "dispersal_multipliers.txt",
    ).fit(num_threads=4)

    scaling = {}
    for case_name, fit_obj, n_maps in [
        ("bsm_64states_2000maps", fit_64, 2000),
        ("bsm_256states_1000maps", fit_256, 1000),
        ("bsm_psychotria_m3_5000maps", fit_m3, 5000),
    ]:
        # Warmup once:
        fit_obj.run_bsm(num_maps=20, num_threads=1, seed=1)
        by_t = {}
        for nt in (1, 2, 4, 8):
            times = []
            for rep in range(5):
                t0 = time.perf_counter()
                fit_obj.run_bsm(num_maps=n_maps, num_threads=nt, seed=1000 + rep)
                times.append(time.perf_counter() - t0)
            by_t[str(nt)] = dict(
                times_sec=times,
                median_sec=float(np.median(times)),
                iqr_sec=[float(np.quantile(times, 0.25)), float(np.quantile(times, 0.75))],
            )
        t1_med = by_t["1"]["median_sec"]
        for nt in ("1", "2", "4", "8"):
            by_t[nt]["speedup_vs_1thread_median"] = float(t1_med / by_t[nt]["median_sec"])
        scaling[case_name] = by_t

    # Trimezieae D2 timing: RustBEARS (8 threads, 20 trees * 6 models + 10,000 maps -> exact & 200-tree full timing)
    # vs unmodified original BioGeoBEARS R oracle on 5 trees * 6 models + 5 * 500 = 2,500 maps!
    d2_trees_all = (d2_dir / "posterior.trees").read_text().splitlines()
    d2_trees_20 = [ln.strip() for ln in d2_trees_all if ln.strip()][:20]
    d2_geog = d2_dir / "geography.data"

    t0 = time.perf_counter()
    mt_20 = b.fit_over_trees(
        d2_trees_20, d2_geog,
        models=["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"],
        likelihood_conditioning="legacy_unconditioned",
        criterion="AICc",
        bsm_maps_per_tree=500,
        num_threads=8,
        seed=20261118,
    )
    rust_20trees_sec = time.perf_counter() - t0

    # Original BioGeoBEARS R oracle on first 5 trees (6 models + 500 BSM maps per tree):
    with tempfile.TemporaryDirectory() as tmp:
        tmp_p = Path(tmp)
        for idx in range(5):
            (tmp_p / f"t_{idx+1}.tre").write_text(d2_trees_20[idx] + "\n")
        r_code = f"""
.libPaths(c({json.dumps(oracle_rlib)}, .libPaths()))
suppressPackageStartupMessages({{
  library(ape)
  library(rexpokit)
  library(cladoRcpp)
  library(BioGeoBEARS)
}})
t0 <- proc.time()[['elapsed']]
for (idx in 1:2) {{
  trfn <- file.path({json.dumps(str(tmp_p))}, paste0('t_', idx, '.tre'))
  geogfn <- {json.dumps(str(d2_geog))}
  best_res <- NULL; best_aicc <- Inf
  for (fam in c('DEC','DEC+J','DIVALIKE','DIVALIKE+J','BAYAREALIKE','BAYAREALIKE+J')) {{
    x <- define_BioGeoBEARS_run(trfn=trfn, geogfn=geogfn, max_range_size=2, include_null_range=TRUE,
                                num_cores_to_use=1, print_optim=FALSE, printlevel=-1)
    p <- x$BioGeoBEARS_model_object@params_table
    base <- sub('\\\\+J$', '', fam)
    if (base == 'DIVALIKE') {{ p['s','type'] <- 'fixed'; p['s',c('init','est')] <- 0; p['ysv','type'] <- '2-j'; p['ys','type'] <- 'ysv*1/2'; p['y','type'] <- 'ysv*1/2'; p['v','type'] <- 'ysv*1/2'; p['mx01v','type'] <- 'fixed'; p['mx01v',c('init','est')] <- 0.5 }}
    if (base == 'BAYAREALIKE') {{ p['s','type'] <- 'fixed'; p['s',c('init','est')] <- 0; p['v','type'] <- 'fixed'; p['v',c('init','est')] <- 0; p['ysv','type'] <- '1-j'; p['ys','type'] <- 'ysv'; p['y','type'] <- 'ysv'; p['mx01y','type'] <- 'fixed'; p['mx01y',c('init','est')] <- 0.9999 }}
    if (grepl('\\\\+J$', fam)) {{ p['j','type'] <- 'free'; p['j',c('init','est')] <- 0.01; p['j','max'] <- if (base=='DEC') 2.99999 else if (base=='DIVALIKE') 1.99999 else 0.99999 }}
    x$BioGeoBEARS_model_object@params_table <- p
    x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
    invisible(capture.output(res <- bears_optim_run(x)))
    k <- if (grepl('\\\\+J$', fam)) 3 else 2
    n <- 57
    aicc <- 2*k - 2*res$total_loglikelihood + 2*k*(k+1)/(n-k-1)
    if (aicc < best_aicc) {{ best_aicc <- aicc; best_res <- res }}
  }}
  invisible(capture.output(inp <- get_inputs_for_stochastic_mapping(best_res)))
  invisible(capture.output(bsm <- runBSM(best_res, stochastic_mapping_inputs_list=inp, nummaps_goal=500, maxtries_per_branch=40000, save_after_every_try=FALSE, savedir=tempdir(), seedval=12345+idx, wait_before_save=0.01, master_nodenum_toPrint=0)))
}}
elapsed <- proc.time()[['elapsed']] - t0
writeLines(sprintf('%.6f', elapsed), file.path({json.dumps(str(tmp_p))}, 'r_elapsed.txt'))
"""
        r_file = tmp_p / "bench_d2_oracle.R"
        r_file.write_text(r_code)
        subprocess.run([rscript, "--vanilla", str(r_file)], check=True)
        oracle_2trees_sec = float((tmp_p / "r_elapsed.txt").read_text().strip())

    rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_rss_mib = rss_kb / (1024.0 * 1024.0) if sys.platform == "darwin" else rss_kb / 1024.0
    load_after = os.getloadavg()

    return dict(
        host_loadavg_before=list(load_before),
        host_loadavg_after=list(load_after),
        clean_thread_scaling_5repeats=scaling,
        trimezieae_d2_timing_comparison=dict(
            rustbears_20trees_6models_10000maps_8threads_sec=rust_20trees_sec,
            rustbears_extrapolated_200trees_100000maps_sec=rust_20trees_sec * 10.0,
            original_biogeobears_2trees_6models_1000maps_1thread_sec=oracle_2trees_sec,
            original_biogeobears_extrapolated_200trees_100000maps_sec=oracle_2trees_sec * 100.0,
            per_tree_speedup_rust8_vs_original1=(oracle_2trees_sec / 2.0) / (rust_20trees_sec / 20.0),
            python_process_peak_rss_mib=peak_rss_mib,
        ),
    )


def run_s4_long_mcmc_chain_calibration() -> dict:
    fits = [make_fit(scale) for scale in (0.6, 0.8, 1.0, 1.2, 1.4)]
    truths_by_tree = [exact_event_counts(f.model, f.params) for f in fits]
    true_pop_mean = {
        kind: float(np.mean([tr[kind] for tr in truths_by_tree]))
        for kind in ("d", "e")
    }
    k_trees = 800
    m_maps = 20
    persistence = 0.8
    n_reps = 500
    tasks = [
        TaskIdentity("v4_long_mcmc", "mcmc800", str(tree), outer_rep=rep, role="draw")
        for rep in range(1, n_reps + 1) for tree in range(k_trees)
    ] + [
        TaskIdentity("v4_long_mcmc", "mcmc800", "selection", outer_rep=rep, role="tree_chain")
        for rep in range(1, n_reps + 1)
    ]
    seeds = derive_task_seeds(20261118, tasks)
    samplers = [PreparedSampler(sampler_buffers(f.model, f.params), f.model.root_sampling_probs, 0) for f in fits]
    metadata = pd.DataFrame(dict(
        map=np.arange(1, k_trees * m_maps + 1),
        tree_id=np.repeat(np.arange(1, k_trees + 1), m_maps),
        model_id=1,
        stratum_weight=1.0 / k_trees,
    ))
    mapped = SimpleNamespace(num_maps=k_trees * m_maps, map_metadata=metadata, tree_order=tuple(range(1, k_trees + 1)))
    design = _Design(mapped, "posterior_population", "mcmc_chain")

    records = {"d": [], "e": []}
    try:
        for rep in range(1, n_reps + 1):
            rng = np.random.default_rng(seeds[TaskIdentity("v4_long_mcmc", "mcmc800", "selection", outer_rep=rep, role="tree_chain")])
            index = rng.integers(0, 5, size=k_trees)
            for i in range(1, k_trees):
                if rng.random() < persistence:
                    index[i] = index[i - 1]
            draws_d, draws_e = [], []
            for i, t in enumerate(index):
                sd = seeds[TaskIdentity("v4_long_mcmc", "mcmc800", str(i), outer_rep=rep, role="draw")]
                raw = samplers[int(t)].draw(1, m_maps, sd, 1, bsm=True)
                ana = raw[1].reshape((-1, 10))
                draws_d.append(np.bincount(ana[ana[:, 7] == 0, 0].astype(int) - 1, minlength=m_maps))
                draws_e.append(np.bincount(ana[ana[:, 7] == 1, 0].astype(int) - 1, minlength=m_maps))
            for kind, arr_list in (("d", draws_d), ("e", draws_e)):
                x = np.concatenate(arr_list)
                s = design.summarize(x, 0.05)
                lo, hi = s["mean_mc_ci_lower"], s["mean_mc_ci_upper"]
                avail = bool(np.isfinite(lo) and np.isfinite(hi))
                records[kind].append(dict(
                    mean=float(s["mean"]),
                    se=float(s["se_target"]),
                    ess=float(s["effective_n_trees"]),
                    avail=avail,
                    hit=bool(lo <= true_pop_mean[kind] <= hi) if avail else False,
                ))
    finally:
        for sp in samplers:
            sp.close()

    res = {"k_trees": k_trees, "maps_per_tree": m_maps, "persistence": persistence, "replicates": n_reps}
    for kind in ("d", "e"):
        recs = records[kind]
        n_avail = int(sum(r["avail"] for r in recs))
        n_hit = int(sum(r["hit"] for r in recs if r["avail"]))
        ess_vals = np.array([r["ess"] for r in recs])
        res[kind] = dict(
            true_population_mean=true_pop_mean[kind],
            n_available=n_avail,
            n_rejected_low_ess=n_reps - n_avail,
            pass_gate_fraction=n_avail / n_reps,
            estimated_ess_median=float(np.median(ess_vals)),
            estimated_ess_iqr=[float(np.quantile(ess_vals, 0.25)), float(np.quantile(ess_vals, 0.75))],
            empirical_sd=float(np.std([r["mean"] for r in recs], ddof=1)),
            rms_estimated_se=float(np.sqrt(np.mean([r["se"] ** 2 for r in recs]))),
            conditional_coverage_when_available=(n_hit / n_avail) if n_avail else None,
            conditional_coverage_cp95=cp95(n_hit, n_avail) if n_avail else None,
        )
    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--d2-dir", type=Path, required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--oracle-rlib", required=True)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output.parent / "stage_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    def run_cached(name, fn):
        p = cache_dir / f"{name}.json"
        if p.exists():
            print(f"[CACHE HIT] {name}", flush=True)
            return json.loads(p.read_text())
        val = fn()
        p.write_text(json.dumps(val, indent=2) + "\n")
        return val

    t0 = time.perf_counter()
    print("[S3] Running clean idle-machine scaling & D2 oracle timing...", flush=True)
    s3 = run_cached("s3", lambda: run_s3_clean_scaling_and_d2_timing(args.d2_dir, args.rscript, args.oracle_rlib))
    print("[S1] Running pure DEC (j=0) extinction & +J false-positive simulations...", flush=True)
    s1 = run_cached("s1", lambda: run_s1_extinction_and_plusj_fp(workers=args.workers))
    print("[S2] Running empirical Psychotria M3 (B=999, e_max scan) & DEC* comparisons...", flush=True)
    s2 = run_cached("s2", lambda: run_s2_empirical_supplements(args.d2_dir, workers=args.workers))
    print("[S4] Running long MCMC tree-chain (K=800, ESS>=50) calibration...", flush=True)
    s4 = run_cached("s4", run_s4_long_mcmc_chain_calibration)

    payload = dict(
        schema="rustbears.reviewer_supplements.v1",
        seed=20261118,
        total_elapsed_sec=time.perf_counter() - t0,
        S1_extinction_and_plusj_false_positives=s1,
        S2_empirical_supplements_emax_and_decstar=s2,
        S3_clean_thread_scaling_and_d2_timing=s3,
        S4_long_mcmc_chain_coverage=s4,
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("ALL_SUPPLEMENTS_DONE")


if __name__ == "__main__":
    main()
