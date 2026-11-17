#!/usr/bin/env python3
"""Cross-language (R vs Python) stratified parity check and empirical Hawaii
Psychotria M3 (time-stratified four-island emergence) end-to-end validation."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import rustbears as b


def run_cross_language_stratified_parity(rscript_bin: str, r_libs: str) -> dict:
    tree_txt = "((T1:0.31234567,T2:0.31234567):0.48765433,(T3:0.52345678,T4:0.52345678):0.27654322);"
    geog = {"T1": "100", "T2": "010", "T3": "001", "T4": "110"}
    times = [0.25, 0.65, 2.0]
    allowed = [
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
        [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
        [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
    ]
    adj = [[[1, 1, 0], [1, 1, 1], [0, 1, 1]]] * 3
    dist = [
        [[1.0, 2.0, 4.5], [2.0, 1.0, 2.2], [4.5, 2.2, 1.0]],
        [[1.0, 2.6, 5.0], [2.6, 1.0, 2.8], [5.0, 2.8, 1.0]],
        [[1.0, 3.2, 6.0], [3.2, 1.0, 3.5], [6.0, 3.5, 1.0]],
    ]
    env = [
        [[1.0, 1.3, 2.1], [1.3, 1.0, 1.5], [2.1, 1.5, 1.0]],
        [[1.0, 1.6, 2.4], [1.6, 1.0, 1.7], [2.4, 1.7, 1.0]],
        [[1.0, 1.9, 2.8], [1.9, 1.0, 2.0], [2.8, 2.0, 1.0]],
    ]
    manual = [
        [[1.0, 0.85, 0.35], [0.90, 1.0, 0.75], [0.25, 0.65, 1.0]],
        [[1.0, 0.55, 0.05], [0.50, 1.0, 0.05], [0.05, 0.05, 1.0]],
        [[1.0, 0.05, 0.05], [0.05, 1.0, 0.05], [0.05, 0.05, 1.0]],
    ]
    area_sz = [[1.6, 1.1, 0.7], [1.4, 0.8, 0.5], [1.2, 0.5, 0.5]]

    tpl = b.model_template("DEC+J+x+w+n+u").with_parameters(
        d=dict(type="fixed", init=0.18),
        e=dict(type="fixed", init=0.11),
        j=dict(type="fixed", init=0.14),
        x=dict(type="fixed", init=-0.4),
        w=dict(type="fixed", init=1.2),
        n=dict(type="fixed", init=-0.3),
        u=dict(type="fixed", init=-0.5),
    )

    py_evals = {}
    for contract in ("tip_survival", "unconditional_complete", "legacy_unconditioned"):
        m = b.BioGeoBEARS(
            tree_txt, geog, tpl,
            timeperiods=times, areas_allowed=allowed, areas_adjacency=adj,
            distances=dist, envdistances=env, dispersal_multipliers=manual,
            area_of_areas=area_sz, likelihood_conditioning=contract,
            allow_null_tips=(contract == "unconditional_complete"),
        )
        py_evals[contract] = m.fit()

    py_sim = b.simulate_datasets(py_evals["tip_survival"], n_reps=200, seed=20261118, num_threads=2)
    py_bsm = b.runBSM(py_evals["tip_survival"], nummaps_goal=100, seedval=909, num_threads=2)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        r_out = tmp_dir / "r_results.json"
        r_sim_csv = tmp_dir / "r_sim.csv"
        r_clado_csv = tmp_dir / "r_clado.csv"
        r_ana_csv = tmp_dir / "r_ana.csv"
        r_code = f"""
suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64')
folder <- {json.dumps(str(tmp_dir))}
writeLines({json.dumps(tree_txt)}, file.path(folder, 'tree.newick'))
writeLines(c('4 3 (A B C)', 'T1 100', 'T2 010', 'T3 001', 'T4 110'), file.path(folder, 'geog.data'))
writeLines(c('0.25', '0.65', '2.0'), file.path(folder, 'times.txt'))
write_blocks <- function(mats, path) {{
    lines <- character()
    for (m in mats) {{
        lines <- c(lines, 'A\\tB\\tC')
        for (r in seq_len(nrow(m))) lines <- c(lines, paste(m[r, ], collapse='\\t'))
        lines <- c(lines, '')
    }}
    writeLines(c(lines, 'END'), path)
}}
write_blocks(list(
    matrix(c(1,1,1,1,1,1,1,1,1),3,3,byrow=TRUE),
    matrix(c(1,1,0,1,1,0,0,0,0),3,3,byrow=TRUE),
    matrix(c(1,0,0,0,0,0,0,0,0),3,3,byrow=TRUE)), file.path(folder, 'allowed.txt'))
write_blocks(replicate(3, matrix(c(1,1,0,1,1,1,0,1,1),3,3,byrow=TRUE), simplify=FALSE), file.path(folder, 'adj.txt'))
write_blocks(list(
    matrix(c(0,2.0,4.5,2.0,0,2.2,4.5,2.2,0),3,3,byrow=TRUE),
    matrix(c(0,2.6,5.0,2.6,0,2.8,5.0,2.8,0),3,3,byrow=TRUE),
    matrix(c(0,3.2,6.0,3.2,0,3.5,6.0,3.5,0),3,3,byrow=TRUE)), file.path(folder, 'dist.txt'))
write_blocks(list(
    matrix(c(0,1.3,2.1,1.3,0,1.5,2.1,1.5,0),3,3,byrow=TRUE),
    matrix(c(0,1.6,2.4,1.6,0,1.7,2.4,1.7,0),3,3,byrow=TRUE),
    matrix(c(0,1.9,2.8,1.9,0,2.0,2.8,2.0,0),3,3,byrow=TRUE)), file.path(folder, 'env.txt'))
write_blocks(list(
    matrix(c(1.0,0.85,0.35,0.90,1.0,0.75,0.25,0.65,1.0),3,3,byrow=TRUE),
    matrix(c(1.0,0.55,0.05,0.50,1.0,0.05,0.05,0.05,1.0),3,3,byrow=TRUE),
    matrix(c(1.0,0.05,0.05,0.05,1.0,0.05,0.05,0.05,1.0),3,3,byrow=TRUE)), file.path(folder, 'manual.txt'))
writeLines(c('A\\tB\\tC','1.6\\t1.1\\t0.7','','A\\tB\\tC','1.4\\t0.8\\t0.5','','A\\tB\\tC','1.2\\t0.5\\t0.5','','END'),
           file.path(folder, 'areas.txt'))

eval_one <- function(contract) {{
    x <- define_BioGeoBEARS_run(
        trfn=file.path(folder,'tree.newick'), geogfn=file.path(folder,'geog.data'),
        timesfn=file.path(folder,'times.txt'), areas_allowed_fn=file.path(folder,'allowed.txt'),
        areas_adjacency_fn=file.path(folder,'adj.txt'), distsfn=file.path(folder,'dist.txt'),
        dispersal_multipliers_fn=file.path(folder,'manual.txt'),
        area_of_areas_fn=file.path(folder,'areas.txt'), max_range_size=3,
        print_optim=FALSE, printlevel=-1, tmpwd=folder, calc_ancprobs=TRUE)
    x$envdistsfn <- file.path(folder,'env.txt')
    x$wd <- folder
    x$likelihood_conditioning <- contract
    x$allow_null_tips <- contract == 'unconditional_complete'
    p <- x$BioGeoBEARS_model_object@params_table
    vals <- c(d=0.18, e=0.11, j=0.14, x=-0.4, w=1.2, n=-0.3, u=-0.5)
    for (k in names(vals)) {{ p[k,'type'] <- 'fixed'; p[k,c('init','est')] <- vals[[k]] }}
    p['ysv','type'] <- '3-j'; p['ys','type'] <- 'ysv*2/3'
    p['y','type'] <- 'ysv*1/3'; p['s','type'] <- 'ysv*1/3'; p['v','type'] <- 'ysv*1/3'
    x$BioGeoBEARS_model_object@params_table <- p
    x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init=TRUE)
    invisible(capture.output(x <- readfiles_BioGeoBEARS_run(x)))
    invisible(capture.output(x <- section_the_tree(inputs=x, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)))
    res <- NULL
    invisible(capture.output(res <- bears_optim_run(x, skip_optim=TRUE, skip_optim_option='return_all')))
    res
}}
ra <- eval_one('tip_survival')
rb <- eval_one('unconditional_complete')
rc <- eval_one('legacy_unconditioned')

sim <- simulate_datasets(ra, n_reps=200L, conditioning='from_model', seed=20261118L, num_threads=2L)
write.csv(sim$tip_states, {json.dumps(str(r_sim_csv))}, row.names=FALSE)

bsm <- runBSM(ra, NULL, nummaps_goal=100L, seedval=909L, num_threads=2L, bsm_engine='native_whole_tree')
write.csv(bsm$clado_events_df[, c('map','ancestral_state','left_state','right_state','time_bp','stratum','event_type','source_area','destination_area')],
          {json.dumps(str(r_clado_csv))}, row.names=FALSE)
write.csv(bsm$ana_events_df[, c('map','stratum','from_state','to_state','time_bp','segment_time','event_type','source_area','destination_area')],
          {json.dumps(str(r_ana_csv))}, row.names=FALSE)

phy <- ape::read.tree(file.path(folder, 'tree.newick'))
members <- RustBEARS:::.bgb_u_tree_members(phy)
clades <- vapply(members, function(z) paste(sort(z), collapse='|'), character(1))
payload <- sprintf(
    '{{"lnL_A": %.17g, "surv_A": %.17g, "lnL_B": %.17g, "lnL_C": %.17g, "clades": [%s], "top_A": [%s], "bot_A": [%s]}}',
    ra$total_loglikelihood, ra$log_tip_survival_prob, rb$total_loglikelihood, rc$total_loglikelihood,
    paste(sprintf('"%s"', clades), collapse=','),
    paste(sprintf('%.17g', as.vector(t(ra$ML_marginal_prob_each_state_at_branch_top_AT_node))), collapse=','),
    paste(ifelse(is.na(as.vector(t(ra$ML_marginal_prob_each_state_at_branch_bottom_below_node))), 'null',
                 sprintf('%.17g', as.vector(t(ra$ML_marginal_prob_each_state_at_branch_bottom_below_node)))), collapse=',')
)
writeLines(payload, {json.dumps(str(r_out))})
"""
        script_file = tmp_dir / "run_cross.R"
        script_file.write_text(r_code)
        env_vars = os.environ.copy()
        env_vars["R_LIBS"] = r_libs
        subprocess.run([rscript_bin, str(script_file)], check=True, env=env_vars)

        r_data = json.loads(r_out.read_text())
        r_sim_df = pd.read_csv(r_sim_csv)
        r_clado_df = pd.read_csv(r_clado_csv)
        r_ana_df = pd.read_csv(r_ana_csv)

    from rustbears.uncertainty.averaging import node_clades
    py_tree = py_evals["tip_survival"].model.tree
    py_clades = ["|".join(sorted(c)) for c in node_clades(py_tree)]
    order_r_to_py = [r_data["clades"].index(c) for c in py_clades]
    n_states = len(py_evals["tip_survival"].model.states)
    r_top = np.array(r_data["top_A"], dtype=float).reshape((-1, n_states))[order_r_to_py]
    r_bot = np.array([np.nan if v is None else v for v in r_data["bot_A"]], dtype=float).reshape((-1, n_states))[order_r_to_py]

    py_top = py_evals["tip_survival"].evaluation.ancestral_probs_branch_top
    py_bot = py_evals["tip_survival"].evaluation.ancestral_probs_branch_bottom
    non_root = [i for i in range(len(py_clades)) if i != py_tree.root]

    diff_lnL_A = abs(r_data["lnL_A"] - py_evals["tip_survival"].lnL)
    diff_surv_A = abs(r_data["surv_A"] - py_evals["tip_survival"].evaluation.log_tip_survival_prob)
    diff_lnL_B = abs(r_data["lnL_B"] - py_evals["unconditional_complete"].lnL)
    diff_lnL_C = abs(r_data["lnL_C"] - py_evals["legacy_unconditioned"].lnL)
    diff_top_A = float(np.max(np.abs(r_top - py_top)))
    diff_bot_A = float(np.max(np.abs(r_bot[non_root] - py_bot[non_root])))

    # Align simulated tip columns by taxon label:
    r_sim_aligned = r_sim_df[list(py_sim.tip_labels)].to_numpy(dtype=np.int64)
    sim_exact_match = bool(np.array_equal(r_sim_aligned, py_sim.tip_states))

    # Compare BSM event counts & timestamps:
    py_clado_sub = py_bsm.clado_events_df[["map", "ancestral_state", "left_state", "right_state", "time_bp", "stratum", "event_type"]]
    r_clado_sub = r_clado_df[["map", "ancestral_state", "left_state", "right_state", "time_bp", "stratum", "event_type"]]
    bsm_clado_max_time_diff = float(np.max(np.abs(py_clado_sub["time_bp"].to_numpy() - r_clado_sub["time_bp"].to_numpy())))
    bsm_clado_exact_states = bool(
        (py_clado_sub.drop(columns=["time_bp"]).reset_index(drop=True)
         == r_clado_sub.drop(columns=["time_bp"]).reset_index(drop=True)).all().all()
    )
    py_ana_sub = py_bsm.ana_events_df[["map", "stratum", "from_state", "to_state", "time_bp", "segment_time", "event_type"]]
    r_ana_sub = r_ana_df[["map", "stratum", "from_state", "to_state", "time_bp", "segment_time", "event_type"]]
    bsm_ana_max_time_diff = float(np.max(np.abs(py_ana_sub["time_bp"].to_numpy() - r_ana_sub["time_bp"].to_numpy())))
    bsm_ana_exact_states = bool(
        (py_ana_sub.drop(columns=["time_bp", "segment_time"]).reset_index(drop=True)
         == r_ana_sub.drop(columns=["time_bp", "segment_time"]).reset_index(drop=True)).all().all()
    )

    assert diff_lnL_A < 1e-13 and diff_surv_A < 1e-13 and diff_lnL_B < 1e-13, (diff_lnL_A, diff_surv_A, diff_lnL_B)
    assert diff_top_A < 1e-13 and diff_bot_A < 1e-13, (diff_top_A, diff_bot_A)
    assert sim_exact_match, "R and Python forward simulation matrices differed!"
    assert bsm_clado_exact_states and bsm_clado_max_time_diff < 1e-13
    assert bsm_ana_exact_states and bsm_ana_max_time_diff < 1e-13

    return dict(
        diff_lnL_A=diff_lnL_A,
        diff_log_survival_A=diff_surv_A,
        diff_lnL_B=diff_lnL_B,
        diff_lnL_C_legacy=diff_lnL_C,
        max_abs_diff_ancestral_top_A=diff_top_A,
        max_abs_diff_ancestral_bottom_A=diff_bot_A,
        simulated_200x4_exact_match=sim_exact_match,
        bsm_100maps_clado_exact_states=bsm_clado_exact_states,
        bsm_100maps_clado_max_time_diff=bsm_clado_max_time_diff,
        bsm_100maps_ana_exact_states=bsm_ana_exact_states,
        bsm_100maps_ana_max_time_diff=bsm_ana_max_time_diff,
    )


def run_empirical_hawaii_psychotria_m3(rscript: str, r_libs: str) -> dict:
    bgb_dir = REPO / "r_package/RustBEARS/inst/extdata/examples/Psychotria_M3strat/BGB"
    tree_fn = bgb_dir / "Psychotria_5.2.newick"
    geog_fn = bgb_dir / "Psychotria_geog.data"
    times_fn = bgb_dir / "timeperiods.txt"
    allowed_fn = bgb_dir / "areas_allowed.txt"
    manual_fn = bgb_dir / "dispersal_multipliers.txt"

    families = ("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
    common = dict(
        timeperiods=times_fn,
        areas_allowed=allowed_fn,
        dispersal_multipliers=manual_fn,
    )

    t0 = time.perf_counter()
    fits_c = {}
    fits_a = {}
    for fam in families:
        fits_c[fam] = b.BioGeoBEARS(tree_fn, geog_fn, fam, likelihood_conditioning="legacy_unconditioned", **common).fit(num_threads=4)
        fits_a[fam] = b.BioGeoBEARS(tree_fn, geog_fn, fam, likelihood_conditioning="tip_survival", **common).fit(num_threads=4)
    fit_elapsed = time.perf_counter() - t0

    # Also cross-check R vs Python Contract-A evaluation on the real near-ultrametric Psychotria M3 tree!
    decj_a = fits_a["DEC+J"]
    with tempfile.TemporaryDirectory() as tmp:
        r_json = Path(tmp) / "hawaii_r.json"
        r_script = Path(tmp) / "hawaii_r.R"
        r_script.write_text(f"""
.libPaths(c(strsplit({json.dumps(r_libs)}, .Platform$path.sep)[[1]], .libPaths()))
library(RustBEARS)
x <- define_BioGeoBEARS_run(
    trfn = {json.dumps(str(tree_fn))},
    geogfn = {json.dumps(str(geog_fn))},
    timesfn = {json.dumps(str(times_fn))},
    areas_allowed_fn = {json.dumps(str(allowed_fn))},
    dispersal_multipliers_fn = {json.dumps(str(manual_fn))},
    max_range_size = 4,
    include_null_range = TRUE,
    num_cores_to_use = 1,
    print_optim = FALSE,
    printlevel = -1,
    calc_ancprobs = TRUE
)
x$likelihood_conditioning <- 'tip_survival'
x$allow_null_tips <- FALSE
p <- x$BioGeoBEARS_model_object@params_table
vals <- c(d = {decj_a.params['d']!r}, e = {decj_a.params['e']!r}, j = {decj_a.params['j']!r})
for (k in names(vals)) {{ p[k, 'type'] <- 'fixed'; p[k, c('init', 'est')] <- vals[[k]] }}
p['ysv', 'type'] <- '3-j'; p['ys', 'type'] <- 'ysv*2/3'
p['y', 'type'] <- 'ysv*1/3'; p['s', 'type'] <- 'ysv*1/3'; p['v', 'type'] <- 'ysv*1/3'
x$BioGeoBEARS_model_object@params_table <- p
x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init = TRUE)
invisible(capture.output(x <- readfiles_BioGeoBEARS_run(x)))
invisible(capture.output(x <- section_the_tree(inputs = x, make_master_table = TRUE, plot_pieces = FALSE, cut_fossils = FALSE)))
res <- NULL
invisible(capture.output(res <- bears_optim_run(x, skip_optim = TRUE, skip_optim_option = 'return_all')))
phy <- ape::read.tree(x$trfn)
members <- RustBEARS:::.bgb_u_tree_members(phy)
clades <- vapply(members, function(z) paste(sort(z), collapse = '|'), character(1))
out <- sprintf('{{"lnL": %.17g, "surv": %.17g, "clades": [%s], "top": [%s]}}',
    res$total_loglikelihood, res$log_tip_survival_prob,
    paste(sprintf('"%s"', clades), collapse = ','),
    paste(sprintf('%.17g', as.vector(t(res$ML_marginal_prob_each_state_at_branch_top_AT_node))), collapse = ','))
writeLines(out, {json.dumps(str(r_json))})
""")
        env = os.environ.copy()
        env["R_LIBS"] = r_libs
        subprocess.run([rscript, "--vanilla", str(r_script)], check=True, env=env)
        hawaii_r = json.loads(r_json.read_text())

    from rustbears.uncertainty.averaging import node_clades
    h_clades = ["|".join(sorted(c)) for c in node_clades(decj_a.model.tree)]
    h_order = [hawaii_r["clades"].index(c) for c in h_clades]
    h_r_top = np.array(hawaii_r["top"], dtype=float).reshape((-1, len(decj_a.model.states)))[h_order]
    hawaii_r_py_lnL_diff = abs(hawaii_r["lnL"] - decj_a.lnL)
    hawaii_r_py_surv_diff = abs(hawaii_r["surv"] - decj_a.evaluation.log_tip_survival_prob)
    hawaii_r_py_top_diff = float(np.max(np.abs(h_r_top - decj_a.evaluation.ancestral_probs_branch_top)))
    assert hawaii_r_py_lnL_diff < 1e-13 and hawaii_r_py_surv_diff < 1e-13 and hawaii_r_py_top_diff < 1e-13, (
        hawaii_r_py_lnL_diff, hawaii_r_py_surv_diff, hawaii_r_py_top_diff
    )

    ma_c = b.model_average(list(fits_c.values()), criterion="AICc")
    ma_a = b.model_average(list(fits_a.values()), criterion="AICc")

    # Parametric bootstrap LRT: DEC vs DEC+J on stratified 4-island emergence under Contract-A!
    t_boot = time.perf_counter()
    boot_dec_vs_decj = b.bootstrap_lrt(
        fits_a["DEC"], fits_a["DEC+J"], n_boot=99,
        num_threads=8, seed=20261118,
    )
    boot_elapsed = time.perf_counter() - t_boot

    # Profile CI on stratified DEC+J under Contract-A:
    t_prof = time.perf_counter()
    prof_decj = b.profile_ci(fits_a["DEC+J"], n_grid=12, tol=1e-4, num_threads=8)
    prof_elapsed = time.perf_counter() - t_prof

    # Akaike-averaged 1,000-map stratified BSM + diagnostics under Contract-A:
    t_bsm = time.perf_counter()
    bsm_1000 = ma_a.run_bsm(num_maps=1000, num_threads=8, seed=20261118)
    diag_1000 = b.bsm_diagnostics(bsm_1000, target_rmcse=0.05)
    bsm_elapsed = time.perf_counter() - t_bsm

    # Verify strict geological emergence invariants on all 1,000 empirical Hawaii maps:
    clado = bsm_1000.clado_events_df
    ana = bsm_1000.ana_events_df
    # Stratum 1: 0..0.5 Ma (K,O,M,H); Stratum 2: 0.5..1.9 Ma (K,O,M; no H!);
    # Stratum 3: 1.9..3.7 Ma (K,O; no M,H!); Stratum 4: 3.7..5.2+ Ma (K only!).
    for stratum_id, forbidden in ((2, {"H"}), (3, {"M", "H"}), (4, {"O", "M", "H"})):
        sc = clado[clado["stratum"] == stratum_id]
        for col in ("ancestral_range", "left_range", "right_range"):
            for r in sc[col]:
                assert not any(f in r for f in forbidden), (stratum_id, col, r)
        sa = ana[ana["stratum"] == stratum_id]
        for col in ("from_range", "to_range"):
            for r in sa[col]:
                assert not any(f in r for f in forbidden), (stratum_id, col, r)

    model_rows = []
    for fam in families:
        fc, fa = fits_c[fam], fits_a[fam]
        wc = float(ma_c.weights_df.loc[ma_c.weights_df["model"] == fam, "weight"].iloc[0])
        wa = float(ma_a.weights_df.loc[ma_a.weights_df["model"] == fam, "weight"].iloc[0])
        model_rows.append(dict(
            model=fam,
            lnL_C=fc.lnL,
            d_C=fc.params["d"], e_C=fc.params["e"], j_C=fc.params["j"],
            AICc_weight_C=wc,
            lnL_A=fa.lnL,
            tip_survival_prob_A=float(np.exp(fa.evaluation.log_tip_survival_prob)),
            d_A=fa.params["d"], e_A=fa.params["e"], j_A=fa.params["j"],
            AICc_weight_A=wa,
        ))

    top_routes = (
        diag_1000.route_summary_df.sort_values("mean", ascending=False)
        .head(6)[["source_area", "destination_area", "mean", "se_target", "rmcse"]]
        .to_dict(orient="records")
    )

    return dict(
        dataset="Hawaii Psychotria M3 (19 taxa, 4 islands K->O->M->H monotone emergence, 4 active strata)",
        timings_sec=dict(
            twelve_fits_sec=fit_elapsed,
            bootstrap_lrt_99reps_sec=boot_elapsed,
            profile_ci_decj_sec=prof_elapsed,
            model_averaged_bsm_1000maps_plus_diag_sec=bsm_elapsed,
        ),
        models_comparison=model_rows,
        bootstrap_lrt_dec_vs_decj=dict(
            d_obs=boot_dec_vs_decj.d_obs,
            pval_asymptotic=boot_dec_vs_decj.pval_asymptotic,
            pval_mc_conservative=boot_dec_vs_decj.pval_mc_conservative,
            n_converged=boot_dec_vs_decj.n_converged,
            n_failed=boot_dec_vs_decj.n_failed,
        ),
        profile_ci_decj=prof_decj.to_dataframe().to_dict(orient="records"),
        bsm_1000_geological_invariants="100% zero violations across all 4 geological strata",
        bsm_1000_top_colonization_routes=top_routes,
        r_python_psychotria_m3_parity=dict(
            diff_lnL_A=hawaii_r_py_lnL_diff,
            diff_log_survival_A=hawaii_r_py_surv_diff,
            max_abs_diff_ancestral_top_A=hawaii_r_py_top_diff,
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-libs", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cross = run_cross_language_stratified_parity(args.rscript, args.r_libs)
    hawaii = run_empirical_hawaii_psychotria_m3(args.rscript, args.r_libs)
    report = dict(cross_language_stratified_parity=cross, empirical_hawaii_m3=hawaii, status="PASSED")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
