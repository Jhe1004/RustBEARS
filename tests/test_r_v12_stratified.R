suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64')

folder <- tempfile('rustbears-v12-strat-r-')
dir.create(folder)
on.exit(unlink(folder, recursive = TRUE), add = TRUE)

# 1. Exact 3-stratum monotone island emergence + adjacency + distances (+x) + env (+n) + manual (+w) + area (+u)
tree_txt <- '((T1:0.31234567,T2:0.31234567):0.48765433,(T3:0.52345678,T4:0.52345678):0.27654322);'
writeLines(tree_txt, file.path(folder, 'tree.newick'))
writeLines(c('4 3 (A B C)', 'T1 100', 'T2 010', 'T3 001', 'T4 110'), file.path(folder, 'geog.data'))
writeLines(c('0.25', '0.65', '2.0'), file.path(folder, 'times.txt'))

write_mat_blocks <- function(mats, path) {
    lines <- character()
    for (m in mats) {
        lines <- c(lines, 'A\tB\tC')
        for (r in seq_len(nrow(m))) lines <- c(lines, paste(m[r, ], collapse = '\t'))
        lines <- c(lines, '')
    }
    lines <- c(lines, 'END')
    writeLines(lines, path)
}

allowed_mats <- list(
    matrix(c(1,1,1, 1,1,1, 1,1,1), 3, 3, byrow = TRUE),
    matrix(c(1,1,0, 1,1,0, 0,0,0), 3, 3, byrow = TRUE),
    matrix(c(1,0,0, 0,0,0, 0,0,0), 3, 3, byrow = TRUE)
)
adj_mats <- replicate(3, matrix(c(1,1,0, 1,1,1, 0,1,1), 3, 3, byrow = TRUE), simplify = FALSE)
dist_mats <- list(
    matrix(c(0, 2.0, 4.5, 2.0, 0, 2.2, 4.5, 2.2, 0), 3, 3, byrow = TRUE),
    matrix(c(0, 2.6, 5.0, 2.6, 0, 2.8, 5.0, 2.8, 0), 3, 3, byrow = TRUE),
    matrix(c(0, 3.2, 6.0, 3.2, 0, 3.5, 6.0, 3.5, 0), 3, 3, byrow = TRUE)
)
env_mats <- list(
    matrix(c(0, 1.3, 2.1, 1.3, 0, 1.5, 2.1, 1.5, 0), 3, 3, byrow = TRUE),
    matrix(c(0, 1.6, 2.4, 1.6, 0, 1.7, 2.4, 1.7, 0), 3, 3, byrow = TRUE),
    matrix(c(0, 1.9, 2.8, 1.9, 0, 2.0, 2.8, 2.0, 0), 3, 3, byrow = TRUE)
)
manual_mats <- list(
    matrix(c(1.0, 0.85, 0.35, 0.90, 1.0, 0.75, 0.25, 0.65, 1.0), 3, 3, byrow = TRUE),
    matrix(c(1.0, 0.55, 0.05, 0.50, 1.0, 0.05, 0.05, 0.05, 1.0), 3, 3, byrow = TRUE),
    matrix(c(1.0, 0.05, 0.05, 0.05, 1.0, 0.05, 0.05, 0.05, 1.0), 3, 3, byrow = TRUE)
)
write_mat_blocks(allowed_mats, file.path(folder, 'allowed.txt'))
write_mat_blocks(adj_mats, file.path(folder, 'adj.txt'))
write_mat_blocks(dist_mats, file.path(folder, 'dist.txt'))
write_mat_blocks(env_mats, file.path(folder, 'env.txt'))
write_mat_blocks(manual_mats, file.path(folder, 'manual.txt'))
writeLines(c('A\tB\tC', '1.6\t1.1\t0.7', '', 'A\tB\tC', '1.4\t0.8\t0.5', '', 'A\tB\tC', '1.2\t0.5\t0.5', '', 'END'),
           file.path(folder, 'areas.txt'))

make_strat_run <- function(contract, fixed_only = TRUE, j_free = FALSE, x_free = FALSE, allowed_file = 'allowed.txt') {
    x <- define_BioGeoBEARS_run(
        trfn = file.path(folder, 'tree.newick'),
        geogfn = file.path(folder, 'geog.data'),
        timesfn = file.path(folder, 'times.txt'),
        areas_allowed_fn = file.path(folder, allowed_file),
        areas_adjacency_fn = file.path(folder, 'adj.txt'),
        distsfn = file.path(folder, 'dist.txt'),
        dispersal_multipliers_fn = file.path(folder, 'manual.txt'),
        area_of_areas_fn = file.path(folder, 'areas.txt'),
        max_range_size = 3,
        print_optim = FALSE,
        printlevel = -1,
        tmpwd = folder,
        calc_ancprobs = TRUE
    )
    x$envdistsfn <- file.path(folder, 'env.txt')
    x$wd <- folder
    x$likelihood_conditioning <- contract
    x$allow_null_tips <- contract == 'unconditional_complete'
    p <- x$BioGeoBEARS_model_object@params_table
    vals <- c(d = 0.18, e = 0.11, j = 0.14, x = -0.4, w = 1.2, n = -0.3, u = -0.5)
    for (k in names(vals)) {
        p[k, 'type'] <- if (fixed_only) 'fixed' else if (k %in% c('d', 'e') || (k == 'j' && j_free) || (k == 'x' && x_free)) 'free' else 'fixed'
        p[k, c('init', 'est')] <- vals[[k]]
    }
    if (!fixed_only && !j_free) p['j', c('init', 'est')] <- 0
    if (j_free || fixed_only) {
        p['ysv', 'type'] <- '3-j'
        p['ys', 'type'] <- 'ysv*2/3'
        p['y', 'type'] <- 'ysv*1/3'
        p['s', 'type'] <- 'ysv*1/3'
        p['v', 'type'] <- 'ysv*1/3'
    }
    x$BioGeoBEARS_model_object@params_table <- p
    x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object, update_init = TRUE)
    invisible(utils::capture.output(x <- readfiles_BioGeoBEARS_run(x)))
    invisible(utils::capture.output(x <- section_the_tree(inputs = x, make_master_table = TRUE, plot_pieces = FALSE, cut_fossils = FALSE)))
    x
}

fit_a <- NULL; fit_b <- NULL; fit_c <- NULL
invisible(utils::capture.output(fit_a <- bears_optim_run(make_strat_run('tip_survival', fixed_only = TRUE), skip_optim = TRUE, skip_optim_option = 'return_all')))
invisible(utils::capture.output(fit_b <- bears_optim_run(make_strat_run('unconditional_complete', fixed_only = TRUE), skip_optim = TRUE, skip_optim_option = 'return_all')))
invisible(utils::capture.output(fit_c <- bears_optim_run(make_strat_run('legacy_unconditioned', fixed_only = TRUE), skip_optim = TRUE, skip_optim_option = 'return_all')))

stopifnot(is.finite(fit_a$total_loglikelihood), is.finite(fit_b$total_loglikelihood), is.finite(fit_c$total_loglikelihood))
stopifnot(is.finite(fit_a$log_tip_survival_prob))
stopifnot(abs(fit_a$total_loglikelihood - (fit_b$total_loglikelihood - fit_a$log_tip_survival_prob)) < 1e-13)
# Root prior in oldest stratum (only A allowed) must put probability 1.0 on state A (index 2 in 1-based R):
stopifnot(identical(fit_a$resolved_root_prior$sampling_probs, c(0, 1, 0, 0, 0, 0, 0, 0)))

# Exhaustive pattern sum over all allowed present-day states in R:
phy <- ape::read.tree(file.path(folder, 'tree.newick'))
states <- RustBEARS:::rcpp_areas_list_to_states_list(areas = c('A', 'B', 'C'), maxareas = 3, include_null_range = TRUE)
mask0 <- RustBEARS:::.bgb_v12_layer_mask(fit_b$inputs, states, 1L, TRUE)
allowed_b <- which(mask0) - 1L
allowed_a <- setdiff(allowed_b, 0L)

# Simulation duality & stratified BSM geological invariants in R:
sim_a <- simulate_datasets(fit_a, n_reps = 2000L, conditioning = 'from_model', seed = 20261118L, num_threads = 2L, return_node_states = TRUE)
stopifnot(all(sim_a$tip_states %in% allowed_a))

bsm_a <- runBSM(fit_a, NULL, nummaps_goal = 300L, seedval = 909L, num_threads = 2L, bsm_engine = 'native_whole_tree')
s3_clado <- bsm_a$clado_events_df[bsm_a$clado_events_df$stratum == 3L, ]
stopifnot(all(s3_clado$ancestral_range == 'A'), all(s3_clado$left_range == 'A'), all(s3_clado$right_range == 'A'))
s2_clado <- bsm_a$clado_events_df[bsm_a$clado_events_df$stratum == 2L, ]
stopifnot(!any(grepl('C', s2_clado$ancestral_range)), !any(grepl('C', s2_clado$left_range)), !any(grepl('C', s2_clado$right_range)))
stopifnot(!('AC' %in% bsm_a$clado_events_df$ancestral_range), !('AC' %in% bsm_a$ana_events_df$to_range))

# Non-monotone area sinking must fail explicitly:
sinking_mats <- list(
    matrix(c(1,0,0, 0,0,0, 0,0,1), 3, 3, byrow = TRUE),
    matrix(c(1,1,0, 1,1,0, 0,0,1), 3, 3, byrow = TRUE),
    matrix(c(1,1,0, 1,1,0, 0,0,0), 3, 3, byrow = TRUE)
)
write_mat_blocks(sinking_mats, file.path(folder, 'sinking.txt'))
bad_run <- make_strat_run('tip_survival', fixed_only = TRUE, allowed_file = 'sinking.txt')
stopifnot(inherits(try(bears_optim_run(bad_run, skip_optim = TRUE, skip_optim_option = 'return_all'), silent = TRUE), 'try-error'))

# End-to-end R uncertainty workflows on stratified island-emergence models:
fit_null <- NULL; fit_alt <- NULL
invisible(utils::capture.output(fit_null <- bears_optim_run(make_strat_run('tip_survival', fixed_only = FALSE, j_free = FALSE), skip_optim = FALSE)))
invisible(utils::capture.output(fit_alt <- bears_optim_run(make_strat_run('tip_survival', fixed_only = FALSE, j_free = TRUE, x_free = TRUE), skip_optim = FALSE)))
stopifnot(is.finite(fit_null$total_loglikelihood), is.finite(fit_alt$total_loglikelihood))

prof <- profile_ci(fit_alt, params = c('j', 'x'), n_grid = 5L, tol = 1e-3, num_threads = 2L)
prof_df <- as.data.frame(prof)
stopifnot(identical(prof_df$param, c('j', 'x')), all(prof_df$endpoint_deviance_residual_max < 1e-3))

boot <- bootstrap_lrt(fit_null, fit_alt, n_boot = 4L, num_threads = 2L, seedval = 42L)
stopifnot(boot$n_planned == 4L, is.finite(boot$d_obs))

ma <- model_average(list(DEC = fit_null, DEC_Jx = fit_alt), criterion = 'AIC', num_maps = 20L, num_threads = 2L, seedval = 55L)
stopifnot(abs(sum(ma$weights) - 1) < 1e-12, nrow(ma$joint_bsm$clado_events_df) == 20L * 3L)

trees <- list(
    ape::read.tree(text = tree_txt),
    ape::read.tree(text = '((T1:0.29,T2:0.29):0.51,(T3:0.50,T4:0.50):0.30);'),
    ape::read.tree(text = '((T2:0.32,T1:0.32):0.49,(T3:0.51,T4:0.51):0.30);')
)
mt_template <- make_strat_run('tip_survival', fixed_only = FALSE, j_free = TRUE)
mt <- fit_over_trees(mt_template, trees, models = c('DEC', 'DEC+J'), likelihood_conditioning = 'tip_survival',
                     criterion = 'AIC', bsm_maps_per_tree = 8L, num_threads = 2L, seedval = 88L)
stopifnot(mt$n_trees_completed == 3L)
diag_mt <- bsm_diagnostics(mt$joint_bsm, tree_target = 'posterior_population', tree_dependence = 'independent')
stopifnot(diag_mt$num_maps_total == 24L, nrow(diag_mt$event_summary_df) > 0L)

cat('R_V12_STRATIFIED_PASS\n')
