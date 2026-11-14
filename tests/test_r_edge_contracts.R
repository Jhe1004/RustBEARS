suppressPackageStartupMessages({
    library(ape)
    library(optimx)
    library(GenSA)
    library(RustBEARS)
})
root <- tempfile("bgb-r-edge-")
dir.create(root)
setwd(root)
writeLines("((T1:1,T2:1):1,T3:2);", "tiny.tree")
writeLines(c("3 3 (A B C)", "T1 100", "T2 010", "T3 001"), "tiny.data")
mask_block <- c("A\tB\tC", "1\t1\t1", "1\t1\t0", "1\t0\t1", "")
mult_block <- c("A\tB\tC", "1\t0.5\t0.2", "0.5\t1\t0.5", "0.2\t0.5\t1", "")
options(BioGeoBEARS.precision = "LegacyF32")
evaluate <- function(times = NULL, as_adjacency = FALSE, optim = FALSE, with_w = FALSE) {
    file <- NA
    if (!is.null(times)) {
        file <- file.path(root, "times.txt")
        writeLines(as.character(times), file)
    }
    mask <- file.path(root, "allowed.txt")
    writeLines(c(rep(mask_block, if (optim && length(times) > 1L) length(times) else 1L), "END"), mask)
    mult <- NA
    if (with_w) {
        mult <- file.path(root, "multipliers.txt")
        writeLines(c(rep(mult_block, if (length(times) > 1L) length(times) else 1L), "END"), mult)
    }
    x <- define_BioGeoBEARS_run(
        trfn = file.path(root, "tiny.tree"),
        geogfn = file.path(root, "tiny.data"),
        timesfn = file,
        max_range_size = 3,
        areas_allowed_fn = if (!as_adjacency) mask else NA,
        areas_adjacency_fn = if (as_adjacency) mask else NA,
        dispersal_multipliers_fn = mult,
        calc_ancprobs = TRUE,
        print_optim = FALSE,
        use_optimx = !optim,
        tmpwd = root
    )
    if (with_w) {
        x$BioGeoBEARS_model_object@params_table["w", "type"] <- "free"
        x$BioGeoBEARS_model_object@params_table["w", "init"] <- 0.8
        x$BioGeoBEARS_model_object@params_table["w", "est"] <- 0.8
    }
    x <- readfiles_BioGeoBEARS_run(x)
    if (length(x$timeperiods) > 1L) x <- section_the_tree(x, make_master_table = TRUE, plot_pieces = FALSE)
    x$wd <- root
    bears_optim_run(x, skip_optim = !optim, skip_optim_option = "return_all")
}
baseline <- evaluate()
stopifnot(ncol(baseline$ML_marginal_prob_each_state_at_branch_top_AT_node) == 7L)
for (times in list(numeric(), 10, c(10, 20, 30), c(2, 10))) {
    other <- evaluate(times)
    stopifnot(abs(other$total_loglikelihood - baseline$total_loglikelihood) <= 1e-12)
}
adj <- evaluate(10, TRUE)
stopifnot(ncol(adj$ML_marginal_prob_each_state_at_branch_top_AT_node) == 6L)

# R2-R-04: +w single-stratum evaluation uses native compiled tree context (0 compatibility fallbacks)
m2_env <- RustBEARS:::.bgb_m2_runtime
before_calls <- m2_env$native_calls
before_fallbacks <- m2_env$compatibility_fallbacks
res_w <- evaluate(with_w = TRUE)
stopifnot(is.finite(res_w$total_loglikelihood))
stopifnot(m2_env$native_calls > before_calls)
stopifnot(m2_env$compatibility_fallbacks == before_fallbacks)

# R2-R-02: get_inputs_for_stochastic_mapping_stratified default cluster_already_open = FALSE works cleanly
writeLines("(((T1:0.4,T2:0.4):0.8,T3:1.2):0.8,T4:2.0);", "strat.tree")
writeLines(c("4 3 (A B C)", "T1 100", "T2 010", "T3 001", "T4 110"), "strat.data")
writeLines(c("0.6", "1.5", "5.0"), "strat.times")
xs <- define_BioGeoBEARS_run(
    trfn = file.path(root, "strat.tree"),
    geogfn = file.path(root, "strat.data"),
    timesfn = file.path(root, "strat.times"),
    max_range_size = 3,
    calc_ancprobs = TRUE,
    print_optim = FALSE,
    tmpwd = root
)
xs <- readfiles_BioGeoBEARS_run(xs)
xs <- section_the_tree(xs, make_master_table = TRUE, plot_pieces = FALSE)
xs$wd <- root
res_strat <- bears_optim_run(xs, skip_optim = TRUE, skip_optim_option = "return_all")
bsm_inputs <- get_inputs_for_stochastic_mapping_stratified(res_strat)
stopifnot(is.list(bsm_inputs), length(bsm_inputs) == 3L)

# R2-R-01: every shipped example R script and 395lab.zip member parses cleanly
ex_root <- system.file("extdata", "examples", package = "RustBEARS")
r_files <- list.files(ex_root, pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
stopifnot(length(r_files) >= 12L)
for (f in r_files) {
    parse(file = f)
}
zip_fn <- file.path(ex_root, "395lab.zip")
if (file.exists(zip_fn)) {
    uz <- tempfile("bgb-zip-check-")
    dir.create(uz)
    utils::unzip(zip_fn, exdir = uz)
    zr_files <- list.files(uz, pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
    stopifnot(length(zr_files) >= 10L)
    for (zf in zr_files) {
        parse(file = zf)
    }
}
cat("R_EDGE_CONTRACTS_PASS\n")
