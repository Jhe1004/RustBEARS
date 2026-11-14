suppressPackageStartupMessages({
    library(ape)
    library(optimx)
    library(GenSA)
    library(RustBEARS)
})
folder <- tempfile("bgb-r-opt-")
dir.create(folder)
setwd(folder)
writeLines("((A:150,B:150):150,C:300);", "tree.tree")
writeLines(c("3 3 (X Y Z)", "A 100", "B 010", "C 001"), "geog.data")
options(RustBEARS.precision = "LegacyF32", BioGeoBEARS.precision = NULL)

# Verify proposal search genuinely enters StrictF64 (1L) even when RustBEARS.precision="LegacyF32"
seen_precisions <- integer()
probe_fn <- function(p) {
    seen_precisions <<- c(seen_precisions, RustBEARS:::.bgb_precision())
    -(sum((p - 0.02)^2))
}
dummy <- RustBEARS:::.bgb_composite_search(c(d = 0.01, e = 0.01), probe_fn, list(), c(1e-6, 1e-6), c(1, 1), maxit = 5L)
stopifnot(
    length(seen_precisions) > 0L,
    all(seen_precisions == 1L),
    identical(getOption("RustBEARS.precision"), "LegacyF32"),
    is.null(getOption("BioGeoBEARS.precision"))
)

results <- list()
for (strata in c(FALSE, TRUE)) {
    for (optimx in c(FALSE, TRUE)) {
        if (strata) writeLines(c("100", "400"), "time.txt")
        x <- define_BioGeoBEARS_run(
            trfn = "tree.tree",
            geogfn = "geog.data",
            timesfn = if (strata) "time.txt" else NA,
            max_range_size = 3,
            num_cores_to_use = 2,
            use_optimx = optimx,
            print_optim = FALSE,
            calc_ancprobs = TRUE,
            tmpwd = folder
        )
        x <- readfiles_BioGeoBEARS_run(x)
        if (strata) x <- section_the_tree(x, make_master_table = TRUE, plot_pieces = FALSE)
        x$wd <- folder
        y <- bears_optim_run(x, skip_optim = FALSE)
        stopifnot(
            is.finite(y$total_loglikelihood),
            is.matrix(y$ML_marginal_prob_each_state_at_branch_top_AT_node),
            !is.null(y$outputs@params_table),
            identical(getOption("RustBEARS.precision"), "LegacyF32")
        )
        if (optimx) {
            stopifnot(
                inherits(y$optim_result, "optimx"),
                all(c("p1", "p2", "value", "fevals", "convcode") %in% colnames(y$optim_result))
            )
        } else {
            stopifnot(is.null(names(y$optim_result$par)), grepl("StrictF64 proposal", y$optim_result$message, fixed = TRUE))
        }
        results[[paste(strata, optimx)]] <- y$total_loglikelihood
    }
}
stopifnot(
    abs(results[["FALSE FALSE"]] - results[["FALSE TRUE"]]) < 1e-4,
    abs(results[["TRUE FALSE"]] - results[["TRUE TRUE"]]) < 1e-4
)
cat("R_OPTIMIZER_AND_THREADS_PASS\n")
