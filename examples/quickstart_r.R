# Quickstart R example using the bundled 19-tip Hawaiian Psychotria dataset
suppressPackageStartupMessages(library(RustBEARS))

trfn   <- system.file("extdata/Psychotria_5.2_ultrametric.newick", package = "RustBEARS")
geogfn <- system.file("extdata/Psychotria_geog.data", package = "RustBEARS")
work   <- tempfile("rustbears-quickstart-")
dir.create(work)

run <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  calc_ancprobs = TRUE,
  num_cores_to_use = 1L,
  print_optim = FALSE,
  printlevel = -1,
  tmpwd = work
)
run$likelihood_conditioning <- "tip_survival"
run$root_prior <- define_root_prior("uniform_nonnull")
run$use_optimx <- FALSE
run <- readfiles_BioGeoBEARS_run(run)
run$wd <- work

fit_dec <- bears_optim_run(run)

run_j <- run
run_j$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
run_j$BioGeoBEARS_model_object@params_table["j", c("init", "est")] <- 0.01
fit_decj <- bears_optim_run(run_j)

avg <- model_average(list(DEC = fit_dec, DECJ = fit_decj), criterion = "AICc", num_maps = 100L, seedval = 42L)
diag <- bsm_diagnostics(avg$joint_bsm, tree_target = "fixed_sample")

cat(sprintf("DEC lnL (tip_survival)   = %.4f\n", fit_dec$total_loglikelihood))
cat(sprintf("DEC+J lnL (tip_survival) = %.4f\n", fit_decj$total_loglikelihood))
cat("Akaike weights:\n")
print(round(avg$weights, 4))
cat("R_QUICKSTART_PASS\n")
