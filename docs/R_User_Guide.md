# RustBEARS R Package User Guide

<p align="left">
  <strong>English</strong> | <a href="R_User_Guide_zh-CN.md">简体中文</a> | <a href="../README.md">Back to README</a>
</p>

`RustBEARS` provides a drop-in R interface compatible with `BioGeoBEARS` (all **545** public R functions, S4 classes, and plotting routines exported), powered by a compiled Rust engine requiring no Fortran (`rexpokit`) or C++ (`cladoRcpp`) packages.

---

## 1. Installation & Setup

### Prerequisites
- **R** `>= 4.0.0`
- **Rust toolchain** (`rustc >= 1.82` and `cargo`): install from <https://rustup.rs/>.
  - *Windows note*: Rtools uses a GNU toolchain; add the GNU target once via `rustup target add x86_64-pc-windows-gnu`.

### Install from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("Jhe1004/RustBEARS", subdir = "r_package/RustBEARS")
```

*(Note: `MultinomialCI` was archived on CRAN; `RustBEARS` declares its official CRAN Archive URL in `Remotes`, so `remotes::install_github()` and `pak::pkg_install("Jhe1004/RustBEARS/r_package/RustBEARS")` install it automatically.)*

### Global R Options

```r
library(RustBEARS)

# Numerical precision: "StrictF64" (64-bit IEEE-754 default) or "LegacyF32" (32-bit historical emulation)
options(RustBEARS.precision = "StrictF64")

# Optional thread budget for native Rust kernels
options(RustBEARS.threads = 4L)
```

---

## 2. Drop-In Migration of Existing BioGeoBEARS Scripts

Any existing `BioGeoBEARS` script runs without modification—just replace `library(BioGeoBEARS)` with `library(RustBEARS)`:

```r
library(RustBEARS)

trfn   <- system.file("extdata/Psychotria_5.2_ultrametric.newick", package = "RustBEARS")
geogfn <- system.file("extdata/Psychotria_geog.data", package = "RustBEARS")

BioGeoBEARS_run_object <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  num_cores_to_use = 4L
)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

res_DEC <- bears_optim_run(BioGeoBEARS_run_object)
plot_BioGeoBEARS_results(res_DEC, plotwhat = "pie")
```

By default, `BioGeoBEARS_run_object` uses `likelihood_conditioning = "legacy_unconditioned"` so legacy scripts produce identical unconditioned log-likelihoods and ancestral reconstructions.

---

## 3. Extant Tip-Survival Conditioning & Model Families

To correct null-range survivorship bias on contemporary ultrametric trees, set `likelihood_conditioning <- "tip_survival"` and specify a root prior via `define_root_prior()` before calling `bears_optim_run()`:

```r
run <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  calc_ancprobs = TRUE,
  print_optim = FALSE
)
run$likelihood_conditioning <- "tip_survival"
run$root_prior <- define_root_prior("uniform_nonnull")  # or custom_weights
run <- readfiles_BioGeoBEARS_run(run)

# 1. DEC
fit_dec <- bears_optim_run(run)

# 2. DEC+J
run_decj <- run
run_decj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
run_decj$BioGeoBEARS_model_object@params_table["j", c("init", "est")] <- 0.01
fit_decj <- bears_optim_run(run_decj)

# 3. DIVALIKE (set s = 0, v = free; add j = free for DIVALIKE+J)
run_diva <- run
run_diva$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
run_diva$BioGeoBEARS_model_object@params_table["s", c("init", "est")] <- 0.0
run_diva$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
run_diva$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
run_diva$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"
run_diva$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
run_diva$BioGeoBEARS_model_object@params_table["mx01v", c("init", "est")] <- 0.5
fit_diva <- bears_optim_run(run_diva)

# 4. BAYAREALIKE (set s = 0, v = 0, y = free; add j = free for BAYAREALIKE+J)
run_bay <- run
run_bay$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["s", c("init", "est")] <- 0.0
run_bay$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["v", c("init", "est")] <- 0.0
run_bay$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
run_bay$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
run_bay$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"
run_bay$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
run_bay$BioGeoBEARS_model_object@params_table["mx01y", c("init", "est")] <- 0.9999
fit_bay <- bears_optim_run(run_bay)
```

### Multi-Stratum Time Slicing (Chronologically Emerging Islands)

`tip_survival` directly supports `timeperiods`, `areas_allowed`, `areas_adjacency`, `distances`, and dispersal multipliers whenever areas emerge monotonically toward the present (every area allowed in an older stratum remains allowed in younger strata):

```r
run_strat <- run
run_strat$timesfn <- system.file("extdata/examples/Psychotria_M3strat/BGB/timeperiods.txt", package = "RustBEARS")
run_strat$areas_allowed_fn <- system.file("extdata/examples/Psychotria_M3strat/BGB/areas_allowed.txt", package = "RustBEARS")
run_strat <- readfiles_BioGeoBEARS_run(run_strat)
run_strat <- section_the_tree(inputs = run_strat, make_master_table = TRUE, plot_pieces = FALSE)

fit_strat_dec <- bears_optim_run(run_strat)
```

---

## 4. Fast Endpoint-Conditioned Uniformization Stochastic Mapping

Use `runBSM()` with `bsm_engine = "native_whole_tree"` to sample endpoint-conditioned Uniformization histories without rejection trial-limit failures:

```r
# Fixed map count (e.g., 500 histories)
maps <- runBSM(
  fit_dec, NULL,
  bsm_engine = "native_whole_tree",
  maxnum_maps_to_try = 500L,
  seedval = 42L
)

# Or adaptive precision mode (stops once target relative MCSE is achieved)
maps_adaptive <- runBSM(
  fit_dec, NULL,
  bsm_engine = "native_whole_tree",
  target_rmcse = 0.05,
  min_maps = 200L,
  batch_size = 200L,
  max_maps = 5000L,
  seedval = 42L
)
```

---

## 5. Integrated Uncertainty Quantification Functions

| Function | Purpose | Key Arguments |
| :--- | :--- | :--- |
| `profile_ci(fit, ...)` | Profile likelihood confidence / support intervals with zero-boundary and upper-ridge diagnostics | `params = c("d", "e")`, `level = 0.95`, `n_grid = 16L`, `include_onesided_zero_ref = FALSE` |
| `profile_surface(fit, ...)` | 2D joint profile likelihood surface | `param_x = "d"`, `param_y = "e"`, `grid_x = 12L`, `grid_y = 12L` |
| `bootstrap_lrt(null_fit, alt_fit, ...)` | Parametric bootstrap LRT with null warm starts, conservative MC p-value, and resumable checkpoints | `n_boot = 199L`, `seedval = 42L`, `checkpoint_file = "lrt.rds"`, `resume = FALSE` |
| `model_average(fits, ...)` | AIC/AICc-weighted ancestral probabilities and stratified Uniformization mapping | `criterion = "AICc"`, `num_maps = 500L`, `seedval = 42L` |
| `fit_over_trees(run, trees, ...)` | Refit models and stochastic maps across posterior trees with clade monophyly & missing-fit bounds | `models = c("DEC", "DEC+J")`, `likelihood_conditioning = "tip_survival"`, `bsm_maps_per_tree = 200L`, `checkpoint_file = "trees.rds"` |
| `bsm_diagnostics(bsm_obj, ...)` | Event & route MCSEs, Student-t / Geyer MCMC ESS (`>= 50`) checks, and equivalence tests | `tree_target = "fixed_sample"` (or `"posterior_population"`), `tree_sampling = "independent"` / `"mcmc"`, `equivalence_delta = 0.25` |

### Complete R Uncertainty Example

```r
# 1. Profile likelihood 95% confidence intervals
prof <- profile_ci(fit_dec, params = c("d", "e"), n_grid = 16L)
print(as.data.frame(prof))
plot(prof)

# 2. Parametric bootstrap LRT (supports interruption & resume = TRUE)
lrt <- bootstrap_lrt(
  fit_dec, fit_decj,
  n_boot = 99L,
  seedval = 2026L,
  checkpoint_file = tempfile(fileext = ".rds")
)
print(lrt)

# 3. Akaike model averaging + stochastic mapping diagnostics
avg <- model_average(
  list(DEC = fit_dec, DECJ = fit_decj),
  criterion = "AICc",
  num_maps = 400L,
  seedval = 42L
)
diag <- bsm_diagnostics(avg$joint_bsm, tree_target = "fixed_sample", equivalence_delta = 0.25)
print(diag$event_summary)

# 4. Multi-tree propagation across a posterior tree sample
# `trees` can be a multiPhylo object, a character vector of Newick strings, or a Newick/Nexus file path
trees_demo <- ape::rmtree(5, N = 4) # replace with ape::read.tree("posterior.trees")
# See `examples/uncertainty_workflow.R` for a self-contained runnable script!
```

Run `?profile_ci`, `?bootstrap_lrt`, `?model_average`, `?fit_over_trees`, `?bsm_diagnostics`, or `vignette("uncertainty_workflow", package = "RustBEARS")` inside R for detailed help pages.
