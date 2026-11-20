# RustBEARS (R Package)

`RustBEARS` (`1.2.0`) is a drop-in R interface compatible with `BioGeoBEARS` (Matzke 2013, 2014), backed by a memory-safe Rust computational engine requiring no external Fortran (`rexpokit`) or C++ (`cladoRcpp`) package compilation.

## Installation

Ensure a Rust toolchain (`rustc >= 1.82`, `cargo`) is installed from <https://rustup.rs/>, then install from GitHub in R:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Jhe1004/RustBEARS", subdir = "r_package/RustBEARS")
```

Or from a local repository clone:

```bash
R CMD INSTALL r_package/RustBEARS
```

## Quick Usage

Replace `library(BioGeoBEARS)` with `library(RustBEARS)` to accelerate existing scripts unchanged, or enable tip-survival conditioning (`run$likelihood_conditioning <- "tip_survival"`) and uncertainty quantification (`profile_ci`, `bootstrap_lrt`, `model_average`, `fit_over_trees`, `bsm_diagnostics`). See the top-level repository guides ([English](../../docs/R_User_Guide.md) | [简体中文](../../docs/R_User_Guide_zh-CN.md)) for complete tutorials.

When publishing analyses with `RustBEARS`, please run `citation("RustBEARS")` and cite Matzke (2013, 2014) alongside `RustBEARS`.
