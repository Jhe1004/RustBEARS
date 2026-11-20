# RustBEARS

<p align="left">
  <strong>English</strong> | <a href="README_zh-CN.md">简体中文</a>
</p>

[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Release-v1.2.0-success.svg)](CHANGELOG.md)
[![Rust Core](https://img.shields.io/badge/Rust-1.82%2B-orange.svg)](rust/Cargo.toml)
[![R Package](https://img.shields.io/badge/R-4.0%2B-276DC3.svg)](r_package/RustBEARS/DESCRIPTION)
[![Python Package](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](pyproject.toml)

**RustBEARS** is a memory-safe, parallel **Rust computational engine** for probabilistic historical biogeography, accompanied by a **BioGeoBEARS-compatible R package** (`RustBEARS`) and a **native Python package & CLI** (`rustbears`). It requires no external Fortran 77 (`rexpokit`) or C++ (`cladoRcpp`) package compilation.

---

## Why RustBEARS?

1. **Drop-In BioGeoBEARS Compatibility & High Speed**
   - Retains all **545 public R functions, S4 classes, and plotting utilities** of [`BioGeoBEARS`](http://phylo.wikidot.com/biogeobears) (Matzke 2013, 2014)—existing R scripts run unchanged simply by replacing `library(BioGeoBEARS)` with `library(RustBEARS)`.
   - Reproduces `BioGeoBEARS` log-likelihoods (`median |ΔlnL| = 9.95e-14`) and ancestral probabilities (`<= 1.50e-14`) across 1,339 configurations from 77 published studies, while accelerating 64- to 256-state optimization **2.8–17.3×** on 1 core and **10.9–205.5×** on 8 cores.
2. **Extant Tip-Survival-Conditioned Likelihoods**
   - Corrects the unconditioned null-range survivorship bias that drives DEC local extinction estimates toward zero and inflates support for founder-event speciation (`+J`), across both unstratified phylogenies and multi-stratum models with chronologically emerging islands.
3. **Trial-Free Uniformization Stochastic Mapping**
   - Replaces trial-limited branch-wise rejection sampling with endpoint-conditioned **Uniformization** Markov bridges (**170.9×** median speedup across 30 empirical studies), resolving rare stepping-stone pathways without trial-limit failures and guaranteeing deterministic per-map histories regardless of thread count.
4. **Integrated Frequentist Uncertainty Quantification**
   - Built-in routines for **profile likelihood confidence intervals**, **parametric bootstrap likelihood-ratio tests** (with resumable checkpoints), **Akaike model averaging**, **external posterior-tree propagation**, and **stochastic mapping Monte Carlo precision / MCMC ESS diagnostics**.

---

## Documentation

| Guide | English | 简体中文 |
| :--- | :--- | :--- |
| **Repository Overview & Quick Start** | **README.md** (this page) | [README_zh-CN.md](README_zh-CN.md) |
| **R Package User Guide** (`library(RustBEARS)`) | [docs/R_User_Guide.md](docs/R_User_Guide.md) | [docs/R_User_Guide_zh-CN.md](docs/R_User_Guide_zh-CN.md) |
| **Python API & CLI Guide** (`import rustbears` / `rustbears`) | [docs/Python_and_CLI_Guide.md](docs/Python_and_CLI_Guide.md) | [docs/Python_and_CLI_Guide_zh-CN.md](docs/Python_and_CLI_Guide_zh-CN.md) |
| **Statistical Formulations & Scope** | [docs/Statistical_Methods_and_Contracts.md](docs/Statistical_Methods_and_Contracts.md) | [docs/Statistical_Methods_and_Contracts.md](docs/Statistical_Methods_and_Contracts.md) |
| **Runnable Example Scripts** | [`examples/quickstart_r.R`](examples/quickstart_r.R) · [`examples/uncertainty_workflow.R`](examples/uncertainty_workflow.R) · [`examples/quickstart_python.py`](examples/quickstart_python.py) | 同左 |

---

## Installation

### Step 0: Prerequisite (Rust Toolchain)

Install a standard Rust toolchain (`rustc >= 1.82` and `cargo`) once via [rustup.rs](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

*(All Rust crates are vendored inside the repository, so compilation requires no external Rust network downloads. On Windows with Rtools, add the GNU target once: `rustup target add x86_64-pc-windows-gnu`.)*

---

### Option A: Install the R Package (`RustBEARS`, R >= 4.0)

**Install directly from GitHub in R** (automatically resolves CRAN and archived `MultinomialCI` dependencies):

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("Jhe1004/RustBEARS", subdir = "r_package/RustBEARS")
```

*(Alternatively, using `pak`: `pak::pkg_install("Jhe1004/RustBEARS/r_package/RustBEARS")`.)*

**Or install from a local clone:**

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
R CMD INSTALL r_package/RustBEARS
```

---

### Option B: Install the Python Package & CLI (`rustbears`, Python >= 3.10)

**Install directly from GitHub:**

```bash
python3 -m pip install "git+https://github.com/Jhe1004/RustBEARS.git"
```

**Or install from a local clone:**

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
python3 -m pip install .
```

Verify the command-line tool:

```bash
rustbears fit --help
```

---

## Quick Start

### 1. R Quick Start (Bundled Hawaiian *Psychotria* Dataset)

```r
library(RustBEARS)

# 1. Define a standard BioGeoBEARS run using bundled example files
trfn   <- system.file("extdata/Psychotria_5.2_ultrametric.newick", package = "RustBEARS")
geogfn <- system.file("extdata/Psychotria_geog.data", package = "RustBEARS")

run <- define_BioGeoBEARS_run(
  trfn = trfn,
  geogfn = geogfn,
  max_range_size = 4,
  num_cores_to_use = 1L,
  print_optim = FALSE
)

# 2. Optional: Enable extant tip-survival conditioning (omit for legacy BioGeoBEARS mode)
run$likelihood_conditioning <- "tip_survival"
run$root_prior <- define_root_prior("uniform_nonnull")
run <- readfiles_BioGeoBEARS_run(run)

# 3. Fit DEC and DEC+J
fit_dec <- bears_optim_run(run)

run_j <- run
run_j$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
run_j$BioGeoBEARS_model_object@params_table["j", c("init", "est")] <- 0.01
fit_decj <- bears_optim_run(run_j)

# 4. Profile likelihood 95% confidence intervals & model averaging
ci_dec <- profile_ci(fit_dec, params = c("d", "e"))
print(ci_dec)

avg <- model_average(list(DEC = fit_dec, DECJ = fit_decj), criterion = "AICc", num_maps = 200L, seedval = 42L)
print(avg$weights)
```

👉 See **[docs/R_User_Guide.md](docs/R_User_Guide.md)** for legacy script migration, multi-stratum island emergence, Uniformization `runBSM()`, `bootstrap_lrt()`, and multi-tree propagation `fit_over_trees()`.

---

### 2. Python & CLI Quick Start

```python
import rustbears as rb

tree = "((K_sp1:2.0,O_sp1:2.0):2.0,((O_sp2:1.0,M_sp1:1.0):1.5,(M_sp2:1.0,H_sp1:1.0):1.5):1.5);"
geog = {
    "K_sp1": "1100", "O_sp1": "0100", "O_sp2": "0110",
    "M_sp1": "0010", "M_sp2": "0011", "H_sp1": "0001",
}

# Fit candidate models under extant tip-survival conditioning
models = ["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"]
fits = [
    rb.BioGeoBEARS(
        tree=tree,
        geog=geog,
        model=m,
        max_range_size=4,
        likelihood_conditioning="tip_survival",
        root_prior=rb.RootPriorSpec("uniform_nonnull"),
    ).fit(num_threads=4)
    for m in models
]

print(rb.compare_models(fits))

# Akaike model-averaged Uniformization stochastic mapping & MCSE diagnostics
avg = rb.model_average(fits, criterion="AIC")
maps = avg.run_bsm(num_maps=200, num_threads=4, seed=42)
diag = rb.bsm_diagnostics(maps)
print(maps.summary_counts["event_means"])
```

**Command-line batch execution:**

```bash
rustbears fit \
  --tree tree.newick \
  --geog geography.data \
  --models DEC,DEC+J,DIVALIKE,DIVALIKE+J,BAYAREALIKE,BAYAREALIKE+J \
  --max-range-size 4 \
  --likelihood-conditioning tip_survival \
  --maps 200 \
  --threads 4 \
  --out-dir results_biogeography
```

👉 See **[docs/Python_and_CLI_Guide.md](docs/Python_and_CLI_Guide.md)** for full Python class parameters, stratified models, resumable bootstrap tests, posterior tree propagation, and CLI flags.

---

## Precision Modes & Likelihood Formulations

- **Numerical Precision**: Defaults to 64-bit IEEE-754 double precision (`StrictF64` in R; `precision="strict_f64"` in Python). To reproduce historical 32-bit `cladoRcpp` / `EXPOKIT` rounding behavior bit-for-bit, set `options(RustBEARS.precision = "LegacyF32")` in R or `precision="legacy_f32"` in Python.
- **Likelihood Formulations**:
  1. `"legacy_unconditioned"` *(default for backward compatibility)*: Unconditioned pruning score matching upstream `BioGeoBEARS` (supports fossil tips, non-monotone masks, and trait-coupled models).
  2. `"tip_survival"`: Conditions observed ranges on all sampled extant tips surviving to the present with non-empty ranges (supports DEC, DIVALIKE, BAYAREALIKE, `+J`, spatial modifiers, and multi-stratum chronological area emergence).
  3. `"unconditional_complete"`: Normalized generative probability distribution over all $S^N$ tip outcomes (including null tips) for complete-space forward simulation and probability-sum verification.

---

## Repository Structure

```text
├── rust/                       # Shared Rust engine (bgb-core, bgb-numerics, bgb-clado, bgb-Step, bgb-bsm, bgb-ffi)
├── r_package/RustBEARS/        # BioGeoBEARS-compatible R package (`library(RustBEARS)`)
├── rustbears/                  # Native Python package & CLI (`import rustbears` / `rustbears`)
├── docs/                       # English & Simplified Chinese guides for R, Python, CLI, and methods
├── examples/                   # Ready-to-run R and Python quickstart & uncertainty scripts
├── tests/                      # Rust, R, and Python regression & statistical test suites
├── scripts/                    # Source synchronization, release builder, and benchmark scripts
├── CITATION.cff                # Co-citation metadata
└── LICENSE                     # GNU General Public License v3.0 or later (GPL-3.0-or-later)
```

---

## Citation & License

`RustBEARS` is distributed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** (see [`LICENSE`](LICENSE), [`NOTICE-Rust-port.md`](NOTICE-Rust-port.md), and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)).

When publishing analyses conducted with `RustBEARS`, please cite **both** the foundational `BioGeoBEARS` articles by Nicholas J. Matzke and `RustBEARS` (run `citation("RustBEARS")` in R or see [`CITATION.cff`](CITATION.cff)):

- Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Frontiers of Biogeography*, 5(4), 242–248. https://doi.org/10.21425/F5FBG19694
- Matzke, N. J. (2014). Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. *Systematic Biology*, 63(6), 951–970. https://doi.org/10.1093/sysbio/syu056
- He, J. (2026). RustBEARS: Survival-conditioned likelihoods, Uniformization stochastic mapping, and uncertainty quantification for historical biogeography (v1.2.0). https://github.com/Jhe1004/RustBEARS

Please report any installation, compilation, or runtime questions via [GitHub Issues](https://github.com/Jhe1004/RustBEARS/issues).
