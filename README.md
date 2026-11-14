# RustBEARS

The `v1.2` branch is a development release (Python `1.2.0.dev0`, R
`1.2.0.9000`). It adds explicit probability contracts, prepared simulation,
profile intervals, parametric bootstrap tests, conditional model/tree summaries,
and mapping precision diagnostics. New workflows support the validated single-layer
six standard model families and distance/manual dispersal modifiers. Time-stratified,
fossil-hook and trait-coupled uncertainty workflows remain explicitly disabled.

R now defaults to `StrictF64`; request `options(RustBEARS.precision="LegacyF32")`
for historical arithmetic. Contract C retains the historical score; A conditions
on nonempty tips and B completes empty-state transmission. Historical scores do
not silently become normalized generating models.

Read the [implementation contract](docs/v1.2-implementation-contract.md),
[Python tutorial](docs/uncertainty-workflow-python.en.md), and
[Chinese Python tutorial](docs/uncertainty-workflow-python.md). Development checks
are separate from the planned large calibration study: neither a working example
nor a CI configuration establishes a completed publication validation.

[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](LICENSE)
[![Rust Core](https://img.shields.io/badge/Rust-1.82%2B-orange.svg)](rust/Cargo.toml)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](pyproject.toml)
[![R](https://img.shields.io/badge/R-4.0%2B-276DC3.svg)](r_package/RustBEARS/DESCRIPTION)

**`RustBEARS`** is a self-contained, memory-safe **Rust computational engine** for the [`BioGeoBEARS`](http://phylo.wikidot.com/biogeobears) probabilistic historical biogeography framework (Matzke 2013, 2014). It provides both a **function-compatible R package** (`RustBEARS`) and a **native Python package & CLI** (`rustbears`) without requiring external Fortran 77 (`rexpokit`/`EXPOKIT`) or C++ (`cladoRcpp`) package compilation.

> **Why `RustBEARS`?** Naming our modernized derivative `RustBEARS` (and `rustbears` in Python) preserves full namespace sovereignty for the original creator's `BioGeoBEARS` packages across CRAN, PyPI, and GitHub, allows `RustBEARS` and legacy `BioGeoBEARS` to be installed side-by-side on the same machine, and keeps 100% of downstream R function signatures (`define_BioGeoBEARS_run`, `bears_optim_run`, `runBSM`, `plot_BioGeoBEARS_results`) drop-in compatible.

---

## Key Capabilities

- **Machine-Precision Equivalence (`10^-14`)**: Supports both **`LegacyF32` (`legacy_f32`)**—which faithfully emulates the historical 32-bit single-precision (`f32`) and decimal step-rounding boundary effects across matrix construction, cladogenetic weights, MaxEnt scaling, and `EXPOKIT` Arnoldi stepping ($9.95 \times 10^{-14}$ median $|\Delta\text{LnL}|$ across 1,339 empirical configurations from 77 published studies)—and **`StrictF64` (`strict_f64`)** pure 64-bit IEEE-754 reference arithmetic.
- **Multi-Basin Cold-Start Optimization**: Tree-scaled interior probes, quadratic infeasible-domain guidance, and scaled absolute-step boundary refinement eliminate flat-penalty stalls on deep-time phylogenies (100% blind non-inferiority across 1,331 empirical fits; rescues 37 suboptimal historical local stalls and solves high-dimensional models up to 848 discrete states).
- **Deterministic Multi-Threaded Stochastic Mapping (BSM)**: Accelerates biogeographical stochastic mapping by a median of **171-fold** single-threaded (up to **2,800-fold**) with an additional **8.59-fold** strong scaling across 8 CPU cores, while guaranteeing **bit-identical per-map event histories** on a given platform across `1, 2, 4, 8, ...` threads via independent four-lane keyed `Xoshiro256PlusPlus` streams derived from `(seed, map_index, stratum, step)`.
- **Zero-Dependency Native Build**: Bundled offline Cargo vendor directory (`rust/vendor`) compiles cleanly with a single modern Rust toolchain across Linux, macOS (`x86_64` and `arm64` Apple Silicon), and Windows.

---

## Repository Structure

```text
├── rust/                       # Unified Rust workspace (bgb-core, bgb-numerics, bgb-clado, bgb-Step, bgb-bsm, bgb-ffi)
├── rustbears/                  # Native Python package & CLI (`import rustbears`)
├── r_package/RustBEARS/        # Function-compatible R package (`library(RustBEARS)`)
├── tests/                      # Self-contained Python & R regression test suites + empirical edge fixtures
├── scripts/                    # Reproducible release synchronization and R source tarball builder
├── CITATION.cff                # Programmatic co-citation metadata (Matzke 2013, 2014 + RustBEARS)
├── NOTICE-Rust-port.md         # Upstream copyright chain, component license matrix, and provenance
├── THIRD_PARTY_LICENSES.md     # Full license texts for vendored MIT crates and reference BLAS/LAPACK
└── LICENSE                     # GNU General Public License v3.0 or later (GPL-3.0-or-later)
```

---

## Installation

### Prerequisites
- **Rust toolchain** (`cargo` and `rustc >= 1.82`): install via [rustup.rs](https://rustup.rs/).
  - *Windows R package note*: Rtools uses a MinGW/GNU linker; if your default Windows Rust host is MSVC, add the GNU target once via `rustup target add x86_64-pc-windows-gnu`.

### 1. Python Package (`rustbears`)

Install into an existing scientific Python (`>= 3.10`) environment containing `numpy`, `scipy`, `pandas`, and `matplotlib`:

```bash
python -m pip install numpy scipy pandas matplotlib
python -m pip install --no-build-isolation --no-deps .
```

### 2. R Package (`RustBEARS`)

Requires R `>= 4.0.0`. No `rexpokit` or `cladoRcpp` Fortran/C++ compilation is required. Build and install directly from `r_package/RustBEARS`:

```bash
R CMD INSTALL r_package/RustBEARS
```

Or generate a reproducible source tarball first via `./scripts/build_release_artifacts.sh` and install `dist/RustBEARS_1.1.4.tar.gz`.

---

## Quick Start

### Python API & CLI

```python
from rustbears import RustBEARS, compare_models

# Fit standard or +J models (precision="legacy_f32" or "strict_f64")
models = ["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"]
fits = [
    RustBEARS("tree.newick", "geography.data", model=m, max_range_size=4, precision="legacy_f32").fit()
    for m in models
]
print(compare_models(fits))

# Deterministic multi-threaded Biogeographical Stochastic Mapping (200+ maps in sub-second runtime)
best_fit = fits[1]
maps = best_fit.run_bsm(num_maps=200, num_threads=8, seed=12345)
print(maps.summary_counts["event_means"])
```

Command-line batch execution:

```bash
rustbears fit --tree tree.newick --geog geography.data \
  --models DEC,DEC+J,DIVALIKE,DIVALIKE+J,BAYAREALIKE,BAYAREALIKE+J \
  --max-range-size 4 --out-dir results
```

### R Compatibility Layer (`library(RustBEARS)`)

Replace `library(BioGeoBEARS)` with `library(RustBEARS)`—all downstream `BioGeoBEARS` modelling, BSM, and plotting functions run **without modification**:

```r
library(RustBEARS)

# Switch between legacy rounding emulation ("LegacyF32", default) and pure 64-bit ("StrictF64")
options(RustBEARS.precision = "LegacyF32")

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- "tree.newick"
BioGeoBEARS_run_object$geogfn <- "geography.data"
BioGeoBEARS_run_object$max_range_size <- 4

res <- bears_optim_run(BioGeoBEARS_run_object)
```

### Frontend Interface Conventions (R vs. Python)

- **Ambiguous Tip Geography (`?`)**: Both frontends reject or ignore unexpanded `?` wildcards by default and expand partial/complete `?` tip encodings identically when explicitly enabled (`BioGeoBEARS_run_object$useAmbiguities <- TRUE` in R; `use_ambiguities=True` or `--use-ambiguities` in Python).
- **Fossil Hook / Zero-Length Branches**: The shared Rust likelihood kernel routes any branch shorter than `min_branchlength` (default `1e-6`) as a direct-ancestor fossil hook. Python (`RustBEARS`) accepts exact `0.0` branch lengths in Newick strings directly, whereas R's top-level `check_BioGeoBEARS_run()` preserves the historical upstream guard against `<= 0` edge lengths (use a positive hook length `< min_branchlength`, such as `1e-8`, in top-level R runs).

---

## Verification & Regression Tests

Run the full release verification script (or execute the Rust and Python tests directly):

```bash
./scripts/build_release_artifacts.sh
```

---

## Community Governance, License & Scholarly Co-Citation

1. **Open-Source License Inheritance (`GPL-3.0-or-later`)**: Newly authored Rust crates are licensed under `GPL-2.0-or-later`; because the combined distribution retains upstream `BioGeoBEARS` R scripts carrying both `GPL (>= 2)` and `GPL-3` / `GPL (>= 3)` headers, the combined repository and packages are distributed under **`GPL-3.0-or-later`** (see [`LICENSE`](LICENSE), [`NOTICE-Rust-port.md`](NOTICE-Rust-port.md), and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)).
2. **Scholarly Co-Citation**: As a matter of academic attribution, when publishing analyses performed with `RustBEARS`, **please cite both the foundational `BioGeoBEARS` publications by Nicholas J. Matzke and `RustBEARS`** (see `citation("RustBEARS")` in R or [`CITATION.cff`](CITATION.cff)):
   - Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Frontiers of Biogeography*, 5(4), 242–248. https://doi.org/10.21425/F5FBG19694
   - Matzke, N. J. (2014). Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. *Systematic Biology*, 63(6), 951–970. https://doi.org/10.1093/sysbio/syu056
   - He, J. (2026). RustBEARS: Unified Rust computational engine for probabilistic historical biogeography (v1.1.4). https://github.com/Jhe1004/RustBEARS
3. **Explicit Derivative Naming & Support Boundary**: `RustBEARS` is an independent modernized derivative maintained by the refactoring authors. **Please report all installation, compilation, or runtime issues regarding `RustBEARS` to [GitHub Issues](https://github.com/Jhe1004/RustBEARS/issues)** so that the original upstream author is not burdened with derivative support requests.
