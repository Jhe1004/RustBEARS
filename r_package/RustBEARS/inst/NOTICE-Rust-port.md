# Attribution, Provenance, and Support Boundary Notice (`RustBEARS`)

## 1. Open-Source License Inheritance and Component License Matrix

The combined `RustBEARS` repository, Python distribution (`rustbears`), and drop-in R package (`RustBEARS`) are distributed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`, i.e. `GPL >= 3` in R)** in accordance with Section 9 of GPLv2:

- **Newly authored Rust workspace crates (`rust/crates/bgb-*`)**: Authored under `GPL-2.0-or-later` so pure-Rust components remain compatible with both `GPLv2+` and `GPLv3+` downstream ecosystems.
- **Upstream `BioGeoBEARS`, `cladoRcpp`, and `rexpokit` (`GPL >= 2` and `GPL-3` / `GPL >= 3`)**: Original foundational statistical methodology, R modelling/plotting scripts, state-ordering conventions, and cladogenetic inheritance formulations by **Nicholas J. Matzke** (2013, *Frontiers of Biogeography*; 2014, *Systematic Biology*) and contributors credited in the upstream R packages. Because specific retained upstream R files (`aa_generics_v1.R`, `BioGeoBEARS-package.R`) carry explicit `GPL-3` / `GPL (>= 3)` headers alongside `GPL (>= 2)` files, the combined work is distributed under **`GPL-3.0-or-later`** while preserving every upstream file's original copyright and license header intact.
- **`EXPOKIT`**: The legacy degree-six Padé and Krylov subspace matrix exponential formulation (`LegacyF32` compatibility path) derives from **Roger B. Sidje**, *EXPOKIT: A software package for computing matrix exponentials* (*ACM Transactions on Mathematical Software*, 24(1), 130–156, 1998), as distributed with permission in GPL `rexpokit`.
- **Maximum-Entropy Iterative Scaling**: Discrete maximum-entropy solver written from scratch in `bgb-core` (`GPL-2.0-or-later`) and wrapped in R via a newly authored `zzz_bgb_maxent.R` (`GPL-3.0-or-later`) directly from the published Improved Iterative Scaling (IIS) equations of **Stephen Della Pietra, Vincent Della Pietra & John Lafferty** (1997, *IEEE TPAMI*, 19(4), 380–393). To avoid any `GPL-2`-only vs `GPL-3.0-or-later` license ambiguity from code originally bundled in the CRAN package `FD` (Laliberté & Legendre 2010; Shipley 2010) and copied into `rexpokit`, **no R or Fortran source lines from `FD` (`maxent.R` / `itscale5.f`) are included in `RustBEARS`**, while preserving full scholarly credit and public API compatibility.
- **Vendored Rust Crates & BLAS/LAPACK Notices**: Permissively licensed Rust crates statically compiled into native binaries (`rayon`, `rayon-core`, `either`, `crossbeam-utils`, `crossbeam-epoch`, `crossbeam-deque` under MIT) and reference BLAS/LAPACK loop conventions retain their full copyright notices and license texts in the accompanying [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) shipped in the repository root, Python wheels/sdists, and R package `inst/` directory.
- **64-bit Reference Matrix Exponential Action (`StrictF64`)**: Independently implemented in Rust from **Awad H. Al-Mohy & Nicholas J. Higham** (2009, DOI `10.1137/09074721X`; 2011, *SIAM Journal on Scientific Computing*, 33(2), 488–511, DOI `10.1137/100788860`). No legacy Fortran or C++ object code from `rexpokit` or `cladoRcpp` is linked at compile or run time.

## 2. Scholarly Co-Citation Request

As a matter of scholarly academic practice (and **not** as an additional legal restriction on rights granted under the GNU GPL), when publishing empirical or methodological research using `RustBEARS`—either the drop-in R package (`RustBEARS`) or the native Python package (`rustbears`)—please **co-cite both** the foundational `BioGeoBEARS` publications and `RustBEARS`:

1. Matzke, N. J. Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Front. Biogeogr.* **5**, 242–248 (2013).
2. Matzke, N. J. Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. *Syst. Biol.* **63**, 951–970 (2014).
3. He, J. RustBEARS: Unified Rust computational engine for probabilistic historical biogeography (v1.1.4, 2026). https://github.com/Jhe1004/RustBEARS

## 3. Derivative Naming, Namespace Respect, and Support Boundary

This project is named **`RustBEARS`** (R package `RustBEARS`; Python package `rustbears`) so that the upstream `BioGeoBEARS` package namespaces on CRAN, PyPI, and GitHub remain completely unencumbered for the original creator, while allowing `RustBEARS` and legacy `BioGeoBEARS` to coexist side-by-side. **To shield the original upstream author from derivative maintenance burdens**, please direct all build, installation, numerical, or API questions regarding `RustBEARS` exclusively to our issue tracker:

- **Issue Tracker**: https://github.com/Jhe1004/RustBEARS/issues
