# RustBEARS (Drop-in R Interface Compatible with BioGeoBEARS)

`RustBEARS` provides a self-contained, memory-safe Rust computational engine (`rustc >= 1.82`, `R >= 4.0.0`, `GPL >= 3`) for probabilistic historical biogeography while keeping all function signatures, S3/S4 return structures, and plotting workflows of `BioGeoBEARS` (Matzke 2013, 2014) 100% compatible.

Named `RustBEARS` to respect and preserve the upstream `BioGeoBEARS` CRAN/GitHub namespace so both packages can coexist cleanly on the same machine.

## Installation

```bash
# On Windows with Rtools, ensure a GNU Rust target is installed once:
# rustup target add x86_64-pc-windows-gnu
R CMD INSTALL r_package/RustBEARS
```

## Usage

Replace `library(BioGeoBEARS)` with `library(RustBEARS)`—all downstream `define_BioGeoBEARS_run()`, `bears_optim_run()`, `runBSM()`, and `plot_BioGeoBEARS_results()` calls run unmodified.

Please run `citation("RustBEARS")` and co-cite both Matzke (2013, 2014) and RustBEARS. Report technical issues to https://github.com/Jhe1004/RustBEARS/issues.
