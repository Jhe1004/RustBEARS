# Cross-platform build notes

The workflow runs installed Python wheels and the complete R package check on
Linux, macOS and Windows. The R check fails on errors or warnings. Its additional
R scripts exercise both optimizer engines and the v1.2 probability, statistical
and summary APIs.

The retained upstream R tests require the optional `rexpokit`, `cladoRcpp`, `FD`,
`snow`, `minqa`, `gdata` and `phylobase` packages; CI installs these explicitly.
They are test dependencies, not replacements for RustBEARS's native backend.
`MultinomialCI` 1.2 is installed from its official CRAN archive because legacy
summary functions retain that dependency.

Windows Rtools links the GNU Rust archive with `synchronization` for the Rust
standard library's `WaitOnAddress` and `WakeByAddress*` synchronization primitives.
The packaged vendor directories are checked out byte-for-byte so that Cargo's
published source checksums remain valid on Windows too.

On all R platforms, the generated Makevars links the shared library before
removing the intermediate `rust/target` directory. This follows the standard
rextendr build ordering (`all: $(SHLIB) rust_clean`) documented in
[rextendr issue 445](https://github.com/extendr/rextendr/issues/445) and
[extendr issue 915](https://github.com/extendr/extendr/issues/915).
R's compiled-code check otherwise attributes symbols from the bundled Rust
standard-library archive to package object files. The cleanup changes no native
code, panic handling or allocation-failure behavior. It does not establish that
all fatal standard-library paths are unreachable; the installed shared library
may retain the standard library's `abort` import. CI preserves the original
check output and does not filter warnings or alter R's symbol-check settings.

The retained `Psychotria_ML_DEC` example has its serialized S4 class package
metadata migrated from `BioGeoBEARS` to `RustBEARS`. The numerical example values
are unchanged. This removes an accidental requirement to install the upstream
package merely to load the bundled data.
