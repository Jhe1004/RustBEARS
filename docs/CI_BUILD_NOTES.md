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

The independent upstream regression oracle is installed after its dependencies.
BioGeoBEARS 1.1.4 retains example data that load `devtools` during installation,
although its DESCRIPTION does not declare that loader dependency. CI installs
`devtools` first, then installs the unmodified pinned BioGeoBEARS source with
no dependency upgrades. The original native oracle dependencies are also pinned
(cladoRcpp `bc5198f075f135eb68913cf987a2f5169021c7ce`, rexpokit
`bfc071ec6d2235a906bbeee122f711f6fea1b78b`). These are optional validation tools;
they do not replace RustBEARS production kernels.

Relative artifact/source-manifest keys use POSIX `/` separators on every OS.
This makes budget fingerprints and exported artifact names portable to Windows;
it does not change absolute runtime paths, numeric results, seeds, or input data.
