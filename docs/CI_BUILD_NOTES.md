# Cross-platform build notes

The current validation scope is Linux and macOS, as explicitly requested on
2026-11-15. Both workflow matrices contain only these platforms. Windows build
code is retained, but Windows is outside this verification and is not reported
as passing. The general viewer-process regression remains mandatory.

The workflow defines installed Python wheel tests and the complete R package
check on both selected platforms. The R check fails on errors or warnings; its
additional R scripts exercise both optimizer engines and the v1.2 probability,
statistical and summary APIs. GitHub-hosted execution remains blocked by the
account billing/spending-limit gate. The scope-change commit skips hosted CI;
no payment settings or repository visibility are changed. Equivalent complete
checks are instead recorded separately on the local macOS host and Tongfang
Linux, with their actual runtime, source hash and results. Such records are
local verification, not successful GitHub Actions jobs.

The retained upstream R tests require the optional `rexpokit`, `cladoRcpp`, `FD`,
`snow`, `minqa`, `gdata` and `phylobase` packages; CI installs these explicitly.
They are test dependencies, not replacements for RustBEARS's native backend.
`MultinomialCI` 1.2 is installed from its official CRAN archive because legacy
summary functions retain that dependency.

The retained Windows Rtools build code links the GNU Rust archive with `synchronization` for the Rust
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

The c7 Windows R check reached the end of its full `testthat` suite, then failed
only the new viewer fixture's positive control: unqualified `open` returned 34
instead of the sentinel's 99 and did not write its marker. The actual executable
was not recorded, so that status does not identify its source. The subsequent
real PDF assertions and scientific tests passed, but the separate installed-R
scripts were skipped after the failed package check. This was not a passing
Windows job. R documents that Windows `system()` launches commands directly,
unlike the shell-based Unix interface ([R system documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/system.html)).

The fixture now runs a real noninteractive Rscript and clones the two installed
plotting functions without changing their bodies or arguments. A private child
environment observes only their `system` boundary and forwards it to an
absolute-path Rscript sentinel, which writes a marker and exits 99. Each actual
plotting function first takes a controlled interactive branch, proving that it
reaches that sentinel. The fixture then restores `base::interactive`, deletes
both positive-control PDFs, and requires fresh complete PDFs with no viewer
invocation. Original namespace environments remain unchanged. This tests the
plotting functions' viewer-call contract, not the operating system's default PDF
association. All positive-control and output assertions remain strict.

CI installs the package and runs this quick regression with
`stop_on_failure=TRUE` before the five existing installed-R scripts and full
`R CMD check`. The complete scientific tests and warning/error gate are retained;
only their order changes so a process-fixture problem is reported earlier.
