# GPL >= 3. Startup banner enforcing programmatic co-citation and derivative support boundary.
.onAttach <- function(libname, pkgname) {
    if (isTRUE(getOption("RustBEARS.quiet_startup", getOption("BioGeoBEARS.quiet_startup", FALSE)))) {
        return(invisible(NULL))
    }
    packageStartupMessage(
        "RustBEARS v1.2.0 (Unified Rust Engine compatible with BioGeoBEARS; GPL >= 3)\n",
        "  - Foundational methodology & R framework: Nicholas J. Matzke (2013, 2014).\n",
        "  - Co-citation policy: Please cite BOTH Matzke (2013, 2014) and RustBEARS (see citation('RustBEARS')).\n",
        "  - Derivative support boundary: Report technical/build issues to https://github.com/Jhe1004/RustBEARS/issues (not the upstream author)."
    )
}
