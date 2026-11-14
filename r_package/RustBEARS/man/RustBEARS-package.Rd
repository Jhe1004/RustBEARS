\name{RustBEARS-package}
\alias{RustBEARS-package}
\alias{RustBEARS}
\docType{package}
\title{
RustBEARS: Rust-Accelerated Probabilistic Historical Biogeography Compatible with 'BioGeoBEARS'
}
\description{
Drop-in function-compatible R interface to the \code{BioGeoBEARS} historical biogeography framework powered by \code{RustBEARS}, a self-contained, memory-safe Rust computational engine (zero external Fortran/C++ package compilation required). Reproduces probabilistic historical biogeography models (\code{DEC}, \code{DIVALIKE}, \code{BAYAREALIKE}, and founder-event speciation \code{+J}) and biogeographical stochastic mapping (BSM) to machine precision while accelerating likelihood optimization and Monte Carlo mapping.
}
\details{
Named \code{RustBEARS} to preserve the upstream \code{BioGeoBEARS} package namespace while keeping all modelling, optimization, stochastic mapping, and plotting function signatures (\code{define_BioGeoBEARS_run}, \code{bears_optim_run}, \code{calc_loglike_sp}, \code{runBSM}, \code{plot_BioGeoBEARS_results}) 100\% compatible.

Please co-cite BOTH the foundational \code{BioGeoBEARS} publications by Nicholas J. Matzke (2013, 2014) and \code{RustBEARS}; run \code{citation("RustBEARS")} for complete bibliographic entries.
}
\author{
Nicholas J. Matzke (foundational \code{BioGeoBEARS} methodology and R modelling framework),
Jian He (unified Rust computational engine, drop-in R bridge, and derivative package maintenance <j.he930724@gmail.com>).
}
\references{
Matzke, N. J. (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. \emph{Frontiers of Biogeography}, 5(4), 242--248.

Matzke, N. J. (2014). Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. \emph{Systematic Biology}, 63(6), 951--970.

He, J. (2026). Modernizing legacy scientific software with verifiable AI agents: a case study on BioGeoBEARS. \url{https://github.com/Jhe1004/RustBEARS}
}
\keyword{package}
