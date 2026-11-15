\name{Psychotria_ML_DEC}
\alias{Psychotria_ML_DEC}
\docType{data}
\title{Retained example DEC fit for Hawaiian Psychotria}
\description{An original BioGeoBEARS example result retained with the compatible
R interface. It illustrates the legacy result structure, not a v1.2 calibration study.}
\usage{data(Psychotria_ML_DEC)}
\format{An object of class calc_loglike_sp_results with likelihood, optimizer,
model input/output and node probability components.}
\details{This stored historical object has no v1.2 input fingerprint or explicit
probability-contract metadata. Refit the corresponding bundled geography and
tree before calling the new uncertainty or prepared sampling functions.
The serialized S4 class package references have been migrated to RustBEARS;
all numerical values from the retained example are unchanged.}
\source{Original BioGeoBEARS example data distributed under GPL-compatible terms
with the BioGeoBEARS framework by Nicholas J. Matzke.}
\references{Matzke, N. J. (2014). Model selection in historical biogeography
reveals that founder-event speciation is a crucial process in island clades.
Systematic Biology, 63, 951--970.}
\keyword{datasets}
