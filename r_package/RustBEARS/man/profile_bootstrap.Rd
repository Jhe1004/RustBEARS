\name{profile_ci}
\alias{profile_ci}
\alias{profile_surface}
\alias{bootstrap_lrt}
\title{Native profile support and parametric bootstrap inference}
\description{Reoptimizes nuisance parameters and calibrates nested model score
differences using the native Rust probability contracts, without a Python bridge.}
\usage{
profile_ci(res, params = NULL, level = 0.95,
  method = c("asymptotic_chi2", "bootstrap_plugin"), n_boot_plugin = 199L,
  include_onesided_zero_ref = FALSE, eval_exact_zero = TRUE,
  multi_basin_guard = TRUE, n_grid = 16L, tol = 1e-4,
  num_threads = getOption("RustBEARS.threads", 1L), seedval = 12345L,
  simulation_conditioning = "from_model", checkpoint_file = NULL,
  resume = FALSE, maxit = 500L)
profile_surface(res, param_x = "d", param_y = "j", grid_size = 25L,
  reoptimize_nuisance = TRUE,
  num_threads = getOption("RustBEARS.threads", 1L),
  grid_x = NULL, grid_y = NULL, maxit = 500L)
bootstrap_lrt(res_null, res_alt, n_boot = 399L,
  target_conditioning = c("match_fits", "tip_survival",
    "unconditional_complete", "legacy_unconditioned"),
  refit_observed = FALSE,
  root_distribution = c("model_prior", "empirical_posterior"),
  calc_surrogate_power = FALSE, n_alt_sim = 199L, alpha = 0.05,
  num_threads = getOption("RustBEARS.threads", 1L), seedval = 12345L,
  checkpoint_file = NULL, simulation_conditioning = "from_model",
  resume = FALSE, maxit = 500L)
}
\arguments{
\item{res, res_null, res_alt}{Finite fitted v1.2 results, with unchanged input
files and explicit probability-contract metadata. Null and alternative models
must share observations, root weights, precision, and likelihood contract.}
\item{params}{Distinct free parameters to profile; defaults to all free parameters.}
\item{level}{Profile reference probability in (0,1).}
\item{method}{Chi-square reference or local bootstrap calibration at the fitted maximum.}
\item{n_boot_plugin, n_boot, n_alt_sim}{Planned replicate counts; failures are never silently discarded.}
\item{include_onesided_zero_ref}{Also report the chi-square reference at probability
\code{2*level-1} when exact physical zero is feasible.}
\item{eval_exact_zero}{Evaluate exact zero for nonnegative rate parameters,
including when the optimizer uses a positive numerical floor.}
\item{multi_basin_guard}{Refit the global model if a profile point improves its score.}
\item{n_grid}{Initial grid size per side; endpoints use adaptive root refinement.}
\item{tol}{Maximum accepted absolute endpoint deviance residual.}
\item{num_threads}{Native inner threads. R outer tasks are processed serially.}
\item{seedval}{Main seed; new tasks use domain-separated SHA-256 identities.
Decimal strings preserve unsigned 64-bit seeds exactly.}
\item{simulation_conditioning}{Explicit generating law. Legacy contract C
requires \code{"tip_survival"} or \code{"unconditional_complete"}; A/B must
match their likelihood contract.}
\item{checkpoint_file, resume}{Complete RDS checkpoint plus JSON manifest.
Use \code{resume=TRUE} only with a matching existing checkpoint.}
\item{maxit}{Per-search iteration budget for deterministic bounded optimization;
failed replicate fits receive one larger-budget retry with null-start injection.}
\item{param_x, param_y}{Two distinct free surface parameters.}
\item{grid_size, grid_x, grid_y}{Surface grid size, or explicit in-box vectors.}
\item{reoptimize_nuisance}{Whether remaining free surface parameters are refitted.}
\item{target_conditioning, refit_observed}{Changing contracts requires explicit
observed-data refitting before computing the observed statistic.}
\item{root_distribution}{Model prior or explicitly requested empirical posterior
simulation variant. The latter is not a parametric-null calibration guarantee.}
\item{calc_surrogate_power, alpha}{Optionally simulate the fitted alternative
against the single fixed null critical threshold at level \code{alpha}.}
}
\details{
Contract C always reports \code{profile_support_interval}; A/B use
\code{asymptotic_confidence_interval} or \code{bootstrap_plugin_interval}.
The latter calibrates at one fitted point and is not pointwise test inversion.
Neither mode guarantees finite-sample coverage at boundaries or under weak identification.
Intervals describe the connected component containing the maximum. A feasible
included zero is marked \code{at_physical_zero}; a structurally impossible zero
has score negative infinity and a positive endpoint marked
\code{infeasible_zero_approached}. Search-box limits without threshold crossing
are NA with a censoring status. Endpoint residuals and detected disconnected
support are retained. Improved maxima trigger explicit rebasing.

The bootstrap alternative is started from its own estimate and the fitted null.
Only historical numerical lower floors are widened to include a nested physical
zero; a custom positive exclusion is an error. Failed replicates remain in the
planned denominator. With r exceedances, f failures and B planned replicates,
failure bounds are [(r+1)/(B+1),(r+f+1)/(B+1)]. The Monte Carlo interval is the
Clopper--Pearson envelope [CP lower(r,B), CP upper(r+f,B)], a separate quantity.
The fixed-threshold surrogate is not nested-bootstrap test power.

Raw null/alternative scores and their signed difference are retained separately
for audit. Complete attempt records preserve initial parameters, both alternative
starts, convergence status, original scores, fit parameters and failure reasons.
Failed replicates still have NA inferential scores and remain in the failure
envelope; a negative failed raw statistic is never clipped into a usable replicate.

Checkpoints include input bytes, root weights, bounds, native-library hash,
task identities, simulated observations, fitted parameters, and complete fit
artifacts including ancestral arrays. Changed manifests are rejected.
}
\value{A \code{rustbears_profile_ci} list with intervals, curves, Hessian diagnostics
and rebase metadata; a \code{rustbears_profile_surface} data frame; or a
\code{rustbears_bootstrap_lrt} list with separate failure/Monte Carlo intervals,
all replicate statuses and complete replicate artifacts.}
\seealso{\code{simulate_datasets}, \code{define_root_prior}}
