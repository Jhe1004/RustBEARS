\name{probability_contracts}
\alias{define_root_prior}
\alias{prepare_tree_sampler}
\alias{close_prepared_sampler}
\alias{simulate_datasets}
\alias{draw_prepared_bsm}
\title{Explicit probability contracts and prepared native sampling}
\description{Construct a root distribution and sample fixed-tree geographical
states or complete stochastic histories using the native Rust engine.}
\usage{
define_root_prior(mode = c("uniform_nonnull", "custom_weights", "legacy_flat_weights"),
                  custom_weights = NULL)
prepare_tree_sampler(fit_or_model, conditioning = "observed",
                     root_distribution = "model_prior")
close_prepared_sampler(sampler)
simulate_datasets(fit_or_model, n_reps = 100L, conditioning = "from_model",
                  root_distribution = "model_prior", num_threads = 1L,
                  seed = 12345, start = 1, return_node_states = FALSE)
draw_prepared_bsm(sampler, num_maps = 100L, start = 1, seed = 12345,
                  num_threads = 1L, rejection = FALSE)
}
\arguments{
\item{mode}{Nonempty uniform probabilities, explicit custom weights, or historical
nonempty all-one likelihood weights. Simulation probabilities are normalized separately.}
\item{custom_weights}{A nonnegative finite vector in the model state order, optionally
named with canonical zero-based area indices joined by commas; the empty state is
\code{"_"}. The total mass must be positive.}
\item{fit_or_model}{A result from \code{bears_optim_run}, a run object, or (for
\code{simulate_datasets}) an existing prepared forward sampler.}
\item{conditioning}{\code{"observed"} for BSM; \code{"tip_survival"} or
\code{"unconditional_complete"} for forward sampling. \code{"from_model"}
inherits the fitted probability contract and rejects historical scores.}
\item{root_distribution}{\code{"model_prior"}, or the explicitly requested
\code{"empirical_posterior"} simulation variant. The latter draws roots directly
from the fitted posterior and is not a parametric-null calibration guarantee.}
\item{sampler}{A prepared native sampler returned by \code{prepare_tree_sampler}.}
\item{n_reps, num_maps}{Number of data sets or stochastic maps to draw.}
\item{num_threads}{Number of native worker threads, from 1 to 256.}
\item{seed}{A nonnegative exactly represented integer up to 2\eqn{^{53}}-1,
or an unsigned 64-bit integer supplied as a decimal or \code{0x}-prefixed string.}
\item{start}{One-based global replicate or map index. Disjoint batches reproduce
the corresponding rows of a single larger draw.}
\item{return_node_states}{Whether to return simulated internal-node states.}
\item{rejection}{Use rejection bridges instead of uniformization bridges for BSM.}
}
\details{Set \code{run$likelihood_conditioning} to \code{"tip_survival"},
\code{"unconditional_complete"}, or \code{"legacy_unconditioned"}, and optionally
set \code{run$root_prior = define_root_prior(...)} before fitting. The complete
unconditional contract requires \code{run$allow_null_tips = TRUE}. Historical
scores retain their original scale and do not define a unique simulator.

Prepared samplers preserve the historical state order. A sampler owns a native
handle, which is released explicitly or by its finalizer. Closed handles cannot
be reused. Input files changed after fitting are rejected. New probability
contracts currently require the documented validated model support scope.}
\value{\code{define_root_prior} returns a root specification.
\code{prepare_tree_sampler} returns an owned \code{PreparedTreeSampler}.
\code{simulate_datasets} returns a \code{SimulatedBatch} with zero-based state
indices and provenance. \code{draw_prepared_bsm} returns a \code{BSMResult} with
cladogenetic and anagenetic event tables. Closing returns invisible \code{NULL}.}
\seealso{\code{bears_optim_run}, \code{define_BioGeoBEARS_run}, \code{runBSM}}
