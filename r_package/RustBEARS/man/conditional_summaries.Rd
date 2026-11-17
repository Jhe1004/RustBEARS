\name{conditional_summaries}
\alias{model_average}
\alias{fit_over_trees}
\alias{bsm_diagnostics}
\alias{plot.rustbears_model_average}
\alias{plot.rustbears_multitree}
\alias{plot.rustbears_bsm_diag}
\title{Conditional model and tree summaries and mapping precision}
\description{Combine compatible likelihood fits, match labelled clades over
external trees, and estimate stochastic mapping uncertainty with an explicit target.}
\usage{
model_average(res_list, criterion = c("AICc", "AIC"), num_maps = 0L,
  num_threads = getOption("RustBEARS.threads", 1L), seedval = 12345L)
fit_over_trees(BioGeoBEARS_run_object, trees, summary_tree = NULL, models = NULL,
  likelihood_conditioning = c("legacy_unconditioned", "tip_survival", "unconditional_complete"),
  root_prior_spec = NULL, match_mode = c("monophyletic", "mrca"),
  init_strategy = c("warm_hybrid", "full_cold_start", "fixed_ref_params"),
  cold_start_audit_fraction = 0.05, model_averaging_within_tree = TRUE,
  max_tree_failure_rate = 0.02,
  on_excess_failure = c("abort", "conservative_bounds"), bsm_maps_per_tree = 0L,
  max_trees = NULL, num_threads = getOption("RustBEARS.threads", 1L),
  seedval = 12345L, checkpoint_file = NULL, criterion = c("AICc", "AIC"))
bsm_diagnostics(bsm_output, res = NULL, target_rmcse = 0.05, tie_margin = 0.25,
  tree_target = c("fixed_sample", "posterior_population"), tree_dependence = NULL,
  top_k_routes = 5L, min_mean_count = 1, alpha = 0.05)
\method{plot}{rustbears_model_average}(x, ...)
\method{plot}{rustbears_multitree}(x, min_clade_support = 0.5, ...)
\method{plot}{rustbears_bsm_diag}(x, ...)
}
\arguments{
\item{res_list}{Complete list of converged fits sharing labelled observations,
tree, state space, root weights, precision and likelihood contract.}
\item{criterion}{AIC or AICc; the number of tips is the stated AICc sample-size convention.}
\item{num_maps, bsm_maps_per_tree}{Total model-stratified maps, or maps per complete tree.
Positive-weight strata need at least two maps. Zero skips mapping.}
\item{num_threads}{Native thread budget. This R implementation processes outer
tasks serially; worker threads are used within a native operation.}
\item{seedval}{Master seed; task identities derive distinct checked substreams.}
\item{BioGeoBEARS_run_object}{Run template with geography, state space and modifiers.}
\item{trees}{A phylo, multiPhylo, list of phylo objects, Newick strings, or a
Newick/Nexus file. Every tree must contain the same uniquely labelled taxa.}
\item{summary_tree}{Reference tree for clade summaries; defaults to the first tree.}
\item{models}{Candidate standard model names or BioGeoBEARS model objects.
Defaults to the model in the run template. Candidate definitions must be unique.}
\item{likelihood_conditioning, root_prior_spec}{Explicit likelihood contract
and root specification; see \code{define_root_prior}.}
\item{match_mode}{Match exact monophyletic clades, or the MRCA with recorded excess tips.}
\item{init_strategy}{Cold fits, fits initialized from the first tree with a
cold-start audit, or evaluation at the first-tree fitted parameters. The last
mode is a parameter-fixed sensitivity analysis, not re-estimation on each tree.}
\item{cold_start_audit_fraction}{Fraction of trees, prioritizing topologically
distant trees, given an additional independent cold fit. Better cold optima warn
and replace that tree's warm fit.}
\item{model_averaging_within_tree}{Average the complete candidate set within each tree,
or use its lowest-criterion model.}
\item{max_tree_failure_rate, on_excess_failure}{After a cold retry, isolate any
tree with an incomplete model set. Abort above the threshold by default. The
error condition retains the result; conservative mode returns labelled bounds.}
\item{max_trees}{Optional positive maximum count, preserving original order.}
\item{checkpoint_file}{Optional full-manifest checkpoint. Changed inputs are rejected.}
\item{bsm_output}{A native BSMResult with complete map identifiers and strata.}
\item{res}{Reserved fitted-result context; not required for native output.}
\item{target_rmcse}{Target relative Monte Carlo standard error among assessed statistics.}
\item{tie_margin}{Practical equivalence margin in events per map, for paired route differences.}
\item{tree_target}{Fixed supplied tree sample by default, or explicitly selected tree population.}
\item{tree_dependence}{For a population target, explicitly choose independent trees
or an ordered MCMC chain. The latter uses a separate initial-positive-sequence ESS
for every scalar statistic and paired route difference.}
\item{top_k_routes}{Number of ranked routes to compare in adjacent pairs.}
\item{min_mean_count}{Minimum mean events before a statistic is assessed for precision.}
\item{alpha}{Per-comparison interval error probability.}
\item{x}{A result of the corresponding class.}
\item{min_clade_support}{Hide low-support clades in the probability heat map.}
\item{...}{Additional base graphics arguments.}
}
\details{These are plug-in conditional ML/Akaike summaries over external trees,
not a joint Bayesian posterior. AICc weights use tip count as an explicit working
convention, which does not establish independent phylogenetic observations.

A failed candidate isolates its entire tree. Complete-case probability means
are accompanied by deterministic lower/upper missing-fit bounds on the original
eligible tree set. Undefined incoming-root branch probabilities remain NA.
Event counts are unbounded, so joint BSM is labelled complete-case and receives
no fictitious bounded-probability missing-mass correction.

Fixed-tree MCSE combines within-stratum variances with squared weights.
Population precision uses between-tree mean variability and statistic-specific
ESS; equal tree weights are currently required. Independent tree-population
intervals use Student-t critical values with K - 1 degrees of freedom; ordered
MCMC tree-chain intervals use Student-t critical values with ESS - 1 degrees of
freedom when effective_n_trees >= 50, and otherwise retain se_target and
effective_n_trees while reporting NA interval endpoints and
precision_target_met = FALSE. The historical 95 percent range of sampled
histories differs from the Monte Carlo interval for their mean.

Route intervals are approximate exploratory per-comparison intervals. They do
not guarantee simultaneous coverage, selected-rank coverage, or coverage after
adaptive stopping. Near-tie requires the whole paired-difference interval inside
the equivalence band; resolved requires it wholly outside that band on one side.
Other cases are unresolved. Rare statistics are explicitly not assessed.
}
\value{A rustbears_model_average, rustbears_multitree or rustbears_bsm_diag list
with probability arrays, full provenance tables, failure bounds or precision
summaries. Plot methods return the input invisibly.}
\seealso{\code{profile_ci}, \code{bootstrap_lrt}, \code{draw_prepared_bsm}}
