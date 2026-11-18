# RustBEARS v1.2 uncertainty workflow in R

This development tutorial exercises the interfaces. Small example runs do not
establish empirical coverage or complete an empirical analysis. Run
`source('examples/uncertainty_workflow.R')` from the repository root after installing
the R package. Its three bootstrap replicates are deliberately a smoke check.

Select `run$likelihood_conditioning <- 'tip_survival'` and
`run$root_prior <- define_root_prior('uniform_nonnull')` before fitting. Contract A
conditions on every tip being nonempty. Contract B uses `'unconditional_complete'`
and requires `run$allow_null_tips <- TRUE`. Contract C, `'legacy_unconditioned'`,
retains the historical score. R now defaults to StrictF64; historical arithmetic
requires `options(RustBEARS.precision='LegacyF32')`.

`profile_ci(dec, params='e')` reoptimizes nuisance parameters and distinguishes
physical zeros, impossible zeros and search-box censoring. Contract C produces
support intervals. Asymptotic and local bootstrap-plugin intervals under A/B
also lack universal small-sample or boundary coverage guarantees.

`bootstrap_lrt(dec, decj, n_boot=399L, checkpoint_file='bootstrap.rds')` refits both
models for every simulated dataset. Failures remain in the planned denominator.
The deterministic `pval_failure_bounds` and finite-sampling `pval_mc_ci95`
Clopper–Pearson envelope are separate quantities. A fixed-threshold surrogate
is not nested-bootstrap test power. Changing the observed scoring contract
requires `refit_observed=TRUE`; C requires explicit simulation conditioning.

`model_average(list(dec,decj), num_maps=500L)` combines a complete set of compatible
fits. Akaike weights are conditional summaries, not Bayesian model probabilities.
AICc explicitly uses tip count as its sample-size convention; choose AIC when
AICc is undefined. Input files changed after fitting must be refitted.

`fit_over_trees(run, trees, models=c('DEC','DEC+J'),
likelihood_conditioning='tip_survival')` accepts Newick/Nexus samples or phylo
objects. Any failed candidate, after retry, isolates its whole tree. The default
aborts above two percent incomplete trees; conservative mode returns labelled
complete-case means and missing-fit probability bounds. Roots have no incoming
branch; undefined branch-bottom probabilities remain NA. Exact clades and MRCA
matching are distinct options, with excess-tip counts reported for the latter.
Unbounded event counts do not receive artificial probability missing-mass bounds.

Checkpoints retain full manifests, probability arrays and complete event histories.
Changed inputs or implementation fingerprints are rejected. The R implementation
currently runs outer tasks serially with a bounded native thread budget.

For prepared mapping, use `runBSM(dec, NULL, bsm_engine='native_whole_tree',
seedval=42L, target_rmcse=.05, min_maps=200L, batch_size=200L, max_maps=5000L)`.
`bsm_diagnostics(maps)` defaults to a fixed supplied tree sample. Population
precision requires explicit independent/MCMC provenance, with scalar-specific
ESS. Route comparisons are near ties only when their entire paired-difference
interval lies within the practical equivalence band, resolved only when wholly
outside on one side, and otherwise unresolved. Rare statistics are not assessed.
These exploratory intervals do not promise simultaneous, selected-rank or
optional-stopping coverage.

New probability contracts and uncertainty workflows support the six standard model
families, spatial/environmental/area modifiers (`+x`, `+n`, `+w`, `+u`,
`areas_adjacency`), and multi-layer time stratification (`timeperiods`) with
chronologically monotone area emergence (`areas_allowed`). Non-monotone
submergence masks, fossil hooks, fixed internal observations and trait-coupled
models remain explicitly rejected outside `legacy_unconditioned`.
