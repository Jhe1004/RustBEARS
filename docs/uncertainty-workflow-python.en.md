# Python uncertainty workflow

These small examples exercise the interfaces. A successful run on a small tree is not a validation of coverage, type-I error, or power.

## Fit a declared probability model

```python
import rustbears as rb

tree = "((T1:1,T2:1):1,((T3:.5,T4:.5):.5,(T5:.5,T6:.5):.5):1);"
geog = {"T1": "10", "T2": "11", "T3": "01",
        "T4": "10", "T5": "01", "T6": "11"}
options = dict(tree=tree, geog=geog,
               likelihood_conditioning="tip_survival",
               root_prior=rb.RootPriorSpec("uniform_nonnull"))
null_fit = rb.BioGeoBEARS(model="DEC", **options).fit()
alt_fit = rb.BioGeoBEARS(model="DEC+J", **options).fit()
```

`tip_survival` (A) conditions on every tip having a nonempty range. `unconditional_complete` (B) includes deterministic empty-to-empty transmission at internal nodes and requires `allow_null_tips=True`. `legacy_unconditioned` (C) preserves the historical score and raw root-weight scale. C simulation requires an explicit generating distribution; its score can still be bootstrap-calibrated under that declared distribution.

Custom roots accept exact range labels, for example `RootPriorSpec("custom_weights", {"A": 2, "B": 1, "AB": .5})`, or integer area tuples `{(0,): 2, (1,): 1, (0, 1): .5}`. A/B normalize these weights. C preserves custom weight scale in the score; sampling uses separate normalized probabilities. Input mappings are copied when the specification is created.

New contracts currently support the six standard models on one time layer, distance exponent `+x`, and manual dispersal multiplier exponent `+w`. Stratification, traits, fixed internal observations, root edges, fossil hooks, and ambiguous geography remain explicitly unsupported.

## Profile parameters

```python
profile = alt_fit.profile_ci(params=["d"], n_grid=8, num_threads=2)
print(profile.to_dataframe())
figure = profile.plot("profile.png")
surface = alt_fit.profile_surface("d", "j", grid_x=4, grid_y=4)
```

A/B return `asymptotic_confidence_interval`; C returns `profile_support_interval`. The default reference is the two-sided chi-square-one threshold. A boundary point estimate does not automatically switch the threshold to a mixture distribution. `include_onesided_zero_ref=True` adds a separate one-sided zero-boundary reference.

`at_physical_zero` means the exact zero belongs to the support set. `infeasible_zero_approached` means zero is structurally impossible and the actual lower endpoint is positive. `censored_at_box_upper` with `upper=None` means the search box ended before a crossing; it does not establish mathematical nonidentifiability. Results describe the connected component containing the maximum. Detected re-entry into support farther along the grid sets `disconnected_support_detected=True`.

`method="bootstrap_plugin"` simulates at the global estimate and refits both the unrestricted model and the model constrained to the original parameter value. This is a local plug-in calibration, not pointwise test inversion. A failed calibration raises an error rather than silently discarding failed datasets before taking a quantile.

## Bootstrap model comparisons and resume

```python
test = rb.bootstrap_lrt(null_fit, alt_fit, n_boot=19, num_threads=2,
                        seed=20261115, checkpoint_file="demo-lrt.jsonl")
print(test.summary())
# After interruption, retain the same inputs, planned count, seed and options:
# test = rb.bootstrap_lrt(null_fit, alt_fit, n_boot=19, num_threads=2,
#                        seed=20261115, checkpoint_file="demo-lrt.jsonl", resume=True)
```

Nineteen replicates only demonstrate the workflow. Every replicate refits both models and includes the fitted null as an alternative-model starting point. The new workflow may extend a recognized historical numerical floor to physical zero. A custom strictly positive bound that excludes the null is rejected. Inspect `test.effective_parameter_bounds` for the actual bounds used.

`pval_failure_bounds` bounds the unknown tail positions of failed replicates. `pval_mc_ci95` is the distinct Monte Carlo envelope `[CP_L(r,B), CP_U(r+failed,B)]`. The latter remains nonzero in width when no replicate fails. Failed replicates remain in `replicates_df` and in the planned denominator.

Changing the fitted target requires `target_conditioning=...` and `refit_observed=True`; both observed-data models are actually rebuilt and refitted. To retain C scoring, also specify `simulation_conditioning="tip_survival"` or `"unconditional_complete"`.

Checkpoints verify inputs, root weights, bounds, precision, seeds, software and native-library identity. A different task requires a different file. `calc_surrogate_power=True` returns `surrogate_power_fixed_threshold`, which reuses the original null threshold. It is not the power of the complete nested-bootstrap test.

## Average models and diagnose mapping precision

```python
average = rb.model_average([null_fit, alt_fit], criterion="AICc")
print(average.weights_df)
maps = average.run_bsm(num_maps=200, seed=7)
diagnostics = rb.bsm_diagnostics(maps)
adaptive = alt_fit.run_bsm_adaptive(min_maps=200, batch_size=200,
                                  max_maps=600, seed=8)
```

AIC/AICc weights are frequentist model weights, not posterior model probabilities. Each model still uses a conditional point estimate. `fit_over_trees` summarizes a supplied external tree sample; geography does not update its tree posterior weights. A tree with a missing candidate-model fit is isolated as a whole, with explicit missing-mass bounds. See `help(rb.fit_over_trees)` for its arguments.

Use `ParallelConfig(workers=..., native_threads=1)` for outer parallelism. Multiple outer workers cannot each start multiple native workers. BLAS is limited to one thread for the entire executor lifetime. Task identities and streams are independent of worker scheduling.

## Predeclare validation runs

The validation script writes a manifest without starting simulation unless `--execute` is supplied. Smoke and formal jobs use separate manifests and directories. Formal execution remains subject to the project's resource and cost gates.

```bash
python scripts/run_uncertainty_validation.py --analysis v1 --mode smoke \
  --scenario interior --output validation-smoke-v1 --execute
python scripts/run_uncertainty_validation.py --analysis v2b --mode smoke \
  --scenario null_e02 --output validation-smoke-v2b --execute
python scripts/run_uncertainty_validation.py --analysis v2c --mode smoke \
  --scenario dominant --output validation-smoke-v2c --execute

# Write a formal preregistration only; do not start the formal computation:
python scripts/run_uncertainty_validation.py --analysis v1 --mode formal \
  --scenario double_zero --n-tips 40 --areas 4 --output planned-v1-double-zero
python scripts/run_uncertainty_validation.py --analysis v2b --mode formal \
  --scenario alt_j003 --finite-b-extrapolation --output planned-v2b-alt-j003
```

Formal V1 uses 200 outer datasets per scenario. V2b uses 200 outer datasets with B=199 for its main test; finite-B comparisons use a shared largest sample with B=399. V2c uses 300 datasets per scenario and retains true internal-node states for Brier scores, log scores, and reliability bins. Fixed-threshold surrogate power and complete nested rejection rates are saved separately. Report paired outer-sample uncertainty for finite-B extrapolation and do not treat an extrapolated intercept as an exact infinite-B result. [Method reference: Boos & Zhang (2000)](https://doi.org/10.1080/01621459.2000.10474226).
