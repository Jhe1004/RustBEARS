# Changelog

## 1.2.0 (2026-11-21)

- **Survival-conditioned & complete unconditional likelihoods**: Introduced exact whole-tree tip-survival conditioning (`tip_survival`) and complete unconditional (`unconditional_complete`) formulations alongside legacy BioGeoBEARS compatibility (`legacy_unconditioned`), supporting DEC, DIVALIKE, BAYAREALIKE, +J models, spatial/environmental/area modifiers (`+x`, `+n`, `+w`, `+u`, `areas_adjacency`), and multi-stratum time slicing (`timeperiods`) with chronologically monotone area emergence (`areas_allowed`).
- **Default 64-bit precision**: Switched default arithmetic to IEEE-754 double precision (`StrictF64` / `strict_f64`) while preserving optional `LegacyF32` / `legacy_f32` mixed-precision emulation for exact literature auditing.
- **Endpoint-conditioned Uniformization stochastic mapping**: Added deterministic multi-threaded whole-tree Uniformization stochastic mapping alongside legacy branch-wise rejection, with reusable prepared tree samplers across R and Python.
- **Integrated frequentist uncertainty quantification**:
  - Profile likelihood intervals (`profile_ci`) and 2D surfaces (`profile_surface`) with physical-zero probing and upper-ceiling flat-ridge diagnostics.
  - Parametric bootstrap likelihood-ratio tests (`bootstrap_lrt`) with null-optimum warm starts, conservative Monte Carlo p-values, Clopper–Pearson envelopes, failure bounds, and resumable checkpoints.
  - Akaike model averaging (`model_average`) across ancestral probabilities and stochastic mapping quotas.
  - External posterior tree propagation (`fit_over_trees`) with clade monophyly support, clade-conditional ancestral reconstructions, and deterministic missing-fit probability bounds.
  - Stochastic mapping precision diagnostics (`bsm_diagnostics`) with Student-t intervals for independent trees, Geyer IPS MCMC effective-sample-size (`ESS >= 50`) checks, and equivalence-band route tests.
- **CLI & documentation**: Added `--likelihood-conditioning` to the `rustbears fit` CLI and expanded bilingual English and Simplified Chinese user guides for R, Python, and CLI workflows.

## 1.1.4

- Initial unified Rust computational engine (`bgb-core`, `bgb-numerics`, `bgb-clado`, `bgb-Step`, `bgb-bsm`, `bgb-ffi`) with 100% R `BioGeoBEARS` export compatibility (545 public symbols) and native Python interface.
