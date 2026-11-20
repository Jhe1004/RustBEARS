# RustBEARS Python API & CLI Guide

<p align="left">
  <strong>English</strong> | <a href="Python_and_CLI_Guide_zh-CN.md">简体中文</a> | <a href="../README.md">Back to README</a>
</p>

The `rustbears` Python package (`Python >= 3.10`) provides clean object-oriented workflows, multi-core orchestration, resumable JSONL checkpoints, and a `rustbears` command-line interface backed by the shared Rust engine.

---

## 1. Installation

Ensure `rustc >= 1.82` and `cargo` are installed from <https://rustup.rs/>, then run:

```bash
python3 -m pip install "git+https://github.com/Jhe1004/RustBEARS.git"
```

Or from a local clone:

```bash
git clone https://github.com/Jhe1004/RustBEARS.git
cd RustBEARS
python3 -m pip install .
```

---

## 2. Fitting Models & Survival Conditioning

`rb.BioGeoBEARS` (aliased as `rb.RustBEARS`) accepts Newick file paths or strings for `tree`, and PHYLIP geography file paths or Python dictionaries for `geog`.

```python
import rustbears as rb

tree = "((K_sp1:2.0,O_sp1:2.0):2.0,((O_sp2:1.0,M_sp1:1.0):1.5,(M_sp2:1.0,H_sp1:1.0):1.5):1.5);"
geog = {
    "K_sp1": "1000", "O_sp1": "0100", "O_sp2": "0100",
    "M_sp1": "0010", "M_sp2": "0010", "H_sp1": "0001",
}

# Built-in model names: "DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"
common_opts = dict(
    tree=tree,
    geog=geog,
    max_range_size=4,
    precision="strict_f64",                 # or "legacy_f32" for 32-bit historical emulation
    likelihood_conditioning="tip_survival", # or "legacy_unconditioned" (default) / "unconditional_complete"
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
)

fit_dec  = rb.BioGeoBEARS(model="DEC", **common_opts).fit(num_threads=4)
fit_decj = rb.BioGeoBEARS(model="DEC+J", **common_opts).fit(num_threads=4)

print(rb.compare_models([fit_dec, fit_decj]))
print(fit_dec.summary())
fit_dec.plot("ancestral_dec.pdf", plot_what="pie")
```

### Spatial Modifiers & Time Stratification

Pass file paths (or arrays) for geological time slices and spatial modifiers directly to `rb.BioGeoBEARS`:

```python
strat_fit = rb.BioGeoBEARS(
    tree="Psychotria_5.2_ultrametric.newick",
    geog="Psychotria_geog.data",
    model="DEC+J",
    max_range_size=4,
    timeperiods="timeperiods.txt",
    areas_allowed="areas_allowed.txt",          # monotone island emergence supported in tip_survival!
    distances="distances.txt",                  # optional +x distance-dependent dispersal
    dispersal_multipliers="multipliers.txt",    # optional manual dispersal multipliers
    likelihood_conditioning="tip_survival",
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
).fit(num_threads=4)
```

---

## 3. Multi-Threaded Uniformization Stochastic Mapping

Every `FitResult` (and `ModelAverageResult`) supports deterministic multi-threaded Uniformization stochastic mapping:

```python
# Fixed number of histories (deterministic per map index regardless of num_threads)
maps = fit_dec.run_bsm(num_maps=500, num_threads=4, seed=42)
print(maps.summary_counts["event_means"])
rb.plot_bsm(maps, "bsm_events.pdf")

# Adaptive stopping rule based on target relative Monte Carlo standard error (rMCSE)
maps_adap = fit_dec.run_bsm_adaptive(
    target_rmcse=0.05,
    min_maps=200,
    batch_size=200,
    max_maps=2000,
    num_threads=4,
    seed=42,
)
```

---

## 4. Frequentist Uncertainty Quantification Workflows

### 4.1 Profile Likelihood Intervals & 2D Surfaces

```python
prof = fit_dec.profile_ci(params=["d", "e"], n_grid=16, num_threads=4)
print(prof.to_dataframe())
prof.plot("profile_dec.png")

surf = fit_dec.profile_surface("d", "e", grid_x=12, grid_y=12, num_threads=4)
```

- Inspect `status_lower` and `status_upper` (e.g., `at_physical_zero`, `infeasible_zero_approached`, `censored_at_box_upper`) in `prof.to_dataframe()` to identify boundary optima or flat high-extinction ridges.

### 4.2 Parametric Bootstrap LRT (with Resumable Checkpoints)

```python
lrt = rb.bootstrap_lrt(
    fit_dec,
    fit_decj,
    n_boot=199,
    seed=2026,
    checkpoint_file="bootstrap_dec_vs_decj.jsonl",
    parallel=rb.ParallelConfig(workers=4, native_threads=1),
)
print(lrt.summary())
# If interrupted, re-run with identical arguments and `resume=True`!
```

### 4.3 Akaike Model Averaging & Mapping Precision Diagnostics

```python
avg = rb.model_average([fit_dec, fit_decj], criterion="AIC")  # use "AICc" when N is large enough
print(avg.weights_df)

avg_maps = avg.run_bsm(num_maps=400, num_threads=4, seed=42)
diag = rb.bsm_diagnostics(avg_maps, tie_margin=0.25, tree_target="fixed_sample")
print(diag.event_summary_df)
```

### 4.4 Posterior Tree Propagation (`fit_over_trees`)

 Propagate topological, model, and parameter uncertainty across an external sample of posterior trees (e.g., from BEAST or RevBayes):

```python
multi = rb.fit_over_trees(
    trees="posterior_sample.trees",   # Newick/Nexus path or list of Newick strings
    geog="geography.data",
    models=["DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J"],
    summary_tree="mcc_reference.newick",
    likelihood_conditioning="tip_survival",
    root_prior=rb.RootPriorSpec("uniform_nonnull"),
    criterion="AICc",
    bsm_maps_per_tree=200,
    checkpoint_file="posterior_propagation.jsonl",
    parallel=rb.ParallelConfig(workers=4, native_threads=1),
)
print(multi.clade_summary_df)

# Evaluate posterior-population mapping diagnostics (with Geyer IPS ESS >= 50 check for MCMC chains)
pop_diag = rb.bsm_diagnostics(
    multi.bsm_result,
    tree_target="posterior_population",
    tree_dependence="mcmc",  # or "independent" for Student-t intervals
    tie_margin=0.25,
)
```

---

## 5. Command-Line Interface (`rustbears`)

Run multi-model maximum-likelihood fitting, ancestral state export, and stochastic mapping directly from the terminal without writing Python or R code:

```bash
rustbears fit \
  --tree tree.newick \
  --geog geography.data \
  --models DEC,DEC+J,DIVALIKE,DIVALIKE+J,BAYAREALIKE,BAYAREALIKE+J \
  --max-range-size 4 \
  --likelihood-conditioning tip_survival \
  --precision strict_f64 \
  --maps 200 \
  --threads 4 \
  --seed 12345 \
  --out-dir results_out
```

### CLI Options (`rustbears fit`)

| Flag | Description |
| :--- | :--- |
| `--tree PATH` | Input Newick phylogeny file (required) |
| `--geog PATH` | Input PHYLIP-style discrete geography file (required) |
| `--models LIST` | Comma-separated models, e.g. `DEC,DEC+J,DIVALIKE,BAYAREALIKE` (default: `DEC`) |
| `--max-range-size INT` | Optional maximum number of areas occupied simultaneously |
| `--likelihood-conditioning` | `legacy_unconditioned` (default), `tip_survival`, or `unconditional_complete` |
| `--precision` | `strict_f64` (default 64-bit) or `legacy_f32` (32-bit historical emulation) |
| `--use-ambiguities` | Expand `?` wildcard tip range ambiguities |
| `--timeperiods PATH` | Optional stratum boundary ages file |
| `--areas-allowed PATH` | Optional per-stratum allowed-areas matrix |
| `--areas-adjacency PATH` | Optional per-stratum area adjacency matrix |
| `--distances PATH` | Optional inter-area geographic distance matrix (`+x` free exponent) |
| `--envdistances PATH` | Optional environmental distance matrix (`+n` free exponent) |
| `--dispersal-multipliers PATH` | Optional manual dispersal multiplier matrix (`+w` exponent) |
| `--maps INT` | Number of Uniformization stochastic maps per model (`0` skips BSM) |
| `--threads INT` | Number of CPU threads (default: `1`) |
| `--plot {pie,text}` | Ancestral reconstruction PDF style (default: `pie`) |
| `--out-dir PATH` | Output directory for `model-comparison.csv`, `ancestral.csv/.pdf`, and BSM tables |
