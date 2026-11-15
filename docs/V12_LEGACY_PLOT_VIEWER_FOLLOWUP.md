# Legacy automatic PDF viewers: isolated release follow-up

This branch starts at scientific candidate `c7afddba9dd6004a24bfbc99270996c6b3ef4906`. It is a pending Gate 4 compatibility change, separate from the frozen scientific A3 candidate, budgets and validation. It is not merged into the scientific candidate.

## Changes

Five automatic viewer sites call `system()` only in an interactive R session:

| Function | Noninteractive behavior after PDF creation |
|---|---|
| `plot_simDEC_DECJ_inferences` | Skips the optional viewer and proceeds to the existing summary. |
| `plot_inf_accuracy_vs_SSEsims_v2` | Closes the PDF; returns invisible `0L`. |
| `plot_inf_accuracy_vs_SSEsims_v3` | Closes the PDF; returns invisible `0L`. |
| `plot_SSEsims_treestats_pt1` | Closes the PDF; returns invisible `0L`. |
| `plot_BSM_fns` | Closes the PDF and returns `pdffn`. |

The three zero values describe successful PDF closure and skipped optional viewing, not successful GUI launch. Their `doPDF=FALSE` paths still return `NULL`. Interactive calls retain the viewer invocation and its result. Explicitly requested `opd`/`openwd` are outside scope.

Actual execution exposed two additional plotting compatibility defects, fixed on this branch:

- `plot_simDEC_DECJ_inferences` now loads its existing `sim_params_Rdata_table` argument instead of ignoring it and reading `/simdata/BGB/dej_params.Rdata`. The default argument is unchanged. No `/simdata` directory is created.
- `plot_BSM_fns` passes its existing `include_null_range` argument to both color-generation calls (`plot_null_range`) and all four state-plot calls (`include_null_range`). This keeps colors and state labels consistent and preserves the existing default `TRUE`. It does not force empty-state inclusion when the caller selects `FALSE`.

No likelihood, optimization, simulation, random-number, statistical or checkpoint implementation changes. All 21 top-level function signatures in the two R files remain unchanged; the same five function bodies differ from c7.

## Bounded execution evidence

`scripts/validate_legacy_plot_viewers.R` accepts four paths: this worktree, a directory containing the saved three-tip DEC/DEC+J fits and BSM tables, a fresh output directory, and an absolute Rscript executable. It loads the five function definitions into a private environment above the installed RustBEARS namespace. It does not unlock or alter that namespace, and checks that the function bodies and formals match the source.

The actual run used the existing macOS A3 R 4.4.3 environment, saved three-tip fitted objects and one saved stratified map. It generated no fit or map. Three summary helpers used explicitly labeled deterministic six-row engineering inputs, not simulated study observations.

Each function was executed through a controlled interactive positive case and a real noninteractive case. The positive case changes only `interactive()` lookup in the private calling environment; the real `system()` boundary runs an absolute Rscript sentinel which writes a marker and exits 99. The actual plotting function must reach that boundary exactly once. The negative case restores `base::interactive`, must generate a new valid PDF, close the graphics device, and make zero viewer calls. No GUI is launched. The final attempt retains both positive and negative PDFs in separate fresh directories. The preceding five-case attempt intentionally removed its positive PDFs; those deleted files are not counted as currently verifiable artifacts.

The successful engineering attempt completed in 4.15 seconds with these retained noninteractive PDFs:

| Function | Pages | Bytes |
|---|---:|---:|
| `plot_simDEC_DECJ_inferences` | 4 | 14291 |
| `plot_inf_accuracy_vs_SSEsims_v2` | 1 | 7155 |
| `plot_inf_accuracy_vs_SSEsims_v3` | 1 | 6595 |
| `plot_SSEsims_treestats_pt1` | 2 | 9906 |
| `plot_BSM_fns`, empty state included | 1 | 5330 |
| `plot_BSM_fns`, empty state excluded (engineering display) | 1 | 5325 |

PDF headers, EOF markers, page counts and SHA-256 hashes were verified. Color-count checks cover eight states with the empty state and seven without it. All six source data files have identical before/after hashes. The saved SSE summary return object is also retained. These are engineering checks, not scientific calibration or visual publication QA.

## Remaining limitations

The original saved GenSA result is **not** compatible with the old SSE summary reader: the four-page plot succeeds, then its summary reads GenSA's positive negative-log-likelihood as a log-likelihood and rejects negative AIC. That failed attempt and its complete PDF remain retained. The passing case explicitly adapts only an in-memory display copy to the older optim container (`use_optimx=FALSE`, `optim_result$value` set to its already stored `total_loglikelihood`). No score is recalculated, no posterior or other probability array is changed, and the original file is untouched. This does not fix or establish success of the default GenSA summary path; that remains a separate Gate 4 defect.

The retained APE coordinate helper also requires five unqualified APE internals. The engineering child process supplied `node_height`, `node_depth`, `node_depth_edgelength`, `node_height_clado` and `.PlotPhyloEnv` from the installed APE namespace. This documents a required calling context, not a product fix. The SSE helpers require four additional explicit summary bindings: `uniq_params_num`, `siminf_stats_subset`, `ancstate_accuracies_DEC_subset`, and `ancstate_accuracies_DECJ_subset`.

One path for each of the five functions and an additional complete `include_null_range=FALSE` BSM display path were exercised. The latter uses an explicitly labeled engineering copy: it drops an asserted all-zero empty-state probability column and renumbers the same saved nonempty state IDs from 2–8 to 1–7. Event times, labels and nonempty probabilities are preserved. It is not a new fitted model or simulated history. The second-model branch of `plot_BSM_fns` was not executed. This work does not claim arbitrary-input support or universal batch safety for all retained legacy plotting helpers.

The raw attempts, exact commands, source snapshots, original data hashes, scripts, results and PDF inventory are retained outside the repository at `evidence/legacy_viewer_execution_a3aec9f/` in the development stage. Earlier parse-only evidence remains at `evidence/legacy_viewer_followup_c7afddb/`. Main-source validation and final release integration are separate.
