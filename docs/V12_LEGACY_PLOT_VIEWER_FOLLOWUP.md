# Legacy automatic PDF viewers: isolated release follow-up

This branch starts at scientific candidate
`c7afddba9dd6004a24bfbc99270996c6b3ef4906`. It is a pending Gate 4 release
compatibility change, separate from that candidate's scientific execution,
budgets and CI. It has not been merged into the scientific candidate.

Only five automatic viewer sites change:

| Source | Function | Batch behavior after PDF creation |
|---|---|---|
| `BioGeoBEARS_SSEsim_makePlots_v1.R` | `plot_simDEC_DECJ_inferences` | Honors `openPDFs` only in an interactive R session; subsequent calculations and the existing result are unchanged. |
| `BioGeoBEARS_SSEsim_makePlots_v1.R` | `plot_inf_accuracy_vs_SSEsims_v2` | Closes the PDF and skips the viewer; returns invisible `0L` for compatibility with the previous successful viewer path. |
| `BioGeoBEARS_SSEsim_makePlots_v1.R` | `plot_inf_accuracy_vs_SSEsims_v3` | Same as the preceding row. |
| `BioGeoBEARS_SSEsim_makePlots_v1.R` | `plot_SSEsims_treestats_pt1` | Same as the preceding row. |
| `BioGeoBEARS_stochastic_mapping_v1.R` | `plot_BSM_fns` | Closes the PDF, skips the viewer and retains `return(pdffn)`. |

The three invisible zero values mean that the optional viewer was skipped
after normal PDF closure. They do **not** assert that a viewer opened. Their
`doPDF=FALSE` paths still return `NULL`. Interactive sessions retain the prior
viewer call and its exit status. Function arguments, PDF device closure,
drawing calls and output paths are unchanged. Explicitly invoked `opd` and
`openwd` are outside this change, as are the two viewer guards already in c7.
No likelihood, optimization, simulation, random-number, statistical or
checkpoint code changes.

## Verification and limits

Both edited R files parse with the existing R runtime. Comparing parsed
definitions against c7 confirms that all 21 top-level function signatures are
unchanged and that exactly the five function bodies above differ. The diff
contains five one-line guard changes; `git diff --check` passes.

The existing c7 regression for `hist_event_counts` and `check_ML_vs_BSM`
demonstrated real noninteractive PDF generation with a native `open` sentinel
and a positive control. Its completed log records eight passing expectations.
That evidence supports the guard strategy for those **two** functions only.
This follow-up did not execute the five functions above, select or regenerate
BSM inputs, or start scientific calculations. It adds no test that merely
repeats the guard implementation. Complete plotting-path behavior for these
five legacy helpers remains untested in this follow-up.

The local audit record, parser transcript and hash of the earlier regression
log are retained outside the source repository under
`evidence/legacy_viewer_followup_c7afddb/` in the development stage. Independent
review and any final release integration/CI remain pending. Neither this
branch nor the earlier two-function regression claims universal batch safety
for all retained legacy plotting functions.
