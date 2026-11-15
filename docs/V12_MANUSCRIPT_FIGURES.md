# Source-bound manuscript figures

The display rules are fixed before formal outcome tables exist. `--output DIR`
alone writes a figure contract; it does not invent data or render placeholders.
The Python backend exports editable PDF/SVG and a 600 dpi PNG, and retains exact
input files, SHA-256 hashes and complete source tables. Eight-point text is the
minimum at the declared final size. Rendered outputs still require PDF text,
collision and final-size visual checks; script completion is not visual approval.

## Figure 1: workflow and distinct computational workloads

`scripts/plot_v12_performance.py --workflow-only --output DIR --audit-scripts QA_DIR`
draws only the non-data methods panel. For the data figure use `--render`, two
`--benchmark HOST=RESULTS_JSON` arguments, and `--workflow-input INVENTORY.json`.
The benchmark grids must match. Formal benchmark records contain every planned cell, including
failures; an unsupported or failed comparison cannot acquire a time ratio.

Panel b has equal-sized sibling axes: **wall-time ratio (left)** and **peak-RSS
ratio (right)**, both external/Rust on a logarithmic scale, with a reference line at one.
All tips, state counts, hosts and repeats in the matched grid appear as raw
points; no pooled average or interval is inferred. One shared legend sits in a
reserved header inside the bounded b panel. Software colors and host symbols
have the same meaning in both axes; small horizontal offsets retain tip counts
and repeated measurements without changing their state category.

Startup, imports and input parsing are included, so the ratios are not kernel
or optimizer speedups, or measures of kernel-allocation efficiency. Peak RSS is
the whole fresh process measured separately with POSIX `wait4`, including its
interpreter/runtime. Each comparator has its own mathematically aligned Rust
convention. Lagrange-ng contributes only the two-area `v=2` stratum.

Failed or timed-out comparisons have no numerical ratio even when their partial
process measurements exist. They are marked in a separate non-numeric status
strip below each axis. A missing or zero RSS is unavailable, never zero memory;
negative, non-finite or non-numeric measurements are rejected. All raw seconds,
RSS bytes, return codes and timeout flags remain in `benchmark_source.csv`.
`benchmark_counts.json` retains planned and observed denominators, timeouts,
other failures and unavailable RSS pairs; unavailable memory can overlap a
failure. Missingness therefore cannot disappear through successful-only ratios.

Panel c shows
actual formal uncertainty API workflow costs from the locked complete suite,
with fresh full calls separated from partial checkpoint resumes. See
[the workflow cost contract](V12_FORMAL_WORKFLOW_COST_FIGURE.md) for its fixed
six workload classes, source checks and missing-task denominators. Optional
`--p0 HOST=RESULTS_JSON` inputs are supplementary only and can never supply panel c.

This preparation is not a completed Figure 1. Actual calibration/case grids
cover their registered state counts (currently 4/16/32); the original proposal's
64/256-state complete-workflow benchmark requires a separate locked experiment.
The layout remains 180 × 225 mm with an 8 pt minimum. Panel b now includes
memory in the main figure, but its code and preflight checks are not a rendered
layout approval. Real final-size text, collision and alignment QA awaits the
full formal inputs; no formal quantitative panel or mock quantitative figure
has been rendered for this revision.

## Figure 3: locked D2 case

`scripts/plot_v12_case.py --render --case CASE_DIR --tree ORIGINAL_REFERENCE_TREE
--budget-json LOCKED_BUDGET --output DIR --audit-scripts QA_DIR` requires a
completed formal case whose manifest, script hash, actual inputs and every
planned task count match the locked budget. It rejects the three-tree smoke
case. The input is the preselected 57-tip Trimezieae dataset with E/N coding.

All taxa, all three +J model profiles and all three nested comparisons under
both contracts are included. The ancestral comparison includes every internal
clade and nonempty range; missing/nonmonophyletic cases are counted. The event
panel selects d/e/j for their predefined biological roles, while the source
table retains all event types. Legacy profile support is never labelled a
confidence interval. Failure-mass bounds, Monte Carlo intervals and empirical
history quantiles have separate meanings. The deposited tree sample is fixed;
the figure makes no whole-chain convergence claim.

## Reproducible rendering checks

The alignment helper is passed explicitly using `--audit-scripts` instead of
hard-coding a particular user's path. This directory must contain
`audit_panel_alignment.py` from the figure-QA toolchain. Record its version or
hash with the rendered artifact. After rendering, run the toolchain's
`audit_pdf_text.py FIGURE.pdf --min-pt 8`,
`audit_figure_collisions.py FIGURE.pdf`, and source preflight. Retain all reports,
inspect each panel, legend and interval at final size, and re-run the checks
after any layout change. The figure directory must distinguish a render waiting
for QA from an artifact that passed all checks.
