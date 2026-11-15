# Source-bound manuscript figures

The display rules are fixed before formal outcome tables exist. `--output DIR`
alone writes a figure contract; it does not invent data or render placeholders.
The Python backend exports editable PDF/SVG and a 600 dpi PNG, and retains exact
input files, SHA-256 hashes and complete source tables. Eight-point text is the
minimum at the declared final size. Rendered outputs still require PDF text,
collision and final-size visual checks; script completion is not visual approval.

## Figure 1: workflow and distinct computational workloads

`scripts/plot_v12_performance.py --workflow-only --output DIR --audit-scripts QA_DIR`
draws only the non-data methods panel. For the full figure use `--render`, two
`--benchmark HOST=RESULTS_JSON` arguments, and two corresponding
`--p0 HOST=RESULTS_JSON` arguments. The same two host names and benchmark grids
are required. Formal benchmark records must contain every planned cell, including
failures; an unsupported or failed comparison cannot acquire a time ratio.

Panel b uses all matched external/Rust fresh-process timing ratios. Startup,
imports and input parsing are included, so these are not kernel or optimizer
speedups. Each comparator has its own mathematically aligned Rust convention.
Lagrange-ng contributes only the two-area `v=2` stratum. Peak process RSS and
unsupported/failed comparisons remain in the full source table. Panel c shows
every single-worker, conditioned P0 fit, explicitly as a four-tip pilot workload;
all other P0 groups remain in the source table and supplementary figure. This
panel cannot establish full-workflow performance for larger trees.

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
