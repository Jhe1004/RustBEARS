# Core-scope reporting over separate immutable budgets

This is a postprocessing adapter. It does not run, resume, fit, simulate, or
alter a scientific job. The original `inputs.v1` readers retain their full-design
requirements. A smaller collection is accepted only through the explicit v2
schemas below, with `formal: true` and `scope: "core_v1"`.

The core scientific design declares eleven fresh jobs: four 40-tip/four-area V1
parameter scenarios spanning interior and boundary values, two V2b scenarios
(`null_e02`, `alt_j003`), both V2c scenarios,
V3 datasets 1–10, and D1/D2. Complete V4 is the twelfth collection member and
retains its original budget, source, seed, records, and task identity. The new
design is a disclosed amendment after the original full execution began. The
reader does not call the unfinished original full suite complete.

The numerical summaries remain unchanged: Clopper–Pearson intervals, planned
failure bounds, dataset-paired V3 scores, and unavailable values are preserved.
Main-panel row identities are checked against the declared Cartesian products,
yielding 16/6/3/4 rows. A correct row count with duplicated or substituted
identities is rejected. The d/e profiles and all declared calibration outcomes
are exported as supplements. D1/D2 are verified as part of the collection and
remain inputs to Figure 3; they are not flattened into Figure 2 calibration
tables. No missing estimate is replaced with zero.

## Figure 2 collection v2

Each reference is `{ "path": "...", "sha256": "..." }`. Paths may be absolute
or relative to the inventory file. Design paths must resolve to the same files
referenced by their own budgets. Run outputs and their full source checkouts
must be locally readable on the rendering host. Copying a summary alone is
insufficient. Preserve the producer's absolute paths when its raw evidence
requires them; the adapter does not rewrite historical manifests.

```json
{
  "schema": "rustbears.figure2.inputs.v2",
  "formal": true,
  "scope": "core_v1",
  "core_design": {"path": "statistical_design.json", "sha256": "ACTUAL_SHA"},
  "panel_selection": "OBJECT_RETURNED_BY_core_panels()",
  "runs": [
    {
      "id": "01_v1_n40_a4_interior",
      "analysis": "v1",
      "origin": "fresh",
      "path": "ACTUAL_CORE_OUTPUT/01_v1_n40_a4_interior",
      "budget": {"path": "CORE_STATISTICAL_BUDGET", "sha256": "ACTUAL_SHA"},
      "design": {"path": "statistical_design.json", "sha256": "ACTUAL_SHA"},
      "source": {
        "root": "ACTUAL_EXECUTION_CHECKOUT",
        "commit": "ACTUAL_40_CHARACTER_COMMIT",
        "sha256": "EXACT_source_sha256_OBJECT_FROM_THIS_BUDGET"
      },
      "expected_task_ids": {"kind": "outer_replicate", "first": 1, "last": 200, "count": 200}
    }
  ]
}
```

The strings prefixed `ACTUAL_` and the two object placeholders are explanatory,
not accepted values. Include every `job_registry` entry exactly once and then
the single `inherited_evidence` V4 entry. Use `verify_registry(design)` from
`scripts/v12_core_reporting.py` to obtain their jobs, origins and exact expected
task IDs. Use `core_panels()` for the panel object. Fresh entries share the new
statistical budget; the V4 entry uses its old budget and old source checkout.
The V4 `expected_task_ids` is a mapping from each of the six sections to the full
list of planned replicate IDs, as returned by `planned_ids(old_settings)`.

V3 also requires `coordinator_evidence` containing hashed `plan` and `ledger`
references. The original V3 producer has no separate first-start timestamp.
The adapter checks its exact processed coordinator task, output path, settings,
budget, plan digest, and recorded start time. All ten dataset ledgers are still
read through the existing V3 verifier, including retained failures. Successful
datasets retain the existing source-to-raw-chain and score checks.

An immutable scientific design can have `source_commit: null` because it was
fixed before source freeze. The final locked budget must provide a nonempty
40-character source identity and complete `source_sha256`; every referenced
source file is rehashed from the separate execution checkout. The design
timestamp must precede lock, and fresh execution must follow specification and
lock. A non-null design source identity must match the budget. Every run is
read with its **own** budget SHA; old V4 cannot borrow the core budget. Gate 2
and Gate 3 must be true in the actual hash-bound gate file for each budget.
That gate must name the same execution source commit and native-library SHA.
Any source-file hashes recorded by the gate must also match the budget.

From the isolated reporting checkout, the non-rendering check is:

```sh
python scripts/plot_v12_calibration.py \
  --inputs /absolute/path/core_figure2_inputs.json \
  --output /absolute/path/new_core_figure2_check
```

Use `--export-tables` to export complete source tables after eligibility checks.
Only after real formal results exist, add `--render --qa-tools-dir
/absolute/path/nature-figure/scripts`. Engineering fixtures and smoke outputs
must never be supplied as manuscript data. The full v1 default template is
unchanged; it is intentionally not a smaller core template.

## Dedicated workflow cost inventory v2

Use schema `rustbears.workflow_benchmark.inputs.v2`, `formal: true`,
`scope: "core_v1"`, and these required fields:

| Field | Required content |
|---|---|
| `budget` | Hash-bound separate core workflow cost budget |
| `workflow_design` | Exact scientific design referenced by that budget |
| `core_design` | Hash-bound eleven-job core statistical design |
| `source` | Actual cost execution root, commit and exact budget source inventory |
| `benchmark` | The complete dedicated benchmark output directory |
| `expected_task_ids` | Ordered `cell_id:phase` list for all 36 planned tasks |
| `workload_labels` | Exact labels derived from the registered payload below |
| `inherited_v4` | The same full original-budget V4 declaration used above |

The six cell IDs are `s16_r1_w1`, `s16_r1_w8`, `s64_r1_w1`, `s64_r1_w8`,
`s256_r1_w1`, `s256_r1_w8`. Within each cell, task phases are
`six_model_average`, `profile_chi2`, `profile_plugin`, `bootstrap_lrt`,
`run_bsm`, `fit_over_trees`, in that order. These are task IDs, not independent
attempt IDs. The labels, in display order, are:

```json
[
  "χ² profiles: d, e, j",
  "Plug-in profiles: B=199",
  "LRT: B=199",
  "Six models + average",
  "5000 single-tree maps",
  "20 trees × 100 maps"
]
```

The cost budget must bind its exact core statistical predecessor budget through
`predecessor_budget` and `predecessor_statistical_budget_sha256`. Both budgets
use the same gate, source commit and source inventory. The predecessor must
contain all eleven prescribed job IDs/settings, a terminal and idle ledger,
and unchanged source evidence. `failed_retained` remains a processed terminal
failure; it is not relabelled success. For processed outcomes, all recorded
result/evidence hashes are rechecked. Inherited V4 is validated independently
under its own budget instead of being invented as a fresh predecessor task.

The actual benchmark plan, full phase grid, input files, attempt directories,
timer records, per-fit events, profile diagnostics and exported tables retain
the existing dedicated reader checks. Only a fresh, returned, complete phase
can supply a full-cost point. Partial resume, timeout, unobserved process exit,
raised call and blocked phases remain unavailable. API return and numerical
diagnostics remain distinct. No timeout or absent observation becomes zero
seconds. All locked states and workers remain in the legend with observed /
planned counts, including 0/N for wholly unavailable levels.

There is one planned repetition per state/worker/workflow. Attempts after
interruption are not independent repetitions. Exports explicitly report
`independent_repetitions: 1` and `between_repeat_variability_estimated: false`;
the caption carries this interpretation. Do not report a between-repeat
interval, stable-scaling claim, or the combined phase times as suite wall time.

```sh
python scripts/plot_v12_performance.py \
  --workflow-benchmark-input /absolute/path/core_workflow_inputs.json \
  --output /absolute/path/new_core_cost_check
```

Without `--render`, this verifies and exports cost source tables only. A final
Figure 1 additionally needs the actual two-host matched scalar benchmark
inputs and the original QA arguments; the 36 cost tasks do not replace them.
No quantitative figure was rendered while this adapter was developed.

## Isolation and QA

The reporting commit changes readers, plots, engineering tests and this
document only. No production Python, Rust or R inference code, scientific
runner, old checkpoint, lock or budget is changed. The execution commit and
postprocessing commit are separate provenance fields. Figure 2 exports record
both plotting and collection-reader SHA; Figure 1 records the collection,
workflow reader, plotting and graphics-helper SHA.

Engineering tests exercise wrong-source and budget links, missing/duplicate
IDs, failed gates, late designs, old V4 relabelling, retained failed jobs,
single-repetition 0/N legends and NA/Clopper–Pearson/failure bounds. They do not
establish scientific calibration. Final real-data rendering still requires
font glyphs at least 8 pt, alignment/collision QA at 1.5 pt tolerance, visual
inspection, source tables and PDF/SVG/PNG plus 600-dpi LZW TIFF exports. Those
checks remain pending until real formal figure inputs are complete.
