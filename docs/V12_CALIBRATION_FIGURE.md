# Figure 2 requires completed formal inputs

`scripts/plot_v12_calibration.py` prepares the calibration figure contract,
checks an explicit input inventory, and exports every registered numerical
result. The default opens no graphics backend. Formal calibration remains
pending; smoke, pilot and synthetic test results cannot supply manuscript panels.

## Fixed comparisons

| Panel | Main comparison fixed before results | Unit |
| --- | --- | --- |
| a | j, 4 areas, 40/150 tips, interior/e=0/j=0/e=j=0; profile, zero-boundary upper reference, local Wald, plug-in bootstrap | Simulated dataset |
| b | All five V2b scenarios; A, C, explicit C-to-A observed refit | Outer simulated dataset |
| c | Multiple-tree minus MCC Brier score in three shared-clade support strata | Paired dataset, after averaging its eligible clades |
| d | Fixed-count and adaptive d/e mapping intervals | Independent mapping experiment |

The zero-boundary reference remains unavailable when the fitted interval does
not define it. A point with zero evaluable outcomes displays NA. The program
does not fill missing estimates with zero, the nominal rate or another method.
All d/e/j interval endpoints, finite-B and surrogate results, V2c predictions,
V3 calibration/log scores and missing-result bounds, and V4 tree-population,
route-grid and stopping results remain in the supplementary source data.

## Prepare and check without drawing

```sh
python scripts/plot_v12_calibration.py --output /absolute/path/figure2_preparation
python scripts/plot_v12_calibration.py --inputs /absolute/path/figure2_inputs.json \
  --output /absolute/path/figure2_checked
```

The first command writes `figure2_contract.json`, `input_inventory_template.json`
and `input_audit.json`. Its template intentionally has `formal=false`. A completed
formal inventory has this shape; add every required job, including the secondary
V2c jobs and finite-B V2b job:

```json
{
  "schema": "rustbears.figure2.inputs.v1",
  "formal": true,
  "budget": {"path": "formal_budget.locked.json", "sha256": "actual-file-sha256"},
  "runs": [
    {"analysis": "v1", "path": "V1_40_interior"},
    {"analysis": "v2b", "path": "V2b_null_e0"},
    {"analysis": "v2c", "path": "V2c_higher_jump"},
    {"analysis": "v3", "path": "V3"},
    {"analysis": "v4", "path": "V4"}
  ]
}
```

This abbreviated example cannot pass because it omits required jobs. Paths are
relative to the inventory file, or absolute. The complete locked budget defines
the required set. Its SHA, actual Gate 2/3 passes and optional frozen-design SHA
must agree. Each run must independently record formal mode and the same budget
payload. Merely changing the inventory flag cannot promote a smoke run.

V1/V2 checkpoints must match the entire preregistered canonical manifest and its
recomputed fingerprint. Every planned outcome ID must occur once. The execution
start must follow the lock. V4 applies the same fingerprint check and requires
all registered section/replicate IDs. V3 verifies its prepared study and retained
per-dataset ledgers, including the score-to-import-to-external-run-to-raw-chain
hash chain and diagnostic requirements for each completed dataset. Failed
datasets remain in the ledger and denominator. An unprocessed task blocks export.

The checker compares recorded V1/V2 and V4 workflow-script digests with the frozen
budget. V3's coordinator records production/native digests and the budget payload;
its own script digest comes through that budget, without a separate execution-time
script field. The plotting audit preserves this distinction. It verifies recorded
provenance and file integrity; it does not recreate past processes or rerun fits.

## Export all source tables before drawing

```sh
python scripts/plot_v12_calibration.py --inputs /absolute/path/figure2_inputs.json \
  --output /absolute/path/new_figure2_export --export-tables
```

`source_data/figure2_[a-d].csv` contains every main-panel point, available N,
planned N and missingness. `all_summary_values.csv` preserves every summary leaf;
`all_run_summaries.jsonl` preserves the complete summary objects. Separate CSVs
hold all V1 parameters, all V2b routes, V2c predictions/weights and V3 shared-clade
and dataset scores. `all_outer_and_mapping_records.jsonl` retains all outer
outcomes and V4 numerical mapping records. JSON values remain explicit inside CSV
cells; an empty CSV field represents missingness, not zero.

`checkpoint_index.json` lists all original JSONL checkpoint paths and SHA256s,
including inner profile/bootstrap records. The original array shards remain in
their run directories. The index retains the checkpoints' shard references and
stored digests; plotting does not read or reverify every large array shard.
The run directories must accompany an archive, not only the flattened CSVs.

For rates, available-outcome estimates and their exact Clopper-Pearson (CP)
intervals use the actual number of evaluable independent experiments. A separate
deterministic failure range assigns all unknown outcomes to failure or success.
The main plot's thin bars use the CP envelope over these possible counts with the
planned denominator; thicker segments mark the failure range. V2b also includes
inner-bootstrap fit failures in its rejection bounds. These envelopes do not
establish nominal coverage of the underlying profile or adaptive-map intervals.
V3 uses the coordinator's whole-dataset paired bootstrap interval. Support-bin
membership for failed datasets is unknown; the source tables preserve that limit.

## Rendering and final-size QA

```sh
python scripts/plot_v12_calibration.py --inputs /absolute/path/figure2_inputs.json \
  --output /absolute/path/new_figure2_render --render \
  --qa-tools-dir /absolute/path/to/nature-figure/scripts
```

The renderer checks the source files again before importing matplotlib. The
working layout is 183 by 230 mm, with all text at least 8 pt at that final size.
These are project design settings, not a verified journal production limit.
Matplotlib exports editable PDF/SVG and 600-dpi TIFF/PNG. The QA directory must
provide `validate_figure.py`, `audit_panel_alignment.py`, `audit_pdf_text.py` and
`audit_figure_collisions.py`. Their availability is an explicit runtime dependency;
the repository does not vendor a private skill installation.

The final data-axis rectangles must pass 1.5-pt alignment checks. The final PDF
must pass the 8-pt glyph and collision gates before candidate files move to the
output root. The state then remains `rendered_pending_final_size_visual_review`,
with `publication_ready=false`, until a person inspects all four panels at final
size. Failed candidates stay under `qa_candidate/` for diagnosis. Use a new output
directory for each attempt so earlier artifacts cannot obscure a failed attempt.

At preparation time, unit tests exercise provenance rejection, arithmetic and
lossless source export without generating images. Source preflight alone does
not establish legibility, collision freedom, statistical validity or final QA.
