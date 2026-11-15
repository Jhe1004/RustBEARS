# Figure 1c: actual API workflow costs

Status: code and nonvisual validation prepared; formal results and rendered QA
pending. This is the actual calibration/case workload view. It does not complete
the original 16/64/256-state cost benchmark: 64/256-state full workflows require
a separate prospectively locked benchmark. Panel b memory currently remains in
the source table and supplementary evidence, not a main-panel memory curve.

## Fixed display contract

The Python figure uses a 180 × 225 mm layout and at least 8-point text. Panel a
is a methods diagram; b contains mathematically matched fresh-process external
likelihood timings; c displays observed API wall times in six predetermined rows:

| Row | Included measurements, selected before formal results |
| --- | --- |
| χ² profiles | Every V1 scenario, tip/state configuration and outer replicate; d/e/j profiled together |
| Plug-in profiles | Every V1 configuration; d/e/j together, with the locked plug-in B retained in source data |
| LRT 199 + 199 | V2b calls with exactly 199 null and 199 surrogate replicates, every scenario and all three routes |
| Six models + average | V2c's six converged initial fits plus subsequent aggregation, within one invocation and outer replicate |
| 5000 maps | Each D1/D2 contract's single-tree model-averaged mapping call with 5000 planned maps |
| 200 trees × 500 maps | D2's full fit/audit/average/map call, with 200 planned trees and 500 planned maps per complete tree |

No scenario, task or attempt is selected because of its measured cost. Alternative
B values, finite-B runs, every initial fit, case profile and every other recorded
phase remain in complete supplementary source tables. V3 and V4 scientific
outcomes are verified as members of the required suite; they are not currently
instrumented by this API timing recorder and supply no panel-c points.

Colors identify range-state counts; markers identify reserved CPUs. These are
reservations, not CPU affinity or measured utilization. V2c's six initial fits
are serial even when the job reserves two CPUs. Single-tree BSM uses that job's
native worker count; profile, bootstrap and multitree outer workers each use one
native thread. All BLAS thread limits are one. A panel accepts one actual execution
host; moving/resuming across hosts requires an explicit separate display design.

API wall time includes internal checkpoint I/O and excludes later caller-side
table export. The V2c composite is the sum of its seven sequential API calls;
it is not elapsed time between their first and last timestamps. Concurrent
jobs overlap, so no sum of these points is called suite wall time. CPU time and
process-lifetime high-water RSS are retained without being renamed phase peak RSS.

## Complete costs and missing observations

`resume_requested=true` is expected even on a first launcher invocation and is
not evidence of cached work. A complete-cost point requires `status=returned`
and `checkpoint_rows_before=0`. Bootstrap checkpoints must contain every null
and surrogate ID; multitree checkpoints must contain every planned tree. The
scientific evidence validator independently checks those IDs and artifacts.

Every fresh complete attempt remains a point, including repeated complete calls
after an interruption; none is selected as fastest or latest. Repeated attempts
are not independent datasets. Counts show unique tasks with at least one complete
cost divided by planned tasks. A returned call can have retained inner failures;
those points receive a cross and retain their failure count. A 200-tree call with
two failed trees processes 200 planned outcomes but does not produce 100,000 maps;
its actual completed/incomplete counts remain explicit.

Raised calls and partial checkpoint resumes remain in the all-attempt table,
outside complete-cost points. Missing calls are NA, never zero. An abrupt kill
can leave an in-flight call with no timing row; its duration cannot be recovered
from this recorder. A truncated JSONL tail blocks verification pending preservation
and investigation. A job with missing timing rows remains an explicit unavailable
cost task, even when all its statistical outcomes were processed.

## Input and execution

Default preparation creates a contract and inventory template without images:

```sh
python scripts/plot_v12_performance.py --output FIGURE_DIRECTORY
```

The inventory is a local JSON object with schema `rustbears.workflow_costs.inputs.v1`,
`formal: true`, a `budget` object containing its `path` and SHA-256, and a `suite`
directory containing `suite_plan.json` and `suite_ledger.json`. Paths are resolved
relative to the inventory; stored execution paths must remain available for audit.
The default template has `formal: false` and pending paths and cannot authorize
or promote an experiment. Inspect an actual inventory without drawing:

```sh
python scripts/plot_v12_performance.py --output FIGURE_DIRECTORY \
  --workflow-input FORMAL_TIMING_INVENTORY.json
```

Before any quantitative backend opens, `v12_workflow_costs.py` verifies the locked
budget, actual hashed Gate 2/3 passes and any upstream design. It requires the
complete seven-analysis job set, immutable worker arguments, exact unique ledger
IDs and `all_workflows_processed`. A launcher-level retained failure or interruption
does not become a completed suite; retained *statistical* failures inside processed
outcomes remain eligible and visible. Every processed job is rechecked through
`load_run` or `verify_case`, including source/budget settings and per-ID checkpoints.
All frozen execution sources, result evidence and timing-helper hashes are checked.
The reader does not import the POSIX launcher and is portable to Windows.

After those checks and two formal matched benchmark inputs, use `--render` plus
`--audit-scripts` as described in `V12_MANUSCRIPT_FIGURES.md`. P0 inputs are optional
supplementary evidence. They cannot substitute for formal workflow timing. Rendering
writes the complete attempts, complete costs, per-job and per-cost-task denominators,
full source hashes, caption and plotting-script bindings. Keep raw statistical
artifacts on the execution host; the figure's exact raw timing ledgers are also
copied into its source-file archive. Hashes establish provenance, not immunity to
deliberate rewriting of an entire evidence chain.

## Required QA when real data exist

The source preflight and nonvisual tests can run now. After actual rendering,
run alignment at 1.5 pt, minimum PDF glyph size at 8 pt, collision audit and a
panel-by-panel final-size visual review. Repair spacing or enlarge the figure
rather than reducing fonts. Until those checks are complete, output remains
`rendered_awaiting_text_collision_alignment_and_visual_QA`. No mock, smoke or
pilot values are rendered as a manuscript data figure for layout testing.
