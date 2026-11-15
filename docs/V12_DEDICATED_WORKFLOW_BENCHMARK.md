# Dedicated complete-workflow cost protocol

This is a separate computational experiment. It does not change the locked
28-task statistical validation suite or replace its results. The default command
writes a plan without fitting a model:

```sh
python scripts/benchmark_v12_workflows.py --output PLAN_DIRECTORY
```

The prospective default has 18 cells: 16/64/256 range states (4/6/8 areas),
1/8 reserved CPUs, and three repetitions per setting. Every cell uses the same
predefined balanced 19-tip tree and cycling singleton ranges. The six standard
models use exact tip-survival conditioning, uniform nonempty root probabilities,
strict double precision, d/e bounds [0,2] initialized at .015/.008, and free-j
lower bound zero initialized at .03. All remaining model parameters and j upper
bounds are recorded explicitly in the plan. The optimizer limit is 500.

Cells and phases run sequentially. Six initial fits may use up to the cell's
outer-worker allocation, with one native thread per fit and BLAS limits of one.
Profile/bootstrap/multitree workers likewise use one native thread. Single-tree
BSM uses the cell's native worker count. A cell cannot reserve more than 40 CPUs
or more logical CPUs than available; the default maximum is eight. Reservations
are not CPU affinity or evidence of actual utilization.

Each cell measures six complete tasks:

1. All six fits and AICc model averaging, including the initial fits.
2. The complete three-parameter d/e/j χ² profile, with 16 grid points.
3. The complete d/e/j plug-in profile, with B=199 and 16 grid points.
4. DEC versus DEC+J bootstrap LRT with B=199; no surrogate simulations.
5. Model-averaged BSM with 5000 maps.
6. Six-model fitting, audits, averaging and mapping over 20 deterministic copies
   of the tree with branch lengths scaled evenly from 0.9 to 1.1; 100 maps per
   complete tree. These copies are a computational workload, not posterior trees.

Worker comparisons use identical task identities and seeds for the same state
count/repetition; workers are deliberately absent from the RNG key. Different
repetitions retain separate seeds. No inputs, seeds or failures are replaced
according to convergence or speed. Profiles and LRT each start from the same
initial fitted objects; improved profile MLEs and every additional refit are
retained independently, not silently propagated to another cost workload.

## Sources, budgets and predecessor ordering

Formal execution uses its own `rustbears.formal_budget.v1` JSON. Its sole job has
`analysis: "workflow_benchmark"` and every field of `DEFAULT` in the script:
`states`, `workers`, `repeats`, `n_tips`, `n_boot`, `n_boot_plugin`, `n_grid`,
`bsm_maps`, `n_trees`, `maps_per_tree`, `maxiter`, `phase_timeout_seconds`, `seed`.
Any prospective amendment remains explicit through `full_proposal_grid=false`;
a smaller or timed-out experiment cannot be called the completed original grid.

The budget must freeze every file from `source_manifest()` and the actual native
library digest, either as top-level `native_sha256` or its repository-relative
entry in `source_sha256`. It references the same hashed, passed Gate 2/3 evidence
as the statistical budget. In addition it requires:

```json
{
  "predecessor_statistical_budget_sha256": "ACTUAL_LOCKED_STATISTICAL_BUDGET_SHA256",
  "predecessor_budget": {
    "path": "ACTUAL_LOCKED_STATISTICAL_BUDGET.json",
    "sha256": "ACTUAL_LOCKED_STATISTICAL_BUDGET_SHA256"
  }
}
```

These placeholders are documentation, not execution authorization. After the
statistical suite has finished, the executable command is:

```sh
python scripts/benchmark_v12_workflows.py --mode formal \
  --budget-json LOCKED_WORKFLOW_COST_BUDGET.json \
  --after-suite COMPLETED_STATISTICAL_SUITE_DIRECTORY \
  --output NEW_WORKFLOW_COST_DIRECTORY --execute
```

The coordinator acquires both its own lock and the predecessor suite's lock,
verifies that all 28 original job keys/settings match the bound statistical
budget and are terminal with zero reserved CPUs, and retains the predecessor
lock throughout the benchmark. Thus its phases cannot overlap a statistical
suite continuation. Children inherit these lock descriptors: killing the parent
does not permit duplicate execution while a child remains alive. Only owned
process groups are signalled. Execution is POSIX-only; planning and data-reading
helpers remain portable. Nice=0 remains required for real execution.
Formal coordinator and worker processes require the full source checkout and
verify that the imported Python core belongs to it. Plan-only installed-wheel
tests cannot omit core source files from an executed formal experiment.

## Phase clocks, failures and interruption

Each phase runs in a fresh process with a default **14,400-second process limit**.
The process clock includes interpreter startup, cache restoration and final
exports. The API clock excludes those operations and includes its own internal
checkpoint I/O. Whole-process peak RSS is explicitly not API-exclusive memory.
API times from concurrent workers are never summed as total suite wall time.

Initial fits are preserved in JSON plus lossless NPZ arrays. Later phases restore
the exact fitted parameters, likelihood, evaluation arrays and model fingerprint;
they do not refit a supposedly cached MLE or use unsafe pickle loading. Every
actual fit call writes a start event and a terminal event with convergence,
evaluations, parameters, likelihood and elapsed time. A kill may leave a start
without a finish; it remains incomplete. Curves, per-parameter plug-in calibration
replicate tables, bootstrap replicate tables,
all mapping events, seven-field map provenance, tree fits, probability bounds
and diagnostics are retained with checksums.

Timeouts remain `timeout_retained`, failures remain `failed_retained`, and a
failed initial six-model fit leaves dependent tasks explicitly blocked. The
independent multitree task is still attempted. No terminal failure or timeout
is retried automatically. A timeout during CSV export cannot turn a completed
API into a censored API: its terminal timer row preserves the actual API wall
time, while only the process/export remains timed out. Without a terminal API
row, elapsed time from its start to the termination request is reported as a
right-censored elapsed observation; kill-grace time is excluded. If the API
never started, its duration is unavailable. Timeout results supply no completed
inference or full-cost point.

An infrastructure interruption can be continued with the identical command plus
`--resume`; the source, input and full plan fingerprints must match. Existing
checkpoint records make a resumed API cost partial. All attempts, logs and
checkpoints are preserved; no fastest-attempt selection is allowed. New output
directories are required for a changed protocol or source. A small `--mode smoke`
execution exists solely for implementation verification and cannot be promoted.
If a parent dies before recording its child's exit, resume retains that old
attempt as `unobserved_process_exit`. Parent return code, elapsed time, CPU and
RSS remain null. Any complete API timer row survives separately but does not
become a certified fresh cost. Previously observed immutable artifacts are
rechecked before continuation; advancing shared checkpoints are not mistaken
for immutable pre-interruption files.

## Reading and plotting

`v12_workflow_costs.inspect_dedicated_benchmark()` accepts an inventory with schema
`rustbears.workflow_benchmark.inputs.v1`, `formal: true`, a hashed `budget` path
and a `benchmark` output directory. It verifies all 18 × 6 planned tasks (or the
exact prospectively amended grid), source/native bindings, predecessor evidence,
phase logs and primary per-ID artifacts. Complete fresh API measurements and
all censored/failed/partial attempts are exported separately. Missing durations
remain NA. This dedicated reader can supply Figure 1c; the calibration/case
reader remains available as a separate, explicitly labelled workload view.

Formal rendering, text/alignment/collision checks and final-size inspection are
pending genuine formal results. No smoke cost is a manuscript data point.
