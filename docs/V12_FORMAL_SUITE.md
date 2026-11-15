# Formal study execution and recovery

`scripts/run_v12_formal_suite.py` constructs all registered V1, V2b, V2c, V3,
V4, D1 and D2 commands directly from the prospective budget. Without `--execute`
it writes a plan and starts no scientific task. The original complete design has
28 workflow jobs; the number of fits inside them is much larger. A planned job
is not evidence that an experiment ran.

The user selected Tongfang for formal execution. Keep immutable source, runtime,
inputs, temporary files, logs, results and archives in numbered subdirectories
under the designated project root. Use a reviewed, committed source snapshot;
do not edit its code while the study runs. Keep the frozen budget and gate
evidence at their recorded paths, because checkpoints bind them by hash and
retain their provenance.

## Before execution

Freeze the statistical `budget_design.json` first, then record Gate 2/3 evidence,
then lock `formal_budget.locked.json` with both hashes. This ordering avoids a
circular dependency between the budget and gate evidence. The final budget
must include exact hashes for all workflow sources and the terminal validators:

- `scripts/run_v12_formal_suite.py`
- `scripts/v12_formal_protocol.py`
- `scripts/plot_v12_calibration.py`
- `scripts/v12_case_evidence.py`

The launcher verifies the locked timestamp, passed gates, source hashes and
explicit resource budget before scientific execution. The V3 adapter verifies
its real external posterior provenance separately. D1/D2 adapters preserve the
same physical parameter bounds across reference-tree and posterior-tree fits.
The suite has no option to relabel a pilot or reduce replicate counts silently.

Example command structure, with actual paths and CPU allocation supplied from
the reviewed deployment:

```sh
python scripts/run_v12_formal_suite.py \
  --budget-json /PROJECT/04_protocol/formal_budget.locked.json \
  --output /PROJECT/08_formal_results \
  --revbayes /PROJECT/02_runtime/revbayes/bin/rb \
  --d2-input-dir /PROJECT/03_inputs/D2_Trimezieae \
  --cpu-budget 32 --workers-per-job 2 --execute
```

`/PROJECT` is a documentation placeholder. It must be replaced with the actual
user-authorized project root. The example CPU allocation is not a measured
speedup or a claim that all cores will be busy.

## Resource and process ownership

The launcher admits jobs in fixed FIFO order while the sum of their reserved
workers stays within `cpu_budget`. `--large-v1-workers` assigns a separate,
prospectively fixed allocation to the four V1 jobs with 150 tips and five areas;
all other jobs use `--workers-per-job`. The Tongfang design reserves at most 40
workers, with eight for each large V1 job and two for each other job. This is a
resource reservation, not CPU affinity or a claim about measured utilization.

Each fitting task uses one native thread and BLAS thread. Single-tree BSM uses
its job's allocated native-worker count without a concurrent outer fitting pool;
these workers are already included in that job's reservation. V4 is serial and
conservatively reserves the regular allocation. The launcher requires normal
process priority. It does not change another process, machine-wide scheduling
policy or existing user job.

A process-held file lock protects both planning and execution. Children inherit
the lock descriptor, so killing the launcher cannot permit another launcher to
overlap still-running child jobs. SIGINT/SIGTERM stop only the launcher's own
process groups, with bounded cleanup. The fixed scheduling order starts empirical
and tree/mapping studies before the larger parameter grids. It does not depend
on a statistical outcome.

## What counts as processed

Exit code zero and a plausible `summary.json` are insufficient. Terminal
verification checks the actual checkpoint configuration and every planned
replicate or dataset ID. It verifies model, contract and parameter-domain
identities, retained failure records, bootstrap counts, tree/model rows, map
identities and event completeness. Evidence hashes are recorded in the suite
ledger and rechecked on recovery.

`processed` means every planned task has a retained outcome. It does not mean
that every fit converged, coverage met a target, or a route was resolved. A
profile marked numerically uncertified and a retained incomplete tree remain
visible. They cannot disappear from the denominator or be retried until they
produce favorable results. An adapter/infrastructure failure is recorded as
`failed_retained`, with its logs; recovery does not automatically rerun it.

## Recovery and artifacts

Restart an interrupted suite using the same command and `--resume`. The plan,
budget, worker allocation, input paths and sources must still match. Completed
evidence is verified, failed terminal tasks remain retained, and interrupted
adapters resume their own checkpoints. Duplicate/unknown ledger IDs, modified
completed results, unowned output directories and a live competing owner are
errors. Investigate infrastructure failures before any separately documented
repair or new run; keep the original failure evidence.

The launcher appends stdout/stderr and atomically updates `suite_ledger.json`.
The final `all_workflows_processed` status additionally requires an outcome for
every registered workflow and unchanged sources. It is not the software's final
scientific release gate.

Keep large NPZ/checkpoint histories and final archives on the execution host.
Render source-bound figures beside those inputs, then copy figures, compact
source tables and hash manifests to the local review machine. An archive digest
does not establish a public deposit, DOI, statistical acceptance or publication.
