# Figure 1c reporting amendment after source freeze

This isolated reporting branch starts at `b0d6e216658edf9282ab56d8eee05540c32621be`.
It changes the workflow cost reader, plotting code, conversion tests and this
document. It does not change the scientific workflow, timer, production library,
native binary, task identities, locked statistical budget or dedicated cost
budget. The main computation checkout and server snapshot remain unchanged.
Use this separate postprocessing checkout only when genuine formal results are
available; record its commit and script hashes alongside the original execution
source hashes. Do not rewrite a computation manifest to identify this reporting
commit as the executed source.

## API return and interval diagnostics

The frozen Python profile implementation raises `ProfileOptimizationError` when
an endpoint fails numerical certification or a plug-in calibration fit fails.
Such calls remain raised/failed attempts and cannot become fresh complete cost
points. A returned profile may instead have a support interval censored by its
parameter box, disconnected support, or an unavailable local Hessian standard
error. Those are different diagnostics, not numerical failures.

For each completely exported returned profile attempt, the reader now checks the
already bound hashes of `intervals.csv` and `profile_result.json` and carries
every d/e/j endpoint value and status, one-sided reference endpoint, critical
value and interval kind, endpoint residual, disconnected-support
flag, local Hessian standard error (including missing values), and MLE-rebase
flag into the source data. An unknown/uncertified returned endpoint status is
rejected for review, not silently mapped to a successful endpoint. Partial
checkpoint resumes can retain their interval diagnostics but remain outside
fresh complete cost points. An interrupted or incomplete export does not acquire
missing interval diagnostics by imputation.

`workflow_complete_costs.csv` retains API `status` separately from
`process_status`, numerical timings, profile diagnostic fields and artifact
hashes. `workflow_profile_diagnostics.csv` provides one row per returned
parameter/attempt, with task/attempt identities and its fresh-cost eligibility.
Empty endpoint or SE fields remain missing, not zero. Bootstrap/tree-failure
crosses continue to mean retained failed bootstrap replicates or incomplete
trees; profile box censoring and support/SE diagnostics are not labelled as
those failures.

## Planned legend levels and denominators

State colors and reserved-CPU symbols now follow the locked `settings.states`
and `settings.workers` order, independently of which complete observations
survive. Every planned level remains in the legend. A level with no complete
observation is explicitly labelled `0/N`, with N counting its unique planned
tasks across the six workflows and the other grid dimensions. A retry does not
increase the number of observed tasks.

`workflow_cost_task_denominators.jsonl` includes states, workers and repetition
for each task. `workflow_stratum_counts.csv` includes every planned combination
of workflow, states and workers, its planned tasks, unique tasks with a complete
cost observation, total complete attempts and unavailable tasks. Thus an entirely
missing 256-state stratum remains visible as planned even when no point can be
drawn. No missing timing is replaced by zero seconds.

## Validation boundary

The tests cover hash tampering, missing artifact bindings, endpoint-status
semantics, NA-preserving exports, absent state/worker strata, fixed mapping order,
duplicate/missing planned tasks, and repeated-attempt versus unique-task counts.
They are conversion/legend fixtures and do not render quantitative figures.
Formal-input, locked-budget, execution-source and complete-task gates remain
unchanged. No smoke observation is eligible for a manuscript panel.

The three directly relevant test modules passed 39 tests on macOS Python 3.13:
`test_v12_postprocessing_reporting.py`, `test_v12_workflow_costs.py` and
`test_v12_performance_memory.py`. A read-only schema check of the previously
retained small smoke's profile CSVs confirmed that the reader preserves its
box-censored endpoints and absent SEs; this was not a formal-input eligibility
check and generated no figure. A second agent independently reviewed the diff,
producer endpoint semantics, artifact-hash binding and planned denominators
without finding a substantive reporting blocker. This is not a new three-OS
CI certification of the computation snapshot.

Actual formal rendering, the minimum 8 pt typography check, 1.5 pt panel
alignment, PDF text/collision audits and inspection at final physical size remain
pending genuine formal data. Legend additions may require layout adjustment at
that stage; they do not justify shrinking text below the project limit.

The initial source preflight of the plot together with its unchanged shared
figure helper returned 20 PASS, one WARN and zero FAIL. The warning is the
static width detector reading the 180 mm contract value as inches; the actual
figure construction divides the millimetre dimensions by 25.4. This source-only
check does not replace any pending rendered audit.
