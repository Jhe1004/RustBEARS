# Measured workflow costs

The V1/V2b/V2c and empirical-case adapters append `workflow_timings.jsonl`.
Each line records one actual public API invocation, including a raised error.
The timer uses a monotonic wall clock; it does not estimate cost from optimizer
iterations or multiply a small pilot into an observed full-workflow result.
The helper's source hash, host, process, invocation, worker allocation, thread
environment, phase and scientific task identifiers accompany each record.

Wall time includes internal checkpoint I/O and excludes the caller's subsequent
CSV/NPZ exports. CPU time is for the process, including its threads. Peak RSS is
the process lifetime high-water mark, not the incremental memory of a phase.
Concurrent calls overlap, so their wall times cannot be summed into suite elapsed
time. A killed process can leave an unfinished call with no timing row; neither
its duration nor a successful completion is imputed.

The suite supplies `--resume` even on its first invocation. Therefore
`resume_requested=true` does not establish that any work was cached. The before
and after checkpoint counts identify completed checkpoint records, not optimizer
calls. A call with existing checkpoint records is a partial/resumed cost and must
not be presented as a fresh complete-workflow measurement. Repeated invocations
remain in the evidence, including their failures.

V1 times all three parameter profiles together, separately for the chi-square
and local plug-in methods. V2b times each complete bootstrap route, including its
registered surrogate simulations. V2c records each of six fits and the subsequent
model-average aggregation separately; only their matching invocation and outer
replicate can define the combined six-model cost. Empirical cases additionally
record single-tree mapping and the entire multiple-tree fit/map call. Initial
fits, profiles and mapping phases have distinct scopes and must not be silently
combined or substituted for one another.

Formal figure input must be bound to the locked budget and verified outcomes.
Smoke and pilot timings remain labelled as such. Failed computations and resumed
partial costs are retained even when they cannot supply a complete-cost point.
Final statistical results may be unavailable despite an API returning normally;
the outcome and artifact validators, rather than timing status alone, decide
what was scientifically processed.
