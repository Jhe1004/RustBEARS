# V3 / V4 execution boundaries

## Hardware pilot before formal simulations

`python scripts/run_v12_p0.py --output P0_OUTPUT --worker-grid 1,2,4`
runs 20 actual complete DEC+J fits for each of 16, 57 and 128 states, with
one native and BLAS thread per outer worker. It also measures the historical
contract serially on the same observations. Each group uses a fresh process
and records peak RSS, full fit outcomes, seeds, source fingerprints and platform.
`--smoke` uses two fits only and explicitly cannot satisfy the P0 gate.
Run the declared pilot on both macOS and AMD Linux before fixing formal budgets.
These four-tip costs must not be reported as large empirical-tree costs.

## Mapping calibration

`validate_v12_v4_protocol.py --mode plan --output-dir PATH` freezes the complete
V4 settings and exact count moments without starting the repeated experiment.
The formal grid is: fixed 400 × 500 maps; fixed 5 trees × 20 maps for 500 repeats;
independent 10 trees and correlated 60-tree chains × 20 maps for 500 repeats;
500 repeats at each route count 50/200/1000/5000; and 500 runs of the actual
200-minimum/200-batch/5000-maximum adaptive algorithm at target RMCSE 0.05.

Use `--mode smoke` only for plumbing checks. Use `--mode formal` after the
intended server gate. `--sections fixed,independent` selects a declared subset;
a subset never receives the full-protocol flag. Each output directory has an
immutable manifest, JSONL commits, NPZ shards and descriptive results.
`--resume` requires exactly the same task fingerprint. The correlated-tree
experiment uses a known stationary finite-state chain with persistence 0.8,
preserving draw order. This is a controlled V4 population, not a claimed
phylogenetic posterior from sequence data.

Exact first and second count moments use marked CTMC block exponentials and
the tree product rule. The independent reference reuses the established local
cladogenetic weight kernel. Fixed-count and actual stopping-time coverage are
reported separately with binomial intervals. The early-stop metric checks all
predeclared true-active event and route quantities using their exact means and
variances. Route comparisons are paired, nominal and exploratory. No sequence
coverage guarantee or forced monotonic empirical trend is asserted.

The earlier `validate_v12_v4_precision.py` experiment is retained separately.
Its completed 400 × 500 run does not replace the complete V4 protocol. Its
additional five-tree population experiment revealed small-K normal-interval
undercoverage, and the adaptive experiment showed lower empirical coverage.
Neither result is removed or repaired by choosing a different seed.

## Genuine external posterior interface

An executable two-chain RevBayes adapter is provided as
`run_v12_revbayes_posterior.py --dataset DATASET --executable /path/to/rb`.
It uses `v3_jc69_coalescent.Rev.template`: JC69, the known fixed clock rate,
and a constant-population coalescent with an Exponential(rate=0.1) population
prior. The molecular model reads the alignment alone. Its tree prior differs
from the root-rescaled simulation distribution and is recorded explicitly.
Install `arviz==0.22.0` in an isolated validation environment to calculate
rank-normalized split R-hat and bulk/tail ESS. This is not a package dependency.

The adapter runs two separately seeded processes, discards 25% of each recorded
chain after tuning, checks numeric traces and all observed nonconstant clade
indicators, and keeps every raw trace. Constant clades supply no mixing evidence.
These diagnostics cannot prove exploration of unvisited topology modes.
The MCC is the sampled tree maximizing the product of pooled clade frequencies;
its sampled node ages are retained. The pooled posterior interleaves the two
chains after their separate burnin. Import it with `collect --burnin-fraction 0`;
the collector rejects additional burnin for this already-trimmed output.

Model/move declarations follow the [official coalescent tutorial](https://revbayes.github.io/tutorials/coalescent/Constant).
Diagnostic definitions follow [ArviZ R-hat](https://python.arviz.org/en/v0.22.0/api/generated/arviz.rhat.html)
and [ESS](https://python.arviz.org/en/v0.22.0/api/generated/arviz.ess.html).
The four-tip, 80-site integration run exercises this path; it does not replace
the planned 100-dataset, 50-tip V3 calibration.

`validate_v12_v3_external_posterior.py prepare --output-dir PATH` prepares the
formal dimensions (100 datasets, 50 tips, 4 areas; 2000 DNA sites by default).
It generates actual coalescent trees, exact JC69 sequences and geographic
histories under the declared probability contract. Smaller dimensions are
explicit smoke inputs. A completed preparation means **Bayesian inference is
pending**; it does not create a posterior tree collection.

`run --dataset DATASET --engine-config ENGINE.json` requires a real Bayesian
executable and a separately validated engine input template. Its JSON fields
are documented in the script header. The script substitutes explicit input
tokens, invokes argument lists without a shell, records executable/template/
alignment hashes and captures logs. Engine templates and versions must be
validated on the deployment host; none is silently invented by this interface.
The configuration must produce actual posterior trees, an MCC tree and an
external diagnostic report. Multiple independent chains, ESS and R-hat remain
the external workflow's responsibility and are labelled as externally supplied
evidence. The current gate is chains ≥ 2, minimum ESS ≥ 200, maximum R-hat ≤ 1.01.

`collect --dataset DATASET` verifies a successful recorded run and output
hashes, removes the declared burn-in, retains 200 actual posterior draws in
their original order, and converts branch units using the known fixed clock
when necessary. It rejects missing posterior output and insufficient draws.
It never substitutes perturbed true trees, bootstrap trees or replicated truth.

`evaluate --dataset DATASET` performs the MCC and multi-tree geographic fits.
Paired ground-truth scoring is restricted to internal clades present in both
the true tree and MCC, with excluded true/MCC clades reported explicitly. Scores
are separated into the requested external clade-support bins. This comparison
does not pretend to score a nonexistent true ancestor for an incompatible
clade. The current local smoke has tested sequence/geography preparation only;
no external posterior engine, convergence claim or formal V3 result is supplied.
