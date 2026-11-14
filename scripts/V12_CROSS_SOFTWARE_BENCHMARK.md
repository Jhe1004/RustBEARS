# Cross-software alignment before performance claims

The executable protocol is `benchmark_v12_alignment.py`. Its two concrete
adapters are RustBEARS and the preserved original BioGeoBEARS 1.1.4 runtime.
The default invocation performs only three small fixed-parameter scalar cases.
It records public-API elapsed time separately from process startup wall time,
whole-process peak RSS, commands, executable/native hashes and stdout/stderr.
It intentionally reports no speedup from this calibration.

| Software | Current executable adapter | Mathematical boundary |
|---|---|---|
| RustBEARS | Python worker, implemented | Explicit tree/geography/state labels, Q, normalized C, root weights, empty row and conditioning are exported. |
| Original BioGeoBEARS | Read-only M0 oracle, implemented | Original `cladoRcpp`/`rexpokit` bindings are verified; the comparison currently targets the legacy Contract-C score. |
| Lagrange-ng | Pending; executable absent from PATH | Official docs distinguish evaluation from optimization and expose rates and worker controls. Root measure, internal empty-state transmission and leaf-survival normalization still need source/runtime verification. |
| RevBayes | Pending; executable absent from PATH | `dnPhyloCTMCClado` documents explicit Q, cladogenetic probabilities and normalized root frequencies. Their exact correspondence to the selected RustBEARS contract remains unverified. |

The capability entries are based on the official [Lagrange-ng repository
documentation](https://github.com/computations/lagrange-ng) and
[RevBayes distribution reference](https://revbayes.github.io/documentation/dnPhyloCTMCClado.html).
Neither pending entry is a generic command template disguised as a working
adapter. Installing a binary will not automatically mark its model aligned.

The alignment gate freezes one input tree and branch-time unit, one geographic
matrix and area order, the full allowed-state ordering, every Q entry and every
normalized C entry, the root measure including its scale, the empty-state split
rule, and the likelihood conditioning event. The concrete scalar check requires
Q/C differences ≤ 1e-12 and score differences ≤ 1e-8. These are small-case
calibration tolerances, not a universal equivalence theorem. In particular,
comparing a leaf-survival-normalized likelihood with an unnormalized legacy
score is prohibited. A parameter-dependent normalization cannot be repaired by
one fitted constant.

Evaluation and optimization are separate experiment types. `--phase evaluate`
measures repeated fixed-parameter public likelihood calls after one warmup.
`--phase optimization` executes each real default optimizer but labels their
work policies as unmatched, so it also reports no optimizer speedup. A future
matched optimization benchmark must predeclare common free parameters, boxes,
initial points, convergence tolerances, endpoint likelihood tolerance, retries,
and CPU/BLAS budgets. Different multi-start policies may be compared as complete
workflows only when their resulting objective quality and failure rates are
reported alongside runtime.

Lagrange-ng and RevBayes need concrete adapters that export or independently
verify these model components, followed by exact small-tree checks, before
their timings enter a ratio. Inference-engine and process-startup timings must
not be mixed. Only runs on the same host and declared thread budget can be
paired. Failure, unsupported input or unknown alignment yields a missing
speedup value with a reason, never a zero time or an inferred performance win.

Historical stages 04–14 are read-only. All new inputs, outputs and logs go under
stage 19. This local work installs no software and launches no large benchmark.
