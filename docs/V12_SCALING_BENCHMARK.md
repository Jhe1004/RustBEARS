# Fixed-parameter scaling benchmark

`scripts/benchmark_v12_scaling.py` implements the same workload on macOS arm64
and Linux: start a fresh process, load one fixed tree and geography, construct
the declared model, evaluate one scalar log likelihood, and exit. Elapsed time
includes interpreter startup, imports and input parsing. Peak RSS is measured
with POSIX `wait4` for each child. These are end-to-end evaluation workloads;
the ratios are neither optimizer speedups nor isolated numerical-kernel speedups.
The P0 experiments separately measure prepared likelihood and complete fit costs.

The predeclared grid is 19, 50 and 150 tips, with 4, 16, 57 and 128 states.
Balanced trees have root age 10, and singleton tip ranges cycle through areas.
Rates are fixed at d=0.015 and e=0.008. Three repetitions are retained without
best-time selection; execution order uses a predeclared PCG64 permutation with
seed 20261115 for each cell and repetition. BLAS and native
threads are fixed at one. Load averages, CPU time, input/source/executable hashes,
timeouts and failures are retained. The wall timer polls at 2 ms; sub-0.1-second
differences are not interpreted as precise kernel timing. A timeout is a censored workload, never an
omitted slow observation. Public benchmarks must report concurrent host load.

| Comparison stratum | Root and precision | Model correspondence |
|---|---|---|
| Original BioGeoBEARS 1.1.4 vs RustBEARS | Raw unit weights on nonempty states; Rust `legacy_f32` | Standard DEC, legacy unconditioned scoring |
| RevBayes 1.4.2 vs RustBEARS | Uniform normalized nonempty root; Rust `strict_f64` | `nullRange="Include"`, `probType="pattern"`, equal subset/allopatry weights, `widespreadAllopatry=FALSE` |
| Lagrange-ng 0.7.2-7 vs RustBEARS | Effective raw nonempty unit weights; Rust `strict_f64` | **Two areas only**, explicit historical Rust v=2, s=y=1 variant; this is not default BioGeoBEARS DEC |

Every comparison uses the same zero empty-state cladogenesis row. None is a
whole-tree survival-conditioned Contract-A comparison. Lagrange-ng assigns
different default split weights; therefore no default-DEC speedup is reported
for that program. The effective null-root contribution is zero under the shared
zero cladogenesis row, so no post-hoc additive likelihood correction is used.

Complete Q/C/state-order/root checks run in a separate preflight, outside the
timed workloads. RevBayes C values are serialized through `mnFile(format="json")`:
its ordinary printed matrix is rounded to four decimals. Its character converter
exposes the full bit-state dictionary even under a range-size cap; the adapter
verifies both the allowed prefix and excluded suffix before applying the state
permutation. Component tolerance is 1e-12; every timed likelihood must agree
within 1e-8. Failed alignment suppresses the corresponding time ratio without
discarding its timings or affecting a different comparison stratum.

The actual adapter check covers all four state counts on a three-tip tree. It
does not substitute for the full size grid. Formal execution requires recorded
Gate 2 and Gate 3 passes. The default command only writes a manifest:

```sh
python scripts/benchmark_v12_scaling.py --output /absolute/path/to/new/output \
  --rscript /absolute/path/to/original/Rscript \
  --oracle-library /absolute/path/to/original/oracle-library \
  --revbayes /absolute/path/to/pinned/rb \
  --lagrange /absolute/path/to/pinned/lagrange-ng \
  --mode formal --gate-status /absolute/path/to/release_gate_status.json
```

Append `--execute` only for the released experiment. `--mode smoke --tips 3
--repeats 1 --execute` checks adapters without being labelled a formal benchmark.
Original packages must be genuine installations: the R worker explicitly rejects
the RustBEARS namespace and checks that its cladogenesis function comes from
cladoRcpp. The Lagrange adapter uses the official checkpoint rate interface and
verifies the rates actually returned, because this pinned release resets direct
evaluation configuration rates to defaults.

Primary implementation sources: [RevBayes DEC cladogenesis](https://revbayes.github.io/documentation/fnDECCladoProbs.html),
[RevBayes 1.4.2](https://github.com/revbayes/revbayes/tree/v1.4.2),
[Lagrange-ng 0.7.2-7](https://github.com/computations/lagrange-ng/tree/v0.7.2-7).
The local source archives and official release metadata are retained with SHA256
digests in the development evidence directory.
