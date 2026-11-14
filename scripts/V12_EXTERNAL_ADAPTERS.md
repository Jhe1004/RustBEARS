# Pinned external scalar adapters

`validate_v12_external_adapters.py` runs actual Lagrange-ng and RevBayes
executables on a two-area, three-tip tree at three fixed rate pairs. It performs
no optimization or MCMC and reports **no speedup**. It is separate from the
original BioGeoBEARS calibration runner.

The official macOS arm64 releases used here are
[Lagrange-ng v0.7.2-7](https://github.com/computations/lagrange-ng/releases/tag/v0.7.2-7)
and [RevBayes v1.4.2](https://github.com/revbayes/revbayes/releases/tag/v1.4.2).
Their binaries, pinned source archives, release metadata and download hashes are
stored under stage 19, `runtime/` and `external_sources/`. No historical runtime
or system installation was changed. Lagrange-ng's required `libomp.dylib` was
copied into its isolated `lib/`; the launcher sets `DYLD_LIBRARY_PATH` to that
directory. RevBayes's official binary is self-contained on this host.

Run from the repository root:

```sh
PYTHONPATH=. /opt/anaconda3/bin/python scripts/validate_v12_external_adapters.py \
  --output-dir ../19_RustBEARS_v1_2_开发与验证_20261115/evidence/external_scalar_adapters
```

Each case retains its generated control/Rev scripts, exact inputs, native
stdout/stderr, executable hashes, model-component NPZ and machine-readable
comparison. The report's `process_seconds` includes startup, parsing and output;
it is not an isolated likelihood benchmark.

## Verified scope

| Adapter/configuration | Result | Maximum scalar difference |
|---|---|---:|
| RevBayes `Include`, pattern-weighted DEC, normalized nonempty root | Aligned with RustBEARS Contract C using the same normalized root | 1.78e-15 |
| Lagrange-ng native DEC versus RustBEARS standard DEC | Unsupported: different cladogenetic probabilities | Not an equivalence comparison |
| Lagrange-ng native DEC versus explicit **two-area legacy scalar `v=2, s=y=1`** variant | Aligned only for this declared configuration | 2.00e-9 |
| RevBayes default `CondSurv` versus whole-tree Contract A | Unsupported: different conditioning events | Observed differences 0.0062–0.0167 |

The fixed pairs are `(d,e) = (0.2,0.1), (0.05,0.02), (0.4,0.2)` on
`((T1:1,T2:1):1,T3:2);`, with T1/T2 in A and T3 in B. No claim extends these
results to time strata, other modifiers, larger state spaces or optimization.

## RevBayes contract

The adapter explicitly requests `fnDECRateMatrix(nullRange="Include")` and
`fnDECCladoProbs(probType="pattern", eventTypes=["s","a"],
eventProbs=simplex(1,1))`. Root frequencies are `simplex(0,1,1,1)`, so the
RustBEARS model receives `RootPriorSpec("uniform_nonnull")` under Contract C.
This chooses the same root measure at model construction; no score offset is
estimated or added afterward.

Q and the unit-time transition matrix are exported through public member
functions at full precision. **All 64 entries of C are independently measured
at runtime:** a zero-length two-tip tree with a point-mass ancestral root and
fixed daughter states has likelihood `C[a,l,r]`. This also verifies that the
entire null ancestral row is zero. Maximum Q/C/P differences were respectively
0, 2.78e-17 and 1.11e-16. A separate SciPy exponential and tensor-pruning
calculation checks the fixed-tree score.

RevBayes's display formatter for a cladogenetic matrix rounds to four decimals;
those display files are retained for readability and **never** used in the
1e-12 matrix gate. Similarly, string concatenation rounds scalar conversions;
the Rev scripts pass numerical objects directly to `write` with CLI
`-o outputPrecision=17`. CLI options do not update the user's `.RevBayes.ini`.

The default `CondSurv` renormalizes each nonempty row of each branch transition
matrix, dividing by that branch's nonextinction probability. This is not the
whole-tree leaf-survival normalizer in Contract A. The adapter executes this
default as a negative control and records its nonzero difference. Source:
[RateMatrix_DECRateMatrix.cpp](https://github.com/revbayes/revbayes/blob/v1.4.2/src/core/datatypes/phylogenetics/ratematrix/RateMatrix_DECRateMatrix.cpp).

## Lagrange-ng contract and limitations

The release's split generator loops over each occupied singleton and appends
both daughter orders. For the AB ancestor it includes each A|B vicariance event
twice: each ordered vicariance has probability 1/4, while each of four ordered
subset events has probability 1/8. Standard RustBEARS DEC gives six probabilities
of 1/6. The maximum C difference is 1/12. This is a model difference, not a
numerical tolerance issue. Source:
[Operation.cpp](https://github.com/computations/lagrange-ng/blob/v0.7.2-7/src/Operation.cpp).

For **two areas only**, explicit legacy scalar weights `v=2, s=y=1` reproduce
this C exactly. This variant uses the preserved historical scalar API and is
not an unlocked Tier-1b uncertainty model. It must not be advertised as a
standard-DEC speed comparison or generalized to three or more areas.

The rate generator is verified against the pinned source: a gain sums the
dispersal rates from occupied source areas, a loss has extinction rate e, the
empty row is absorbing, and each diagonal is minus its off-diagonal row sum.
Native state IDs use the opposite two-area bit order; the explicit permutation
to Rust order is `[0,2,1,3]`. The normalized split tensor is source-derived and
stored alongside the native output; unlike the RevBayes C audit, it is not a
separate native matrix-export API.

Native root weights are raw ones on all four states. The null root conditional
likelihood is identically zero because the null cladogenetic row is zero.
Thus the effective nonempty measure is exactly the same raw unit measure used
by the historical Rust scalar. This is an exact zero contribution, not a fitted
constant correction. Contract B's completed empty row is not claimed.

This release resets direct evaluation rates to 0.01 in
[`Workspace::getDefaultParams`](https://github.com/computations/lagrange-ng/blob/v0.7.2-7/src/Workspace.hpp).
The adapter therefore supplies the rates through the release's official
[checkpoint resume interface](https://github.com/computations/lagrange-ng/blob/v0.7.2-7/src/Checkpoint.hpp),
records that mechanism, and verifies both rates in native JSON before comparing.
It does not modify the external executable or source.

Standard stdout prints only seven significant digits. Full precision total
likelihood is recovered from the root-state JSON as `joint_log_likelihood -
log(posterior_state_probability)`; every exported root state must give the same
normalizer within 1e-12. This identity follows the release's
[JSON normalization code](https://github.com/computations/lagrange-ng/blob/v0.7.2-7/src/IO.cpp)
and is also checked against independent SciPy pruning.

The observed Lagrange differences up to 2.00e-9 pass the predeclared 1e-8 scalar
gate. Neither this small numerical check nor any unsupported default
configuration authorizes a performance ratio.
