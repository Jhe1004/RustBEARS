# Inherited R examples: scientific regression contracts

The full inherited tests run in CI (`NOT_CRAN=true`). `R CMD check` continues to
fail on errors or warnings. The original failures and assertions remain available
in commit `225409f` and [CI run 43669798991](https://github.com/Jhe1004/RustBEARS/actions/runs/43669798991).
This change repairs the tests; it does not change the probability kernels, model
parameterization, production optimizer, random streams, or tolerances of the
existing v1.2 scientific gates.

## Three-epoch example

`areas_allowed_noC2.txt` contains three masks, but the example supplied the
two-epoch `timeperiods.txt`. The existing `timeperiods_OLD.txt` contains the three
endpoints (0.1, 1.5, 100) assumed by the example's own middle-period assertions.
The test now pairs these matching inputs. No masks are discarded or invented.

Running the corrected input with every old assertion unchanged produced 102
passes and nine failures in both an independent upstream BioGeoBEARS 1.1.4
process and RustBEARS. The failures required a particular realization of 50
random maps: six directional means and the mean, standard deviation and sum of
dispersal counts. For example, the archived sum was 475, whereas these upstream
and Rust runs returned 451 and 457. These counts do not define an expected
probability distribution, and equality is not a valid random-map test.

The replacement retains all 100 original area-exclusion checks and adds checks
of every event's time, allowed ranges on both sides of the transition, one-area
addition/removal, and the two genuine bifurcations per map. It independently
recounts every event type and dispersal direction from the event tables, then
compares all reported means, standard deviations and sums to those counts.
This tests the stratified history and its summaries, without asserting that
unrelated RNG implementations produce the same paths. Independent distributional
simulation gates remain separate from this structural regression.

A fixed DEC+J point (d=0.4, e=0.2, j=0.1) evaluates to
`-18.410726319261713` in both original native kernels and Rust `LegacyF32`.
The test uses an absolute error bound of `1e-10`; input hashes and conventions
are recorded in `inst/validation/three_epoch_fixed_oracle.json`.

## Trait and six-model examples

The inherited assertions required archived optimization coordinates to match
after decimal rounding. Different valid local optima, an improved objective, or
rounding across a threshold could fail. The two BAYAREALIKE+J trait archives
themselves have different coordinates and log likelihoods. The six-model
Psychotria table similarly stored rounded coordinates and rounded objectives,
which are not an exact scalar reference pair.
For example, the archived DIVALIKE point d=0.045, e=0 evaluates to
`-33.151842841127149` in both independent upstream and Rust `LegacyF32` runs;
that value does not round to the stored `-33.1`. This fixed-point comparison
does not prove a global optimum.

All historical parameter constants remain in the examples as fixed candidates.
These are legacy compatibility evaluations; the normalization and simulation
consistency gates for contracts A and B remain separate.
The replacement uses an independent upstream R process with genuine cladoRcpp
and rexpokit kernels, pinned in CI to BioGeoBEARS commit
`876a48c4b3d6ada678e4de4d130b5b3bafea2b2c`. The runner accepts only plain input
data and a configuration whitelist. It reconstructs the upstream model and
trait S4 objects, rereads the source data, and rejects a loaded RustBEARS
namespace. No cached tree likelihoods, Q/C matrices, or native contexts cross
the process boundary. CI checks that the source and installed runner copies
are identical.

For each fitted and archived point, the gates require:

1. The stored fit score equals a fresh evaluation at its actual parameters
   (absolute error below `1e-8`), all free parameters are finite and in bounds,
   and the original search did at least as well as its own initial point.
2. A genuine optimizer result is present. Optim/optimx convergence codes are
   checked; GenSA, which does not provide such a code, must have a real search
   trace and positive evaluation count. No claim of global convergence follows.
3. Rust `LegacyF32` and the independent original process agree at exactly the
   same parameters: log likelihood within absolute `1e-8`, top and bottom
   ancestral matrices within absolute `1e-7`, with identical shapes and missing
   entries and finite non-missing entries.

An archive candidate outside the declared optimizer box is never clipped. For
fixed evaluation only, each engine's local input container expands its box to
include that exact point. The original box and eligibility flag are retained;
an ineligible point cannot become a lower bound for an optimization constrained
to the original box. Negative rate/multiplier candidates are rejected.

The original heuristic fit, including its precision-specific numerical
shortfalls from feasible archive points, is preserved in the evidence. It is never replaced by a
polished fit and is not claimed to reproduce an old optimization trajectory.
An additional, separate bounded regression starts BOBYQA at every eligible
fitted/archive point under `StrictF64`. Its coordinate scale is
`max(abs(start), 0.01)`, initial trust radius is the smaller of 0.1 and one fifth
of the narrowest scaled box, final radius is one millionth of that radius,
and each start has a 5,000-evaluation cap. All outcomes are saved; the best
actually converged result must reach the best starting objective within
absolute `1e-7`. Both sides of this optimization comparison use `StrictF64`.
The selected point is then independently checked against upstream under
`LegacyF32`, including ancestral matrices. No global-optimum guarantee is made.

Earlier diagnostic attempts near fitted points returned line-search codes with
L-BFGS-B under `LegacyF32` and false-convergence codes with nlminb under
`StrictF64`. The latter failures are not attributed to legacy quantization.
Those failures are retained in the local diagnostic evidence and are not
treated as convergence. BOBYQA is a declared derivative-free regression
protocol, not a production optimizer change or a reinterpretation of those
failed runs.

Every CI oracle comparison writes input/output RDS, subprocess logs, package
versions, source/native hashes, original and separate optimizer results, and
candidate eligibility under `legacy-oracle-evidence/` in the R check artifacts.
Nested-model and independent-trait checks remain in the original examples.
The hard-coded 23-worker requests are now one explicit worker so that the
test configuration is valid on two-core CI runners. Interactive PDF viewers
are not launched by noninteractive checks.
