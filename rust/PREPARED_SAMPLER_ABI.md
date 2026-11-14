# Prepared whole-tree sampler ABI

The v1.2 sampler API is additive. The existing `bgb_kernel(26, ...)` adapter
retains the v1.1.4 state order, random stream, and event layout. New uncertainty
workflows use the prepared API to avoid rebuilding geometry, generators,
cladogenetic rows, branch messages, and transition-row caches for every batch.

## Preparation and ownership

```c
SamplerHandle *bgb_sampler_prepare(
    const double **buffers, const size_t *lengths, size_t count,
    const double *root_weights, size_t root_count, uint32_t conditioning);
const char *bgb_sampler_error(const SamplerHandle *);
void bgb_sampler_release(SamplerHandle *);
```

All preparation buffers and root weights are copied before preparation returns;
callers may release or modify their arrays afterward. A handle owns its tree,
model matrices, downward messages, branch targets, and lazy transition rows.
The returned handle must be checked with `bgb_sampler_error`, including after
failed preparation, and released once. Error strings remain borrowed from the
handle until release. A null release is safe.

Wrappers own handles uniquely, serialize use and close with the same lock, and
clear their pointer when closing. Repeated wrapper `close()` calls are safe;
calling the raw release function twice on the same non-null pointer is invalid.
The Rust mutex serializes operations on a live handle but cannot protect a
pointer after a caller has released it. Handles cannot be copied or pickled.

Input buffers are the existing whole-tree BSM geometry with a two-item header:

| Index | Contents |
|---|---|
| 0 | `[tip_count, root_node]`, positive tip count and zero-based root index |
| 1, 2 | Packed states and lengths, in the existing enumerator order |
| 3 | Tip likelihoods, column-major `tip_count × state_count` |
| 4 | Pruning-order triples `[parent, left, right]` |
| 5 | Node ages |
| 6 | Layer index for each node |
| 7 | Fossil-hook flags |
| 8 | Node likelihood constraints, column-major `node_count × state_count` |
| 9 | Segment records `[child_node, layer, duration, start_age]` |
| 10 onward | Per layer: column-major Q, allowed-state mask, flat `[parent,left,right,weight]` cladogenesis records, column-major dispersal matrix |

New prepared modes currently support the tested single-layer fixed state space,
without fossil hooks or internal-node constraints. Empty-state transmission is
added internally as `C(0;0,0)=1` and must not also be supplied in the event buffer.
This mathematical transmission is omitted from biological cladogenesis counts.

## Conditioning and root probabilities

| Mode | Meaning | Root sampling weights |
|---|---|---|
| 0 | BSM conditional on the supplied observed tips | `root_weights × observed downward message` |
| 1 | Simulate datasets conditional on every tip being nonempty | `root_weights × root survival message` |
| 2 | Simulate the unconditional complete absorbing process | Supplied root probabilities; all-one observation messages |
| 4 | Simulate surviving tips from an already conditioned root distribution | Supplied probabilities directly; descendant messages still enforce survival |

Mode 3 belongs to the private legacy adapter and is rejected by public
preparation. Modes 1, 2, and 4 require normalized, finite, nonnegative root
probabilities and the empty state at index zero. Modes 1 and 4 require zero
empty-root mass. Mode 4 also rejects positive mass on states with zero survival
probability. It is used for an explicitly requested observed-data root
posterior; multiplying that posterior by survival again would change the
requested root distribution.

## Drawing and results

```c
Reply *bgb_sampler_draw_bsm(SamplerHandle *, uint64_t start, size_t count,
    uint64_t seed, size_t threads, uint32_t rejection);
Reply *bgb_sampler_simulate_tips(SamplerHandle *, uint64_t start, size_t count,
    uint64_t seed, size_t threads);
```

`start` is one-based and identifies the first logical map/replicate. Calls allow
1–1,000,000 samples, 1–256 threads, and final IDs no greater than `2^53−1`, since
event tables represent IDs as doubles. Seeds are full unsigned 64-bit values.
R wrappers must use decimal strings for values that cannot be represented
exactly by an R double. Every `Reply` is checked and freed with the existing
result API, including error replies.

BSM returns the existing flat cladogenesis and anagenesis tables. Simulation
returns two row-major arrays: `count × tip_count` and `count × node_count`.
State indices are zero-based and follow the original state order. Simulation
handles reject BSM draws; observed-data handles reject dataset simulation.

Logical replicate IDs and the seed determine the random stream. Splitting a
draw into batches or changing thread count preserves the same individual
replicates. The simulation domain is separate from the legacy BSM stream.
New multi-model/tree/bootstrap workflows derive task seeds in
`rustbears.uncertainty.rng`: SHA-256 covers all semantic task fields and a
versioned domain. A 64-bit result is not guaranteed collision-free; callers
declare batches with duplicate-key and truncated-seed collision checks. A
checkpoint retry reuses its original identity and seed.

Prepared handles retain cached transition rows across calls, but outputs are
owned independently by each reply. Concurrent draws on separate handles are
independent. Worker threads call no Python/R interpreter APIs.

## Validation

Native tests exercise already-conditioned root frequencies under unequal
survival probabilities, buffer ownership, mode rejection, and reply cleanup.
The independent Python probability gate enumerates all 256 four-tip patterns
for the complete two-area model and all 81 nonempty patterns for the survival
contract across six model families and the tested `x`/`w` modifiers. Legacy
event fixtures separately check unchanged events, thread counts, and batching.
