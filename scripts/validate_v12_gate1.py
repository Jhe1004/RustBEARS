#!/usr/bin/env python3
"""Independent Tier-1a probability-contract checks and simulation experiments.

The reference uses SciPy matrix exponentials and an independent dense pruning
recursion.  It reuses ONLY the already validated v1.1.4 cladogenetic weight
operator, not the new likelihood context or sampler.  This is deliberately
stated in the result manifest: it is not an independent rewrite of every DEC
cladogenetic convention.

Full settings are six configurations, N=1,000,000 and R=9,999, with ONE Holm
family of twelve tests.  A maximum-TVD summary is diagnostic only.  Fault
experiments use declared test-side mutant distributions, quantify their exact
distance first, and report binomial intervals without selecting random seeds.
No historical artifact is written.  Smoke output cannot claim Gate 1 passed.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.linalg import expm
from scipy.stats import beta

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

import rustbears as b
from rustbears._ffi import kernel
from rustbears.operators import _states


@dataclass(frozen=True)
class Case:
    name: str
    tree: str
    areas: int
    contract: str
    d: float
    e: float
    j: float

    @property
    def model_name(self):
        return "DEC+J" if self.j > 0 else "DEC"


TREE4 = "((A:1.1,B:1.1):1.2,(C:0.8,D:0.8):1.5);"
TREE5 = "(((A:0.7,B:0.7):0.6,C:1.3):1.1,(D:0.9,E:0.9):1.5);"
CORE_CASES = tuple(
    Case(f"{short}_e{e:g}_{'DECJ' if j else 'DEC'}", tree, areas, contract, 0.11, e, j)
    for short, tree, areas, contract in (
        ("B_2area4tip", TREE4, 2, "unconditional_complete"),
        ("A_3area5tip", TREE5, 3, "tip_survival"),
    )
    for e, j in ((0.0, 0.0), (0.05, 0.15), (0.2, 0.0))
)


def make_model(case: Case, *, root_weights=None):
    tree = b.read_newick(case.tree)
    geog = b.TipRanges(tuple("XYZ"[:case.areas]), {name: (0,) for name in tree.tip_labels})
    prior = b.RootPriorSpec() if root_weights is None else b.RootPriorSpec(
        mode="custom_weights", custom_weights=root_weights)
    model = b.BioGeoBEARS(
        tree, geog, case.model_name, likelihood_conditioning=case.contract,
        root_prior=prior, allow_null_tips=case.contract == "unconditional_complete",
        precision="strict_f64",
    )
    params = model.model.resolve({"d": case.d, "e": case.e, "j": case.j})
    return model, params


class ReferenceTree:
    """Small-tree probability law, separate from the production downpass."""

    def __init__(self, model, params, root_weights=None, *, complete_empty=False,
                 empty_restart=False, d_multiplier=1.0):
        self.states = tuple(tuple(s) for s in model.states)
        self.nstates = len(self.states)
        self.ntips = model.tree.tip_count
        self.root = int(model.tree.root)
        self.steps = model.tree.edges.reshape(-1, 2, 2)[:, :, 1]
        self.parents = model.tree.edges.reshape(-1, 2, 2)[:, 0, 0]
        self.children = {int(parent): (int(left), int(right))
                         for parent, (left, right) in zip(self.parents, self.steps)}
        self.nnode = len(model.tree.labels)
        self.nonnull = np.array([bool(s) for s in self.states])
        self.null = self.states.index(())
        if root_weights is None:
            root_weights = self.nonnull.astype(float)
            root_weights /= root_weights.sum()
        self.root_weights = np.asarray(root_weights, float).copy()
        if self.root_weights.shape != (self.nstates,):
            raise ValueError("reference root vector shape mismatch")

        # DEC anagenetic transitions are written here independently of Rust.
        dmat, elist, amat, weights, jump, *priors = model._layer_arrays(0, params)
        if np.any(amat != 0):
            raise ValueError("reference is restricted to DEC/DEC+J without a")
        lookup = {frozenset(s): i for i, s in enumerate(self.states)}
        q = np.zeros((self.nstates, self.nstates))
        for i, source in enumerate(self.states):
            if not source:
                continue
            occupied = frozenset(source)
            for area in range(model.na):
                if area in occupied:
                    target = lookup.get(occupied - {area})
                    rate = elist[area]
                else:
                    target = lookup.get(occupied | {area})
                    rate = d_multiplier * dmat[list(source), area].sum()
                if target is not None:
                    q[i, target] += rate
            q[i, i] = -q[i].sum()
        self.q = q
        self.transition = {}
        for _, child in model.tree.edges:
            child = int(child)
            transition = expm(q * model.tree.lengths[child])
            if np.min(transition) < -1e-13:
                raise ArithmeticError("negative reference transition probability")
            transition[transition < 0] = 0
            self.transition[child] = transition

        # Use the old local cladogenetic event operator only.  The new contract
        # implementation, root correction and whole-tree recursion are not used.
        nonnull_indices = np.flatnonzero(self.nonnull)
        local_states = tuple(self.states[i] for i in nonnull_indices)
        events = kernel(4, *_states(local_states), weights, jump, *priors,
                        model.na, [], precision="strict_f64")
        clado = np.zeros((self.nstates,) * 3)
        for parent, left, right, value in zip(*events):
            clado[nonnull_indices[int(parent)], nonnull_indices[int(left)],
                  nonnull_indices[int(right)]] += value
        for i in nonnull_indices:
            total = clado[i].sum()
            if total <= 0:
                raise ArithmeticError("nonempty cladogenetic row has zero mass")
            clado[i] /= total
        if complete_empty:
            clado[self.null, self.null, self.null] = 1.0
        if empty_restart:
            # Explicit test-side mutant: an empty internal lineage is wrongly
            # restarted in a nonempty singleton instead of remaining empty.
            clado[self.null] = 0.0
            singleton = next(i for i, s in enumerate(self.states) if s == (0,))
            clado[self.null, singleton, singleton] = 1.0
        self.clado = clado

    def downpass(self, tips):
        tips = np.asarray(tips, float)
        if tips.shape != (self.ntips, self.nstates):
            raise ValueError("reference tip likelihood shape mismatch")
        down = np.zeros((self.nnode, self.nstates))
        down[:self.ntips] = tips
        for parent, (left, right) in zip(self.parents, self.steps):
            lmsg = self.transition[int(left)] @ down[int(left)]
            rmsg = self.transition[int(right)] @ down[int(right)]
            down[int(parent)] = np.einsum("ijk,j,k->i", self.clado, lmsg, rmsg)
        return down

    def pattern_likelihoods(self, patterns, *, conditioning="none", batch_size=2048,
                            omit_survival_root_weight=False):
        patterns = np.asarray(patterns, dtype=np.int64)
        if patterns.ndim != 2 or patterns.shape[1] != self.ntips:
            raise ValueError("pattern array shape mismatch")
        root_weights = self.root_weights.copy()
        survival = self.downpass(np.tile(self.nonnull, (self.ntips, 1)))[self.root]
        normalizer = float(root_weights @ survival) if conditioning == "tip_survival" else 1.0
        if conditioning == "tip_survival" and normalizer <= 0:
            raise ValueError("survival probability is zero")
        if omit_survival_root_weight:
            if conditioning != "tip_survival":
                raise ValueError("root-weight mutant requires survival conditioning")
            active = root_weights > 0
            if np.any(survival[active] <= 0):
                raise ValueError("fault root support includes impossible survival")
            root_weights[active] /= survival[active]
            normalizer = 1.0
        answer = np.empty(len(patterns))
        eye = np.eye(self.nstates)
        for begin in range(0, len(patterns), batch_size):
            block = patterns[begin:begin + batch_size]
            down = {i: eye[block[:, i]] for i in range(self.ntips)}
            for parent, (left, right) in zip(self.parents, self.steps):
                lmsg = down[int(left)] @ self.transition[int(left)].T
                rmsg = down[int(right)] @ self.transition[int(right)].T
                down[int(parent)] = np.einsum("ijk,rj,rk->ri", self.clado, lmsg, rmsg,
                                               optimize=True)
            answer[begin:begin + len(block)] = down[self.root] @ root_weights / normalizer
        return answer


def enumerate_patterns(allowed_states, ntips):
    return np.asarray(list(itertools.product(tuple(allowed_states), repeat=ntips)), dtype=np.int64)


def pattern_counts(tip_states, allowed_states):
    values = np.asarray(tip_states)
    allowed = np.asarray(allowed_states, dtype=np.int64)
    if values.ndim != 2 or not np.issubdtype(values.dtype, np.integer):
        raise ValueError("sampler must return an integer tip-state matrix")
    lookup = {int(state): rank for rank, state in enumerate(allowed)}
    flat = values.ravel()
    if not set(map(int, np.unique(flat))).issubset(lookup):
        raise ValueError("sampler produced a state outside the tested observation space")
    digits = np.array([lookup[int(v)] for v in flat], dtype=np.int64).reshape(values.shape)
    code = np.zeros(len(values), dtype=np.int64)
    for column in digits.T:
        code = code * len(allowed) + column
    return np.bincount(code, minlength=len(allowed) ** values.shape[1])


def holm_adjust(pvalues):
    p = np.asarray(pvalues, float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("p values must be a finite probability vector")
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum(1.0, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    result = np.empty_like(p)
    result[order] = adjusted
    return result


def clopper_pearson(successes, trials, alpha=0.05):
    if not 0 <= successes <= trials or trials <= 0:
        raise ValueError("invalid binomial counts")
    lower = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, trials - successes + 1))
    upper = 1.0 if successes == trials else float(beta.ppf(1 - alpha / 2, successes + 1, trials - successes))
    return [lower, upper]


def bin_mapping(probabilities, samples, min_expected=10):
    """Deterministically combine low-probability cells using only null p."""
    p = np.asarray(probabilities, float)
    mapping = np.full(len(p), -1, dtype=np.int64)
    positive = np.flatnonzero(p > 0)
    groups = []
    current, mass = [], 0.0
    for index in positive[np.argsort(p[positive], kind="stable")]:
        current.append(int(index))
        mass += p[index]
        if samples * mass >= min_expected:
            groups.append(current)
            current, mass = [], 0.0
    if current:
        if groups:
            groups[-1].extend(current)
        else:
            groups.append(current)
    for group_id, indices in enumerate(groups):
        mapping[indices] = group_id
    return mapping, np.asarray([p[indices].sum() for indices in groups])


def goodness_statistics(counts, p, bins, bin_probabilities, samples):
    counts = np.asarray(counts)
    one = counts.ndim == 1
    matrix = counts[None, :] if one else counts
    if matrix.shape[1] != len(p) or np.any(matrix < 0) or np.any(matrix.sum(axis=1) != samples):
        raise ValueError("counts do not match expected size")
    tvd = np.abs(matrix / samples - p).sum(axis=1) / 2
    structural = (matrix[:, bins < 0] > 0).any(axis=1)
    grouped = np.zeros((len(matrix), len(bin_probabilities)))
    for row, values in enumerate(matrix):
        grouped[row] = np.bincount(bins[bins >= 0], weights=values[bins >= 0],
                                   minlength=len(bin_probabilities))
    expected = samples * bin_probabilities
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(grouped > 0, grouped * np.log(grouped / expected), 0)
    g2 = 2 * terms.sum(axis=1)
    g2[structural] = np.inf
    result = np.column_stack((tvd, g2))
    return result[0] if one else result


def ideal_reference(p, samples, repetitions, rng, batch_size=32):
    bins, bin_p = bin_mapping(p, samples)
    stats = np.empty((repetitions, 2))
    for start in range(0, repetitions, batch_size):
        size = min(batch_size, repetitions - start)
        draws = rng.multinomial(samples, p, size=size)
        stats[start:start + size] = goodness_statistics(draws, p, bins, bin_p, samples)
    return bins, bin_p, stats


def monte_carlo_p(observed, reference):
    return (1 + (np.asarray(reference) >= observed).sum(axis=0)) / (len(reference) + 1)


def native_pattern_probabilities(model, params, patterns):
    """Exercise the production context while replacing only observed tip data."""
    result = np.empty(len(patterns))
    eye = np.eye(len(model.states))
    with model._context(1) as context:
        for i, pattern in enumerate(patterns):
            lnlik = model._evaluate(context, params, False, tip_likelihoods=eye[pattern])
            result[i] = math.exp(lnlik)
    return result


def assert_probability_vector(p, tolerance=1e-12):
    if not np.isfinite(p).all() or np.min(p) < 0 or abs(float(np.sum(p)) - 1) > tolerance:
        raise AssertionError(f"invalid probability vector: min={np.min(p)} sum={np.sum(p)}")


def pure_extinction_check():
    case = Case("pure_extinction", "((A:1,B:1):3,C:4);", 1,
                "unconditional_complete", 0.0, 0.2, 0.0)
    result = {}
    for contract in ("legacy_unconditioned", "unconditional_complete", "tip_survival"):
        local = Case(**{**asdict(case), "contract": contract})
        model, params = make_model(local, root_weights={(0,): 1.0})
        allowed = range(len(model.states)) if contract != "tip_survival" else [1]
        patterns = enumerate_patterns(allowed, model.tree.tip_count)
        observed = native_pattern_probabilities(model, params, patterns)
        expected_total = math.exp(-0.6) if contract == "legacy_unconditioned" else 1.0
        error = abs(float(observed.sum()) - expected_total)
        if error >= 1e-14:
            raise AssertionError(f"pure extinction {contract}: error={error}")
        result[contract] = {"sum": float(observed.sum()), "expected_sum": expected_total,
                            "absolute_error": error}
        if contract == "unconditional_complete":
            all_occupied = np.all(patterns != 0, axis=1)
            survival = float(observed[all_occupied].sum())
            if abs(survival - math.exp(-1.8)) >= 1e-14:
                raise AssertionError("pure extinction nonempty-tip survival mismatch")
            result[contract]["tip_survival_probability"] = survival
    return result


def fault_distribution(case, model, params, patterns, mode):
    complete = case.contract == "unconditional_complete"
    ref = ReferenceTree(model, params, complete_empty=complete,
                        empty_restart=mode == "empty_internal_restart",
                        d_multiplier=1.02 if mode == "d_plus_two_percent" else 1.0)
    return ref.pattern_likelihoods(
        patterns, conditioning="tip_survival" if case.contract == "tip_survival" else "none",
        omit_survival_root_weight=mode == "omit_root_survival_weight")


def run_fault_experiments(case_records, samples, repeats, seed):
    """Power of the fixed twelve-test gate against predeclared mutant laws.

    Each repetition substitutes one mutant configuration while the other five
    follow their exact null laws.  This evaluates the actual Holm family; the
    reference null arrays are reused, not regenerated until a seed passes.
    """
    targets = (
        ("d_plus_two_percent", "B_2area4tip_e0.05_DECJ"),
        ("omit_root_survival_weight", "A_3area5tip_e0.2_DEC"),
        ("empty_internal_restart", "B_2area4tip_e0.2_DEC"),
    )
    reports = []
    for domain, (mode, name) in enumerate(targets, 100):
        selected = next(r for r in case_records if r["case"].name == name)
        mutant = fault_distribution(selected["case"], selected["model"], selected["params"],
                                    selected["patterns"], mode)
        assert_probability_vector(mutant)
        exact_tvd = float(np.abs(mutant - selected["p"]).sum() / 2)
        if exact_tvd <= 1e-10:
            raise AssertionError(f"predeclared fault {mode} is not observable in {name}")
        rng = np.random.default_rng(np.random.SeedSequence([seed, domain]))
        rejected = 0
        for _ in range(repeats):
            pvalues = []
            for record in case_records:
                generating = mutant if record is selected else record["p"]
                counts = rng.multinomial(samples, generating)
                stat = goodness_statistics(counts, record["p"], record["bins"], record["bin_p"], samples)
                pvalues.extend(monte_carlo_p(stat, record["ideal"]))
            rejected += bool(np.any(holm_adjust(pvalues) <= 0.01))
        reports.append({"fault": mode, "configuration": name,
                        "implementation": "independent_test_side_mutant_probability_law",
                        "exact_tvd": exact_tvd, "trials": repeats, "rejections": rejected,
                        "observed_power": rejected / repeats,
                        "power_ci95_clopper_pearson": clopper_pearson(rejected, repeats),
                        "hard_power_threshold_applied": False,
                        "interpretation": "Sensitivity estimate; no rerun-until-pass or claim of guaranteed 99% power."})
    return reports


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(path)


def run(args):
    destination = Path(args.output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    native_path = b._ffi._resolve_native_library_path()
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPOSITORY,
                              check=True, text=True, capture_output=True).stdout.strip()
    source_files = sorted([
        *REPOSITORY.joinpath("rustbears").rglob("*.py"),
        *REPOSITORY.joinpath("rust", "crates").rglob("*.rs"),
        Path(__file__).resolve(), REPOSITORY / "tests" / "test_v12_contracts.py",
    ])
    source_hashes = {str(path.relative_to(REPOSITORY)): sha256(path) for path in source_files}
    report = {"schema": "rustbears-tier1a-gate1-v1", "status": "running",
              "settings": vars(args), "git_head": revision, "python": sys.version,
              "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__,
              "source_sha256": source_hashes,
              "thread_environment": {key: os.environ.get(key) for key in (
                  "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
              "harness_sha256": sha256(__file__), "native_library": str(native_path),
              "native_sha256": sha256(native_path),
              "reference_boundary": "Independent DEC Q, SciPy expm and dense pruning; reused v1.1.4 local clado weights.",
              "family": "ONE Holm adjustment over six configurations times TVD/G2 (12 tests), alpha=0.01",
              "maximum_tvd_role": "diagnostic_only_no_additional_rejection_gate",
              "historical_freeze_checks": "Owned by root execution manifest, not inferred here.",
              "cases": []}
    output = destination / "gate1_results.json"
    started = time.monotonic()
    write_json(output, report)
    records = []
    try:
        report["pure_extinction"] = pure_extinction_check()
        for index, case in enumerate(CORE_CASES):
            case_started = time.monotonic()
            model, params = make_model(case)
            allowed = np.arange(len(model.states))
            if case.contract == "tip_survival":
                allowed = np.flatnonzero([bool(s) for s in model.states])
            patterns = enumerate_patterns(allowed, model.tree.tip_count)
            ref = ReferenceTree(model, params, complete_empty=case.contract == "unconditional_complete")
            p = ref.pattern_likelihoods(patterns, conditioning="tip_survival" if case.contract == "tip_survival" else "none")
            assert_probability_vector(p)
            native = native_pattern_probabilities(model, params, patterns)
            assert_probability_vector(native)
            error = float(np.max(np.abs(native - p)))
            if error >= 1e-12:
                raise AssertionError(f"{case.name}: native/reference max error {error}")
            # Only absorb floating-point summation roundoff before multinomial.
            p = p / p.sum()
            simulation = b.simulate_datasets(model, params=params, n_reps=args.samples,
                                            seed=args.seed + index, num_threads=args.threads)
            counts = pattern_counts(simulation.tip_states, allowed)
            if counts.sum() != args.samples:
                raise AssertionError("native sample count mismatch")
            rng = np.random.default_rng(np.random.SeedSequence([args.seed, 1, index]))
            bins, bin_p, ideal = ideal_reference(p, args.samples, args.reference_reps, rng)
            observed = goodness_statistics(counts, p, bins, bin_p, args.samples)
            pvalues = monte_carlo_p(observed, ideal)
            np.savez_compressed(destination / f"{case.name}.npz", probabilities=p,
                                native_probabilities=native, counts=counts,
                                ideal_statistics=ideal, allowed_states=allowed)
            record = {"case": case, "model": model, "params": params, "patterns": patterns,
                      "p": p, "bins": bins, "bin_p": bin_p, "ideal": ideal}
            records.append(record)
            report["cases"].append({**asdict(case), "patterns": len(p),
                "native_probability_sum": float(native.sum()),
                "reference_probability_sum": float(p.sum()), "max_native_reference_error": error,
                "samples": int(counts.sum()), "tvd": float(observed[0]), "g2": float(observed[1]),
                "p_tvd": float(pvalues[0]), "p_g2": float(pvalues[1]),
                "seconds": time.monotonic() - case_started})
            write_json(output, report)
            print(json.dumps({"configuration": case.name, "max_error": error,
                              "p_values": pvalues.tolist(), "seconds": time.monotonic() - case_started}), flush=True)
        adjusted = holm_adjust([r[k] for r in report["cases"] for k in ("p_tvd", "p_g2")])
        for row, pair in zip(report["cases"], adjusted.reshape(-1, 2)):
            row["holm_p_tvd"], row["holm_p_g2"] = map(float, pair)
        report["simulator_rejected"] = bool(np.any(adjusted <= 0.01))
        report["maximum_tvd_diagnostic"] = max(r["tvd"] for r in report["cases"])
        if args.fault_repeats:
            report["fault_sensitivity"] = run_fault_experiments(records, args.samples, args.fault_repeats, args.seed)
        else:
            report["fault_sensitivity"] = {"status": "not_run"}
        full = args.samples >= 1_000_000 and args.reference_reps >= 9_999 and args.fault_repeats >= 100
        report["full_requested_gate_sizes_completed"] = full
        report["source_files_unchanged_during_run"] = all(
            sha256(REPOSITORY / name) == digest for name, digest in source_hashes.items())
        report["native_library_unchanged_during_run"] = sha256(native_path) == report["native_sha256"]
        report["status"] = "rejected" if report["simulator_rejected"] else (
            "full_experiments_completed_no_distribution_rejection" if full else "smoke_only_not_gate1_acceptance")
        report["seconds"] = time.monotonic() - started
        write_json(output, report)
        return 1 if report["simulator_rejected"] else 0
    except Exception as error:
        report["status"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        report["seconds"] = time.monotonic() - started
        write_json(output, report)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--samples", type=int, default=1_000_000)
    parser.add_argument("--reference-reps", type=int, default=9_999)
    parser.add_argument("--fault-repeats", type=int, default=100)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20261115)
    args = parser.parse_args()
    if min(args.samples, args.reference_reps, args.threads) < 1 or args.fault_repeats < 0 or args.seed < 0:
        parser.error("counts/threads must be positive, fault repeats and seed nonnegative")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
