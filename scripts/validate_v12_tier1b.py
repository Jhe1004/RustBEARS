"""Native Tier-1b gate before enabling its public model support matrix.

The model is constructed through the legacy geometry front end, then the
additive native contract API is exercised directly. This deliberately tests
the native capability before the high-level guard is relaxed. The reference
uses SciPy expm and independent full-tree pruning from the Gate-1 harness.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import rustbears as b
from rustbears._ffi import CompiledContext, PreparedSampler, _resolve_native_library_path
from rustbears.uncertainty.rng import TaskIdentity, derive_task_seeds, task_fingerprint
from rustbears.uncertainty.simulate import sampler_buffers
from scripts.validate_v12_gate1 import ReferenceTree, enumerate_patterns, pattern_counts

FAMILIES = ("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
MODIFIERS = ("base", "x", "w", "xw")
TREE = "((A:0.3,B:0.7):0.4,(C:0.6,D:0.4):0.5);"
GEOG = {"A": "10", "B": "01", "C": "11", "D": "10"}


def run(samples=100_000):
    identities = [TaskIdentity("tier1b_validation", family+":"+modifier, "four_tip",
                               role=contract)
                  for family in FAMILIES for modifier in MODIFIERS
                  for contract in ("tip_survival", "unconditional_complete")]
    seeds = derive_task_seeds(20261115, identities)
    records = []
    for task in identities:
        family, modifier = task.model_id.split(":")
        changes = {}
        if "x" in modifier:
            changes["x"] = dict(type="free", init=-0.7)
        if "w" in modifier:
            changes["w"] = dict(type="free", init=1.3)
        template = b.model_template(family).with_parameters(**changes)
        model = b.BioGeoBEARS(TREE, GEOG, template,
            distances=[[1.0, 2.0], [4.0, 1.0]],
            dispersal_multipliers=[[1.0, 0.7], [1.6, 1.0]], precision="strict_f64")
        params = model.model.resolve({"d": 0.13, "e": 0.2, "j": 0.17 if "+J" in family else 0.0})
        root = np.array([0.0, 0.1, 0.3, 0.6])
        reference = ReferenceTree(model, params, root_weights=root, complete_empty=True)
        allowed = np.arange(1, 4) if task.role == "tip_survival" else np.arange(4)
        patterns = enumerate_patterns(allowed, model.tree.tip_count)
        expected = reference.pattern_likelihoods(patterns,
            conditioning="tip_survival" if task.role == "tip_survival" else "none")
        geometry = list(model._geometry)
        geometry[7] = root
        arguments = model._layer_buffers(0, params, False, False)
        native = []
        with CompiledContext(geometry, num_threads=1) as context:
            for pattern in patterns:
                result = context.evaluate(arguments, "strict_f64", contract=task.role,
                                          tip_likelihoods=np.eye(4)[pattern])
                native.append(np.exp(result[0][0]))
        native = np.asarray(native)
        error = float(np.max(np.abs(native - expected)))
        assert abs(float(native.sum()) - 1) < 1e-12, task
        assert error < 1e-12, (task, error)
        with PreparedSampler(sampler_buffers(model, params), root,
                             1 if task.role == "tip_survival" else 2) as sampler:
            draws = sampler.draw(1, samples, seeds[task], 2)[0].reshape((-1, 4)).astype(int)
        counts = pattern_counts(draws, allowed)
        # Simultaneous Bernstein bound over every cell in every declared case.
        # The total predeclared type-I upper bound is 1e-6; no seed search.
        logterm = np.log(2 * len(patterns) * len(identities) / 1e-6)
        limit = np.sqrt(2 * expected * (1-expected) * logterm / samples) + 2 * logterm / (3 * samples)
        deviation = np.abs(counts / samples - expected)
        assert not np.any(counts[expected == 0]), (task, "structural zero sampled")
        assert np.all(deviation <= limit), (task, float(np.max(deviation / limit)))
        records.append(dict(task=task.to_dict(), seed=str(seeds[task]),
            task_fingerprint=task_fingerprint(20261115, task), patterns=len(patterns),
            probability_sum=float(native.sum()), max_probability_error=error,
            samples=samples, max_bernstein_fraction=float(np.max(deviation / limit))))
    return dict(native_sha256=hashlib.sha256(_resolve_native_library_path().read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                reference_script_sha256=hashlib.sha256(Path(__file__).with_name("validate_v12_gate1.py").read_bytes()).hexdigest(),
                reference="independent Q, SciPy expm, and full-tree pruning; legacy local cladogenesis operator reused",
                cases=records, family_count=len(FAMILIES), modifiers=list(MODIFIERS),
                alpha_upper_bound=1e-6, status="passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=100_000)
    args = parser.parse_args()
    report = run(args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"Passed {len(report['cases'])} model/modifier/contract cases; output {args.output}")
