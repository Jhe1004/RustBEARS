#!/usr/bin/env python3
"""Link completed Gate 1 evidence to an updated source fingerprint.

Rechecks every exact probability and every native simulated pattern count with
the SAME seeds and settings.  Previously recorded multinomial reference draws
are reused only after the probability vector is independently revalidated.
The original report is never edited or promoted silently.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

import rustbears as b
from scripts.validate_v12_gate1 import (
    CORE_CASES, ReferenceTree, enumerate_patterns, goodness_statistics,
    holm_adjust, make_model, monte_carlo_p, native_pattern_probabilities,
    pattern_counts, pure_extinction_check, sha256, write_json, bin_mapping,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--allow-native-change", action="store_true",
                        help="Revalidate a new native build against every original probability and count")
    args = parser.parse_args()
    original = Path(args.original_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    original_report = original / "gate1_results.json"
    prior = json.loads(original_report.read_text())
    settings = prior["settings"]
    native = b._ffi._resolve_native_library_path()
    if sha256(native) != prior["native_sha256"] and not args.allow_native_change:
        raise RuntimeError("native library changed; this source-only association cannot be used")
    source_names = set(prior["source_sha256"])
    source_names.update(str(path.relative_to(REPOSITORY)) for path in
                        (REPOSITORY / "rustbears").rglob("*.py"))
    source_names.update(str(path.relative_to(REPOSITORY)) for path in
                        (REPOSITORY / "tests").glob("test_v12*.py"))
    source_hashes = {name: sha256(REPOSITORY / name) for name in sorted(source_names)}
    report = {
        "status": "running", "original_report": str(original_report),
        "original_report_sha256": sha256(original_report), "source_sha256": source_hashes,
        "source_changes_since_original": [name for name, digest in source_hashes.items()
                                           if digest != prior["source_sha256"].get(name)],
        "revalidation_script_sha256": sha256(__file__), "native_sha256": sha256(native),
        "original_native_sha256": prior["native_sha256"],
        "native_changed_since_original": sha256(native) != prior["native_sha256"],
        "native_change_revalidation_explicitly_enabled": args.allow_native_change,
        "protocol": "Same-seed native probability and count replay; cached null draws only after exact-law validation.",
        "settings": settings, "cases": [],
    }
    destination = output / "source_revalidation.json"
    started = time.monotonic()
    write_json(destination, report)
    try:
        report["pure_extinction"] = pure_extinction_check()
        pvalues = []
        for index, case in enumerate(CORE_CASES):
            artifact = original / f"{case.name}.npz"
            cached = np.load(artifact)
            model, params = make_model(case)
            patterns = enumerate_patterns(cached["allowed_states"], model.tree.tip_count)
            reference = ReferenceTree(model, params, complete_empty=case.contract == "unconditional_complete")
            probability = reference.pattern_likelihoods(
                patterns, conditioning="tip_survival" if case.contract == "tip_survival" else "none")
            native_p = native_pattern_probabilities(model, params, patterns)
            reference_error = float(np.max(np.abs(probability / probability.sum() - cached["probabilities"])))
            native_identical = bool(np.array_equal(native_p, cached["native_probabilities"]))
            if reference_error > 1e-15 or not native_identical:
                raise AssertionError(f"{case.name}: exact probabilities differ from original evidence")
            batch = b.simulate_datasets(model, params=params, n_reps=settings["samples"],
                                       seed=settings["seed"] + index, num_threads=settings["threads"])
            counts = pattern_counts(batch.tip_states, cached["allowed_states"])
            if not np.array_equal(counts, cached["counts"]):
                raise AssertionError(f"{case.name}: same-seed native counts changed")
            bins, bin_p = bin_mapping(cached["probabilities"], settings["samples"])
            observed = goodness_statistics(counts, cached["probabilities"], bins, bin_p, settings["samples"])
            pvalues.extend(monte_carlo_p(observed, cached["ideal_statistics"]))
            row = {"configuration": case.name, "artifact_sha256": sha256(artifact),
                   "reference_max_error": reference_error,
                   "native_probabilities_bit_identical": native_identical,
                   "native_pattern_counts_identical": True, "samples": int(counts.sum())}
            report["cases"].append(row)
            write_json(destination, report)
            print(json.dumps(row), flush=True)
        adjusted = holm_adjust(pvalues)
        report["holm_adjusted_pvalues"] = adjusted.tolist()
        report["simulator_rejected"] = bool(np.any(adjusted <= 0.01))
        report["source_files_unchanged_during_revalidation"] = all(
            sha256(REPOSITORY / name) == digest for name, digest in source_hashes.items())
        report["native_unchanged_during_revalidation"] = sha256(native) == report["native_sha256"]
        if not report["source_files_unchanged_during_revalidation"] or not report["native_unchanged_during_revalidation"]:
            raise RuntimeError("source/native changed during source revalidation")
        report["fault_reference_evidence"] = {
            "reused_from_original": True,
            "reason": "Exact null probability laws and harness are unchanged; original mutant trials were not selectively rerun.",
            "results": prior["fault_sensitivity"],
        }
        report["status"] = "verified_same_evidence_under_current_source"
        report["seconds"] = time.monotonic() - started
        write_json(destination, report)
    except Exception as error:
        report["status"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        report["seconds"] = time.monotonic() - started
        write_json(destination, report)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
