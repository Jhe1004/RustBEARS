"""Tier-1a tests with analytical and independently pruned probability laws."""
from dataclasses import replace
import math
import unittest

import numpy as np

import rustbears as b
from scripts.validate_v12_gate1 import (
    Case, CORE_CASES, ReferenceTree, TREE4, assert_probability_vector,
    bin_mapping, clopper_pearson, enumerate_patterns, fault_distribution,
    goodness_statistics, holm_adjust, ideal_reference, make_model,
    monte_carlo_p, native_pattern_probabilities, pattern_counts,
    pure_extinction_check,
)


class GateStatisticsTests(unittest.TestCase):
    """Deterministic logic checks do not assert random tests always pass."""

    def test_one_holm_family_of_twelve(self):
        p = np.array([0.0001, 0.001, 0.002, *([0.5] * 9)])
        adjusted = holm_adjust(p)
        self.assertAlmostEqual(adjusted[0], 0.0012)
        self.assertAlmostEqual(adjusted[1], 0.011)
        self.assertEqual(np.count_nonzero(adjusted <= 0.01), 1)
        np.testing.assert_array_equal(holm_adjust(p[::-1]), adjusted[::-1])

    def test_mc_tail_count_includes_observation_and_ties(self):
        self.assertEqual(float(monte_carlo_p(5, np.arange(5))), 1 / 6)
        self.assertEqual(float(monte_carlo_p(0, np.arange(5))), 1)
        self.assertEqual(float(monte_carlo_p(2, np.array([2, 2, 3]))), 1)

    def test_binomial_interval_does_not_turn_99_of_100_into_guarantee(self):
        lo, hi = clopper_pearson(99, 100)
        self.assertLess(lo, 0.99)
        self.assertGreater(hi, 0.99)
        self.assertEqual(clopper_pearson(0, 100)[0], 0)
        self.assertEqual(clopper_pearson(100, 100)[1], 1)

    def test_structural_zero_is_never_hidden_in_rare_bin(self):
        p = np.array([0.0, 0.001, 0.009, 0.09, 0.9])
        bins, bp = bin_mapping(p, 1000)
        self.assertEqual(bins[0], -1)
        self.assertTrue(np.all(1000 * bp >= 10))
        invalid = goodness_statistics(np.array([1, 0, 9, 90, 900]), p, bins, bp, 1000)
        self.assertTrue(np.isinf(invalid[1]))
        valid = goodness_statistics(np.array([0, 1, 9, 90, 900]), p, bins, bp, 1000)
        np.testing.assert_allclose(valid, [0, 0], atol=1e-13)

    def test_count_order_matches_enumeration_with_noncontiguous_states(self):
        patterns = enumerate_patterns([1, 3, 5], 3)
        counts = pattern_counts(patterns, [1, 3, 5])
        np.testing.assert_array_equal(counts, np.ones(27))
        with self.assertRaises(ValueError):
            pattern_counts(np.array([[0, 1, 3]]), [1, 3, 5])

    def test_ideal_controls_are_reproducible_and_no_seed_search(self):
        p = np.array([0.02, 0.08, 0.2, 0.7])
        first = ideal_reference(p, 1000, 31, np.random.default_rng(123))
        second = ideal_reference(p, 1000, 31, np.random.default_rng(123))
        for a, c in zip(first, second):
            np.testing.assert_array_equal(a, c)


class ContractProbabilityTests(unittest.TestCase):
    def test_three_tip_extinction_analytic_counterexample(self):
        results = pure_extinction_check()
        self.assertAlmostEqual(results["legacy_unconditioned"]["sum"], math.exp(-0.6), places=14)
        self.assertAlmostEqual(results["unconditional_complete"]["sum"], 1.0, places=14)
        self.assertAlmostEqual(results["unconditional_complete"]["tip_survival_probability"],
                               math.exp(-1.8), places=14)

    def test_two_area_full_space_matches_independent_reference(self):
        for case in CORE_CASES[:3]:
            with self.subTest(case=case.name):
                model, params = make_model(case)
                reference = ReferenceTree(model, params, complete_empty=True)
                patterns = enumerate_patterns(range(len(model.states)), model.tree.tip_count)
                expected = reference.pattern_likelihoods(patterns)
                observed = native_pattern_probabilities(model, params, patterns)
                assert_probability_vector(expected)
                assert_probability_vector(observed)
                np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0)

    def test_conditional_likelihood_equals_completed_probability_over_survival(self):
        case = Case("conditional_identity", TREE4, 2, "tip_survival", 0.13, 0.2, 0.17)
        roots = {(0,): 0.1, (1,): 0.3, (0, 1): 0.6}
        model_a, params = make_model(case, root_weights=roots)
        model_b, _ = make_model(replace(case, contract="unconditional_complete"), root_weights=roots)
        weights = np.array([roots.get(tuple(s), 0.0) for s in model_a.states])
        reference = ReferenceTree(model_a, params, root_weights=weights, complete_empty=True)
        nonnull = np.flatnonzero([bool(s) for s in model_a.states])
        patterns = enumerate_patterns(nonnull, model_a.tree.tip_count)
        pa = native_pattern_probabilities(model_a, params, patterns)
        pb = native_pattern_probabilities(model_b, params, patterns)
        np.testing.assert_allclose(pa, pb / pb.sum(), atol=1e-12, rtol=0)
        expected = reference.pattern_likelihoods(patterns, conditioning="tip_survival")
        np.testing.assert_allclose(pa, expected, atol=1e-12, rtol=0)
        self.assertAlmostEqual(model_a.evaluate(params).log_tip_survival_prob,
                               math.log(float(pb.sum())), places=12)

    def test_labelled_root_order_and_probability_scale(self):
        case = Case("root_order", TREE4, 2, "unconditional_complete", 0.13, 0.2, 0.0)
        integer = {(): 0.1, (0,): 0.2, (1,): 0.3, (0, 1): 0.4}
        labels = {"XY": 4.0, "Y": 3.0, "_": 1.0, "X": 2.0}
        m1, p = make_model(case, root_weights=integer)
        m2, _ = make_model(case, root_weights=labels)
        self.assertAlmostEqual(m1.evaluate(p).lnL, m2.evaluate(p).lnL, places=13)
        patterns = enumerate_patterns(range(4), 4)
        assert_probability_vector(native_pattern_probabilities(m1, p, patterns))
        for weights in ({"unknown": 1}, {(0,): -1}, {(0,): np.nan}, {(0,): 0}):
            with self.subTest(weights=weights), self.assertRaises(ValueError):
                make_model(case, root_weights=weights)

    def test_empty_ancestor_transmits_empty_in_completed_model(self):
        case = Case("empty_root", TREE4, 2, "unconditional_complete", 0.13, 0.2, 0.15)
        model, params = make_model(case, root_weights={(): 1.0})
        batch = b.simulate_datasets(model, params=params, n_reps=101, seed=31)
        self.assertTrue(np.all(batch.tip_states == 0))
        self.assertTrue(np.all(batch.node_states == 0))
        pattern = np.zeros((1, model.tree.tip_count), dtype=np.int64)
        self.assertAlmostEqual(native_pattern_probabilities(model, params, pattern)[0], 1, places=14)

    def test_sampler_contract_defaults_and_thread_determinism(self):
        for case in (CORE_CASES[1], CORE_CASES[4]):
            with self.subTest(case=case.name):
                model, params = make_model(case)
                one = b.simulate_datasets(model, params=params, n_reps=251, seed=71, num_threads=1)
                four = b.simulate_datasets(model, params=params, n_reps=251, seed=71, num_threads=4)
                np.testing.assert_array_equal(one.tip_states, four.tip_states)
                np.testing.assert_array_equal(one.node_states, four.node_states)
                self.assertEqual(one.tip_states.shape, (251, model.tree.tip_count))
                if case.contract == "tip_survival":
                    self.assertTrue(np.all(one.tip_states != 0))
                else:
                    self.assertTrue(np.any(one.tip_states == 0))
        legacy, params = make_model(replace(CORE_CASES[0], contract="legacy_unconditioned"))
        with self.assertRaises(ValueError):
            b.simulate_datasets(legacy, params=params, n_reps=2)

    def test_prepared_simulation_batches_reproduce_the_same_replicates(self):
        model, params = make_model(CORE_CASES[4])
        with b.PreparedSimulation(model, params) as prepared:
            all_maps = prepared.draw(123, seed=91, num_threads=1)
            first = prepared.draw(43, start_rep=1, seed=91, num_threads=4)
            second = prepared.draw(80, start_rep=44, seed=91, num_threads=2)
        np.testing.assert_array_equal(all_maps.tip_states,
                                      np.concatenate((first.tip_states, second.tip_states)))
        np.testing.assert_array_equal(all_maps.node_states,
                                      np.concatenate((first.node_states, second.node_states)))

    def test_same_context_tip_replacement_is_order_independent(self):
        model, params = make_model(CORE_CASES[1])
        patterns = np.array([[1, 1, 2, 2], [0, 0, 0, 0], [1, 1, 2, 2]])
        with model.prepared_likelihood(params) as prepared:
            lnliks = prepared.evaluate_tip_states(patterns)
        self.assertEqual(lnliks[0], lnliks[2])
        self.assertTrue(np.isfinite(lnliks).all())

    def test_predeclared_faults_have_nonzero_exact_distribution_effect(self):
        combinations = (
            (CORE_CASES[1], "d_plus_two_percent"),
            (CORE_CASES[5], "omit_root_survival_weight"),
            (CORE_CASES[2], "empty_internal_restart"),
        )
        for case, fault in combinations:
            with self.subTest(fault=fault):
                model, params = make_model(case)
                allowed = np.arange(len(model.states))
                if case.contract == "tip_survival":
                    allowed = allowed[allowed != 0]
                patterns = enumerate_patterns(allowed, model.tree.tip_count)
                ref = ReferenceTree(model, params, complete_empty=case.contract == "unconditional_complete")
                correct = ref.pattern_likelihoods(patterns, conditioning="tip_survival" if case.contract == "tip_survival" else "none")
                wrong = fault_distribution(case, model, params, patterns, fault)
                assert_probability_vector(correct)
                assert_probability_vector(wrong)
                self.assertGreater(float(np.abs(correct - wrong).sum()) / 2, 1e-10)


if __name__ == "__main__":
    unittest.main()
