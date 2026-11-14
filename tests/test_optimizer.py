from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np

from rustbears import BioGeoBEARS
from rustbears._optimization import SearchResult, optimize_likelihood
from rustbears._refinement import refine_completed_search


class OptimizerRegression(unittest.TestCase):
    def search(self, fn, start=(0.01, 0.01), height=150):
        return optimize_likelihood(
            fn, ["d", "e"], start, [(1e-12, 5), (1e-12, 5)], height, "legacy_f32", "L-BFGS-B", 300, 1
        )

    def test_deep_tree_numerical_wall_has_recoverable_direction(self):
        def loglike(theta, precision):
            d, e = theta
            if e > 1:
                raise ValueError("invalid action input")
            return -((d - 0.0003) / 0.0002) ** 2 - ((e - 0.0001) / 0.0002) ** 2

        fitted = self.search(loglike, start=(0.01, 5), height=300)
        self.assertGreater(loglike(fitted.theta, "legacy_f32"), -1e-6)

    def test_escape_lower_bound_false_stationarity(self):
        def loglike(theta, precision):
            return -((theta[0] - 0.004) / 0.002) ** 2 - ((theta[1] - 0.0002) / 0.001) ** 2

        fitted = self.search(loglike, start=(1e-12, 0.01))
        self.assertAlmostEqual(fitted.theta[0], 0.004, places=7)
        self.assertGreater(loglike(fitted.theta, "legacy_f32"), -1e-6)

    def test_zero_extinction_competes_with_positive_extinction_basin(self):
        def loglike(theta, precision):
            d, e = theta
            return max(
                1 - ((d - 0.00037) / 0.00015) ** 2 - ((e - 4e-7) / 0.00015) ** 2,
                -((d - 0.002) / 0.001) ** 2 - ((e - 0.006) / 0.002) ** 2,
            )

        fitted = self.search(loglike, height=173)
        self.assertGreater(loglike(fitted.theta, "legacy_f32"), 0.9999)
        self.assertLess(fitted.theta[1], 1e-5)

    def test_programming_error_is_not_converted_to_numeric_penalty(self):
        def loglike(theta, precision):
            raise ValueError("unrelated shape error")

        with self.assertRaisesRegex(ValueError, "unrelated shape error"):
            self.search(loglike)

    def test_returned_score_is_target_precision_and_probe_is_reference_free(self):
        def loglike(theta, precision):
            center = 0.002 if precision == "strict_f64" else 0.0021
            return -float(np.sum((theta - center) ** 2)) * 1e5

        fitted = self.search(loglike)
        self.assertGreater(loglike(fitted.theta, "legacy_f32"), -1e-7)
        self.assertEqual(len(fitted.proposal["probes"]), 4)
        self.assertIn("probe_scale", fitted.proposal)

    def test_absolute_difference_resolves_tiny_boundary_gradient(self):
        def ll(theta, precision):
            return -1e5 - ((theta[0] - 0.005) / 0.01) ** 2 - 1000 * (theta[1] - 0.002) ** 2

        theta = np.array([0.005, 1e-12])
        base = SearchResult(
            theta,
            1,
            True,
            "boundary plateau",
            dict(theta=theta.copy(), lnL=ll(theta, "strict_f64"), message="own search", traces=[]),
        )
        out = refine_completed_search(ll, ["d", "e"], [0.01, 0.01], [(1e-12, 5)] * 2, "legacy_f32", 300, base)
        self.assertGreater(out.theta[1], 0.0019)
        self.assertLess(abs(ll(out.theta, "legacy_f32") + 1e5), 1e-8)

    def test_completed_phase_resume_matches_public_composite_search(self):
        captured = []

        def capture(*args):
            value = optimize_likelihood(*args)
            captured.append(deepcopy(value))
            return value

        model = BioGeoBEARS("((A:3,B:3):3,C:6);", {"A": "10", "B": "10", "C": "01"}, model="DEC")
        with patch("rustbears._optimization.optimize_likelihood", side_effect=capture):
            full = model.fit()
        free = list(model.model.free)
        start = model.model.resolve()
        x0 = [start[k] for k in free]
        bounds = [(model.model.params[k].min, model.model.params[k].max) for k in free]
        with model._context(1) as ctx:
            def ll(theta, precision):
                return model._evaluate(ctx, model.model.resolve(dict(zip(free, theta))), False, precision)

            resumed = refine_completed_search(ll, free, x0, bounds, "legacy_f32", 500, captured[0])
        np.testing.assert_array_equal(resumed.theta, [full.params[k] for k in free])
        self.assertGreaterEqual(resumed.evaluations, captured[0].evaluations)
        self.assertGreaterEqual(full.n_evals, captured[0].evaluations)
        self.assertAlmostEqual(model.evaluate(dict(zip(free, resumed.theta))).lnL, full.lnL, places=12)

    def test_public_fit_on_tall_tree_improves_default_feasible_point(self):
        model = BioGeoBEARS("((A:150,B:150):150,C:300);", {"A": "100", "B": "010", "C": "001"}, model="DEC+J")
        initial = {"d": 0.01, "e": 0.01, "j": 0.0001}
        baseline = model.evaluate(initial).lnL
        fit = model.fit(initial={"j": 0.0001}, num_threads=2)
        self.assertGreater(fit.lnL, baseline + 1e-3)
        self.assertTrue(np.isfinite(fit.evaluation.ancestral_probs_branch_top).all())
        self.assertAlmostEqual(fit.lnL, model.evaluate(fit.params).lnL, places=12)


if __name__ == "__main__":
    unittest.main()
