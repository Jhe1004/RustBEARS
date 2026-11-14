"""Compatibility guards exercise genuine fitted targets, not relabelled scores."""
from dataclasses import replace
import unittest
from unittest.mock import patch

import numpy as np

import rustbears as b
from rustbears.uncertainty import prepare_comparison_fits, validate_compatible_fits


TREE = "((T1:1,T2:1):1,T3:2);"
GEOG = b.TipRanges(("A", "B"), {"T1": (0,), "T2": (0,), "T3": (1,)})


def fitted(*, contract="legacy_unconditioned", tree=TREE, geog=GEOG,
           root=None, precision="strict_f64", d=.2, j=0):
    template = b.model_template("DEC+J" if j else "DEC").with_parameters(
        d=dict(type="fixed", init=d), e=dict(type="fixed", init=.1),
        j=dict(type="fixed", init=j))
    return b.BioGeoBEARS(tree, geog, template, likelihood_conditioning=contract,
                         root_prior=root, precision=precision,
                         allow_null_tips=contract == "unconditional_complete").fit()


class FitContractGuardTests(unittest.TestCase):
    def test_same_target_retains_fits_without_recomputation(self):
        null, alt = fitted(), fitted(j=.1)
        with patch.object(b.BioGeoBEARS, "fit", side_effect=AssertionError("unexpected refit")):
            prepared = prepare_comparison_fits(null, alt)
        self.assertIs(prepared.null_fit, null)
        self.assertIs(prepared.alt_fit, alt)
        self.assertFalse(prepared.recomputed_observed_under_new_contract)

    def test_different_input_contracts_never_enter_a_comparison(self):
        null, alt = fitted(), fitted(contract="tip_survival", j=.1)
        with self.assertRaisesRegex(b.ContractMismatchError, "contracts"):
            prepare_comparison_fits(null, alt, target_conditioning="tip_survival", refit_observed=True)

    def test_changed_target_requires_explicit_refit_and_never_relabels(self):
        null, alt = fitted(), fitted(j=.1)
        with patch.object(b.BioGeoBEARS, "fit", side_effect=AssertionError("unexpected refit")):
            with self.assertRaisesRegex(b.ContractMismatchError, "refit_observed=True"):
                prepare_comparison_fits(null, alt, target_conditioning="tip_survival")

    def test_explicit_switch_rebuilds_and_refits_both_models(self):
        null, alt = fitted(), fitted(j=.1)
        original_scores = (null.lnL, alt.lnL)
        calls = []
        real_fit = b.BioGeoBEARS.fit

        def traced_fit(model, *args, **kwargs):
            calls.append(model)
            return real_fit(model, *args, **kwargs)

        with patch.object(b.BioGeoBEARS, "fit", traced_fit):
            prepared = prepare_comparison_fits(null, alt, target_conditioning="tip_survival",
                                               refit_observed=True)
        self.assertEqual(len(calls), 2)
        self.assertTrue(prepared.recomputed_observed_under_new_contract)
        self.assertEqual(prepared.original_conditioning, "legacy_unconditioned")
        self.assertEqual(prepared.likelihood_conditioning, "tip_survival")
        for old, new in zip((null, alt), prepared.fits):
            self.assertIsNot(old.model, new.model)
            self.assertIsNot(old.model.tree, new.model.tree)
            self.assertIsNot(old.model.geog, new.model.geog)
            self.assertEqual(new.likelihood_conditioning, "tip_survival")
            self.assertAlmostEqual(new.lnL, new.model.evaluate(new.params).lnL, places=13)
            self.assertGreater(abs(new.lnL - old.lnL), .1)
        self.assertEqual((null.lnL, alt.lnL), original_scores)
        self.assertEqual(null.model.likelihood_conditioning, "legacy_unconditioned")

    def test_custom_root_measure_survives_conversion(self):
        prior = b.RootPriorSpec("custom_weights", {"A": 2, "B": 1, "AB": .5})
        null, alt = fitted(root=prior), fitted(root=prior, j=.1)
        prepared = prepare_comparison_fits(null, alt, target_conditioning="unconditional_complete",
                                           refit_observed=True)
        np.testing.assert_array_equal(null.model.root_likelihood_weights, [0, 2, 1, .5])
        for fit in prepared.fits:
            np.testing.assert_array_equal(fit.model.root_likelihood_weights, np.array([0, 2, 1, .5])/3.5)
            self.assertTrue(fit.model.config.allow_null_tips)

    def test_rotated_identical_tree_is_accepted(self):
        null = fitted()
        alt = fitted(tree="(T3:2,(T2:1,T1:1):1);", j=.1)
        self.assertEqual(validate_compatible_fits([null, alt]), (null, alt))

    def test_tree_geography_precision_and_root_scale_must_match(self):
        null = fitted()
        variants = (
            (fitted(tree="((T1:1,T2:1):2,T3:3);"), "trees"),
            (fitted(geog=b.TipRanges(("A", "B"), {"T1": (1,), "T2": (0,), "T3": (1,)})), "geography"),
            (fitted(precision="legacy_f32"), "precision"),
            (fitted(root=b.RootPriorSpec("custom_weights", [0, 2, 2, 2])), "weights"),
            (fitted(root=b.RootPriorSpec("custom_weights", [0, 1, 2, 3])), "weights"),
        )
        for alt, label in variants:
            with self.subTest(label=label), self.assertRaisesRegex(b.ContractMismatchError, label):
                validate_compatible_fits((null, alt))

    def test_empty_observation_cannot_be_switched_to_survival(self):
        geog = b.TipRanges(("A", "B"), {"T1": (), "T2": (0,), "T3": (1,)})
        null = fitted(contract="unconditional_complete", geog=geog)
        alt = fitted(contract="unconditional_complete", geog=geog, j=.1)
        with patch.object(b.BioGeoBEARS, "fit", side_effect=AssertionError("unexpected partial refit")):
            with self.assertRaisesRegex(ValueError, "null tip"):
                prepare_comparison_fits(null, alt, target_conditioning="tip_survival", refit_observed=True)

    def test_invalid_fit_and_mutated_metadata_fail(self):
        fit = fitted()
        with self.assertRaisesRegex(ValueError, "converged"):
            validate_compatible_fits((replace(fit, converged=False),))
        with self.assertRaisesRegex(ValueError, "finite"):
            validate_compatible_fits((replace(fit, lnL=np.nan),))
        bad = replace(fit, evaluation=replace(fit.evaluation, precision_mode="legacy_f32"))
        with self.assertRaisesRegex(b.ContractMismatchError, "metadata"):
            validate_compatible_fits((bad,))
        with self.assertRaises(TypeError):
            prepare_comparison_fits(fit, fit, refit_observed="yes")

    def test_empirical_posterior_root_mode_is_explicitly_disabled(self):
        fit = fitted(contract="tip_survival")
        with self.assertRaisesRegex(b.UnsupportedModelContractError, "already-conditioned"):
            b.PreparedSimulation(fit, root_distribution="empirical_posterior")


if __name__ == "__main__":
    unittest.main()
