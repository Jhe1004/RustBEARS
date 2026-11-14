import copy
import io
from pathlib import Path
import pickle
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import rustbears as b
from rustbears.__main__ import _run_cli
from rustbears._ffi import kernel

TREE = "((A:1,B:1):1,C:2);"
GEOG = {"A": "10", "B": "10", "C": "01"}


class AuditRegressionTests(unittest.TestCase):
    def test_compiled_context_forbids_copy_and_pickle(self):
        model = b.RustBEARS(TREE, GEOG)
        with model._context(1) as ctx:
            with self.assertRaises(TypeError):
                copy.copy(ctx)
            with self.assertRaises(TypeError):
                copy.deepcopy(ctx)
            with self.assertRaises(TypeError):
                pickle.dumps(ctx)
            res = model._evaluate(ctx, model.model.resolve(), False)
            self.assertTrue(np.isfinite(res))

    def test_rejection_bsm_refuses_truncated_high_rate_branch(self):
        long_tree = "((A:4000,B:4000):4000,C:8000);"
        spec = b.model_template("DEC").with_parameters(
            d=b.Parameter("fixed", 0.0, 0.0, 5.0),
            e=b.Parameter("fixed", 0.0, 0.0, 5.0),
            a=b.Parameter("fixed", 4.5, 0.0, 5.0),
        )
        fit = b.RustBEARS(long_tree, GEOG, model=spec, include_null_range=False).fit()
        with self.assertRaisesRegex(ValueError, "exceeded 10000 transitions"):
            fit.run_bsm(num_maps=2, mode="rejection", seed=7)

    def test_coo_kernels_reject_negative_nan_and_fractional_indices(self):
        for bad in (-1.0, -0.25, 0.5, float("nan"), float("inf"), 3.0):
            with self.assertRaises(ValueError):
                kernel(6, [bad], [1.0], [3.0])
            with self.assertRaises(ValueError):
                kernel(7, [1.0, 1.0], [1.0, 1.0], [bad], [0.0], [0.0], [1.0], [1.0, 1.0, 1.0])
            with self.assertRaises(ValueError):
                kernel(8, [1.0, 1.0], [1.0, 1.0], [bad], [0.0], [0.0], [1.0], [1.0, 1.0, 1.0])

    def test_fractional_and_boolean_discrete_indices_rejected(self):
        states = b.rcpp_areas_list_to_states_list(["A", "B"])
        q = b.rcpp_states_list_to_DEmat(["A", "B"], states, np.ones((2, 2)), [0.01, 0.01])
        self.assertEqual(q.shape, (4, 4))

        for bad_state in ((0.9,), (-0.1,), (-2,), (True,), (-1, 0)):
            with self.assertRaises(ValueError):
                b.rcpp_states_list_to_DEmat(["A", "B"], [(), bad_state, (1,)], np.ones((2, 2)), [0.01, 0.01])

        for bad_fix in (4.5, 3.9, float("nan"), True):
            with self.assertRaises(ValueError):
                b.RustBEARS(TREE, GEOG, fixnode=bad_fix, fixlikes=np.ones((1, 4)))

        model = b.RustBEARS(TREE, GEOG)
        for bad_threads in (0, -1, 2**64, 100_000, True, False):
            with self.assertRaises(ValueError):
                model.evaluate(num_threads=bad_threads)

        fit = model.fit()
        for bad_bsm_arg in (True, False, 0, -2, 2**64):
            with self.assertRaises(ValueError):
                fit.run_bsm(num_maps=2, num_threads=bad_bsm_arg)

    def test_dataframe_geography_rejects_non_binary_values(self):
        valid = pd.DataFrame({"X": [1, 0.0, "1"], "Y": [0, 1.0, "0"]}, index=["A", "B", "C"])
        parsed = b.read_geog(valid)
        self.assertEqual(parsed.ranges["A"], (0,))
        self.assertEqual(parsed.ranges["B"], (1,))

        for bad in (0.9, 1.9, -0.1, 2, True, float("nan"), "2", "-"):
            df = pd.DataFrame({"X": [bad, 0], "Y": [0, 1]}, index=["A", "B"])
            with self.assertRaises(ValueError):
                b.read_geog(df)

    def test_tip_ambiguities_flag_and_likelihood_expansion(self):
        ambig_geog = {"A": "1?", "B": "01", "C": "??"}
        with self.assertRaisesRegex(ValueError, "use_ambiguities=True"):
            b.RustBEARS(TREE, ambig_geog)

        m = b.RustBEARS(TREE, ambig_geog, use_ambiguities=True)
        # States order with null range: (), (0,), (1,), (0, 1)
        np.testing.assert_array_equal(m.tips[0], [0.0, 1.0, 0.0, 1.0])
        np.testing.assert_array_equal(m.tips[1], [0.0, 0.0, 1.0, 0.0])
        np.testing.assert_array_equal(m.tips[2], [1.0, 1.0, 1.0, 1.0])
        self.assertTrue(np.isfinite(m.evaluate().lnL))

        df_ambig = pd.DataFrame({"X": ["1", "0", "?"], "Y": ["?", "1", "?"]}, index=["A", "B", "C"])
        m_df = b.RustBEARS(TREE, df_ambig, use_ambiguities=True)
        self.assertAlmostEqual(m.evaluate().lnL, m_df.evaluate().lnL, places=14)

    def test_deep_pectinate_newick_and_stratified_sectioning_without_recursion_error(self):
        n_tips = 1300
        newick = "T0:1"
        for k in range(1, n_tips):
            newick = f"({newick},T{k}:1):1"
        newick += ";"
        tree = b.read_newick(newick)
        self.assertEqual(tree.tip_count, n_tips)
        self.assertEqual(tree.edges.shape, (2 * n_tips - 2, 2))
        table = tree.node_table()
        self.assertEqual(len(table), 2 * n_tips - 1)
        sections = b.section_the_tree(tree, [250.5, 650.5, 1050.5, 2000.0])
        self.assertEqual(len(sections.pieces), 4)

    def test_newick_file_with_parentheses_and_overflow_branch_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree_file = Path(tmp) / "tree (replicate 1).nwk"
            tree_file.write_text("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);\n")
            parsed = b.read_newick(str(tree_file))
            self.assertEqual(parsed.tip_count, 5)

        with self.assertRaisesRegex(ValueError, "overflowed|finite"):
            b.read_newick("((A:1e308,B:1):1e308,C:1);")

    def test_cli_formats_expected_errors_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            stderr_buf = io.StringIO()
            with patch("sys.stderr", stderr_buf):
                with self.assertRaises(SystemExit) as cm:
                    _run_cli(
                        [
                            "fit",
                            "--tree",
                            str(Path(tmp) / "missing (1).nwk"),
                            "--geog",
                            str(Path(tmp) / "missing.data"),
                            "--out-dir",
                            str(Path(tmp) / "out"),
                        ]
                    )
            self.assertEqual(cm.exception.code, 1)
            err_text = stderr_buf.getvalue()
            self.assertTrue(err_text.startswith("rustbears: error:"), err_text)
            self.assertNotIn("Traceback", err_text)

    def test_fit_handles_equal_bounds_on_free_parameters(self):
        tree5 = "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);"
        geog5 = {"A": "10", "B": "10", "C": "01", "D": "01", "E": "11"}

        all_equal = b.model_template("DEC").with_parameters(
            d=b.Parameter("free", 0.02, 0.02, 0.02),
            e=b.Parameter("free", 0.01, 0.01, 0.01),
        )
        all_fixed = b.model_template("DEC").with_parameters(
            d=b.Parameter("fixed", 0.02, 0.0, 5.0),
            e=b.Parameter("fixed", 0.01, 0.0, 5.0),
        )
        fit_eq = b.RustBEARS(tree5, geog5, model=all_equal).fit()
        fit_fx = b.RustBEARS(tree5, geog5, model=all_fixed).fit()
        self.assertEqual(fit_eq.n_params, 0)
        self.assertEqual(fit_fx.n_params, 0)
        self.assertAlmostEqual(fit_eq.lnL, fit_fx.lnL, places=12)
        self.assertAlmostEqual(fit_eq.aic, fit_fx.aic, places=12)
        self.assertAlmostEqual(fit_eq.aicc, fit_fx.aicc, places=12)

        mixed_equal = b.model_template("DEC").with_parameters(
            d=b.Parameter("free", 0.02, 0.02, 0.02),
            e=b.Parameter("free", 0.01, 1e-12, 1.0),
        )
        mixed_fixed = b.model_template("DEC").with_parameters(
            d=b.Parameter("fixed", 0.02, 0.0, 5.0),
            e=b.Parameter("free", 0.01, 1e-12, 1.0),
        )
        fit_mix_ref = b.RustBEARS(tree5, geog5, model=mixed_fixed).fit()
        self.assertEqual(fit_mix_ref.n_params, 1)
        for opt in ("L-BFGS-B", "differential_evolution", "dual_annealing"):
            fit_mix = b.RustBEARS(tree5, geog5, model=mixed_equal).fit(optimizer=opt, maxiter=15)
            self.assertEqual(fit_mix.n_params, 1)
            self.assertTrue(np.isfinite(fit_mix.lnL))
            self.assertAlmostEqual(fit_mix.params["d"], 0.02)
        fit_mix_lbfgsb = b.RustBEARS(tree5, geog5, model=mixed_equal).fit()
        self.assertAlmostEqual(fit_mix_lbfgsb.lnL, fit_mix_ref.lnL, places=10)
        self.assertAlmostEqual(fit_mix_lbfgsb.aicc, fit_mix_ref.aicc, places=10)

    def test_plotting_closes_figures_on_save_error(self):
        fit = b.RustBEARS(TREE, GEOG).fit()
        bsm = fit.run_bsm(num_maps=4, seed=3)
        before = set(plt.get_fignums())
        with self.assertRaises(OSError):
            fit.plot("/nonexistent_rustbears_dir_9999/ancestral.png")
        with self.assertRaises(OSError):
            bsm.plot("/nonexistent_rustbears_dir_9999/bsm.png")
        self.assertEqual(set(plt.get_fignums()), before)


if __name__ == "__main__":
    unittest.main()
