import gc
from pathlib import Path
import tempfile
import unittest

import numpy as np
from pandas.testing import assert_frame_equal

import rustbears as b

TREE = "((A:1,B:1):1,C:2);"
GEOG = {"A": "10", "B": "10", "C": "01"}


class PackageTests(unittest.TestCase):
    def test_alias_module_parity(self):
        self.assertEqual(b.__version__, "1.1.4")
        self.assertIs(b.RustBEARS, b.BioGeoBEARS)

    def test_tree_order_and_fossil(self):
        tree = b.read_newick("(((A:1,'B b':1):1,C:2):1,(D:1,E:0.8):2)base:0.25;")
        self.assertEqual(tree.root_edge, 0.25)
        self.assertEqual(tree.tip_labels, ("A", "B b", "C", "D", "E"))
        self.assertEqual(tree.labels[tree.root], "base")
        np.testing.assert_array_equal(
            tree.edges + 1,
            [[8, 1], [8, 2], [9, 4], [9, 5], [7, 8], [7, 3], [6, 7], [6, 9]],
        )
        self.assertTrue(tree.node_table().loc[4, "fossils"])
        self.assertEqual(tree.node_table().loc[0, "parent_br"], 3)

    def test_areas_and_null(self):
        with self.assertRaisesRegex(ValueError, "null tip range"):
            b.BioGeoBEARS(TREE, {**GEOG, "A": "00"})
        m = b.BioGeoBEARS(TREE, {**GEOG, "A": "00"}, allow_null_tips=True)
        self.assertEqual(m.tips[0, 0], 1)
        with self.assertRaisesRegex(ValueError, "taxon sets differ"):
            b.BioGeoBEARS(TREE, {**GEOG, "D": "10"})

    def test_layer_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "test-layers.txt"
            p.write_text("# two layers\n1 0.5\n0.5 NA\nEND\n1 1\n1 1\n")
            from rustbears.io import read_layer_matrices

            a = read_layer_matrices(p, 2, 2, missing_diagonal=True)
            self.assertEqual(a[0][1, 1], 1)
            self.assertEqual(a[1][0, 1], 1)

    def test_single_and_stratified_equal_if_layers_identical(self):
        for p in ("legacy_f32", "strict_f64"):
            a = b.BioGeoBEARS(TREE, GEOG, precision=p).evaluate()
            z = b.BioGeoBEARS(TREE, GEOG, timeperiods=[0.5, 1.5, 10], precision=p).evaluate()
            self.assertLess(abs(a.lnL - z.lnL), 1e-7 if p == "legacy_f32" else 1e-12)
            np.testing.assert_allclose(
                a.ancestral_probs_branch_top,
                z.ancestral_probs_branch_top,
                atol=1e-7,
                rtol=0,
            )
        with self.assertRaisesRegex(ValueError, "too close"):
            b.section_the_tree(TREE, [1, 10])

    def test_b_free_reuses_context(self):
        template = b.model_template("DEC").with_parameters(b=dict(type="free"))
        m = b.BioGeoBEARS(TREE, GEOG, template)
        before = b.context_stats()
        with m._context(1) as context:
            first = m._evaluate(context, m.model.resolve(), True)
            second = m._evaluate(context, m.model.resolve({"b": 0.5}), True)
            third = m._evaluate(context, m.model.resolve(), True)
        after = b.context_stats()
        self.assertNotEqual(first.lnL, second.lnL)
        self.assertEqual(first.lnL, third.lnL)
        self.assertEqual(after["created"] - before["created"], 1)
        self.assertEqual(after["dropped"] - before["dropped"], 1)

    def test_parameter_formula_is_not_executable(self):
        with self.assertRaisesRegex(ValueError, "only names, numbers and arithmetic"):
            b.model_template("DEC").with_parameters(ysv={"type": "print(3)"})

    def test_context_lifecycle(self):
        m = b.BioGeoBEARS(TREE, GEOG)
        before = b.context_stats()
        for _ in range(4):
            with m._context(1) as context:
                v = m._evaluate(context, m.model.resolve(), True)
                self.assertTrue(np.isfinite(v.lnL))
            context.close()
        context = m._context(1)
        del context
        gc.collect()
        after = b.context_stats()
        self.assertEqual(after["created"] - before["created"], 5)
        self.assertEqual(after["dropped"] - before["dropped"], 5)

    def test_native_bsm_threads_and_rejection(self):
        model = b.BioGeoBEARS(TREE, GEOG, timeperiods=[0.5, 1.5, 10])
        p = model.model.resolve({"d": 0.2, "e": 0.1})
        ev = model.evaluate(p)
        fit = b.FitResult(ev.lnL, p, 2, 0, 0, 0, True, ev, model)
        base = fit.run_bsm(100, num_threads=1)
        for n in (4, 8):
            other = fit.run_bsm(100, mode="parallel_substream", num_threads=n)
            assert_frame_equal(base.clado_events_df, other.clado_events_df, check_exact=True)
            assert_frame_equal(base.ana_events_df, other.ana_events_df, check_exact=True)
            for k, v in base.summary_counts.items():
                self.assertTrue(v.equals(other.summary_counts[k]))
        self.assertEqual(len(base.clado_events_tables), 100)
        self.assertEqual(len(base.ana_events_tables), 100)
        self.assertEqual(len(fit.run_bsm(10, mode="rejection").clado_events_df), 20)

    def test_traits_and_fixed_likelihood(self):
        trait = b.TraitModel(("a", "b"), [1, 2], [[0, 0.04], [0.02, 0]])
        for times in (None, [0.5, 1.5, 10]):
            m = b.BioGeoBEARS(
                TREE,
                GEOG,
                timeperiods=times,
                trait_values={"A": "a", "B": "a", "C": "b"},
                trait_model=trait,
                fixnode=5,
                fixlikes=np.ones(8),
            )
            v = m.evaluate()
            self.assertTrue(np.isfinite(v.lnL))
            self.assertEqual(v.ancestral_probs_branch_top.shape, (5, 8))
            p = m.model.resolve()
            f = b.FitResult(v.lnL, p, 2, 0, 0, 0, True, v, m)
            self.assertEqual(len(f.run_bsm(5).clado_events_df), 10)

    def test_statistic_contract(self):
        self.assertAlmostEqual(b.lrttest(-2, -4, 3, 2), 0.04550026389635857)
        self.assertEqual(b.lrttest(-2, -2, 2, 2), 1)
        self.assertEqual(b.lrttest(-2, -3, 2, 2), 0)
        model = b.BioGeoBEARS(TREE, GEOG)
        result = model.fit()
        self.assertTrue(result.converged)
        self.assertEqual(result.aic, 2 * result.n_params - 2 * result.lnL)
        self.assertEqual(b.compare_models([result]).loc[0, "AICweight"], 1)


if __name__ == "__main__":
    unittest.main()
