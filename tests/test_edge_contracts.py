import csv
import json
import math
from pathlib import Path
import unittest

import numpy as np

import rustbears as b
from rustbears.model import BioGeoBEARSModelObject, Parameter, model_template

FIXTURES = Path(__file__).resolve().parent / "fixtures"
TREE = "((T1:1,T2:1):1,T3:2);"
GEOG = {"T1": "100", "T2": "010", "T3": "001"}


def _number(x):
    try:
        n = float(x)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def _records():
    return json.loads((FIXTURES / "empirical_edge_cases.json").read_text())


def _state_names(m):
    return ["|".join(m.areas[i] for i in s) if s else "_" for s in m.states]


def _make_model(r, precision="legacy_f32"):
    opts = {k: FIXTURES / v for k, v in r["files"].items()}
    for k in ("include_null_range", "min_branchlength", "force_sparse"):
        if r["run_settings"].get(k) is not None:
            opts[k] = r["run_settings"][k]
    template = model_template(r["model"])
    supplied = r.get("initial_parameters")
    if supplied:
        p = {
            x["_row"]: Parameter(x["type"], float(x["init"]), float(x["min"]), float(x["max"]))
            for x in supplied
        }
        for name in template.free:
            x = p.get(name)
            if x and x.type == "free":
                p[name] = Parameter(x.type, template.params[name].init, x.min, x.max)
        template = BioGeoBEARSModelObject(template.name, p)
    return b.BioGeoBEARS(model=template, max_range_size=r["max_range_size"], precision=precision, **opts)


def _frozen_params(r):
    return {
        k: _number(r["reference"].get(k))
        for k in ("d", "e", "j")
        if _number(r["reference"].get(k)) is not None
    }


def _root_prob_max_error(m, ev, r):
    rel_file = r["reference"].get("root_probability_file")
    names = {name: i for i, name in enumerate(_state_names(m))}
    p = ev.ancestral_probs_branch_top
    peak = 0.0
    with (FIXTURES / rel_file).open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            ref = _number(row["probability"])
            if ref is None:
                continue
            node = int(row["node"]) - 1 if "node" in row else m.tree.root
            idx = names.get(row["range"])
            v = p[node, idx] if idx is not None else 0.0
            peak = max(peak, abs(v - ref))
    return peak


class EmpiricalEdgeContractRegression(unittest.TestCase):
    def test_full_spatial_submatrix_and_fixed_likelihood_width(self):
        allowed = np.ones((3, 3))
        allowed[1, 2] = allowed[2, 1] = 0
        for k in ("areas_allowed", "areas_adjacency"):
            m = b.BioGeoBEARS(TREE, GEOG, **{k: allowed}, fixnode=4, fixlikes=np.ones((1, 8)))
            expected = 7 if k == "areas_allowed" else 6
            self.assertEqual(len(m.states), expected)
            self.assertEqual((0, 1, 2) in m.states, k == "areas_allowed")
            self.assertNotIn((1, 2), m.states)
            self.assertEqual(m.fixlikes.shape, (1, expected))
            self.assertTrue(np.isfinite(m.evaluate().lnL))

    def test_times_extension_and_trim(self):
        args = (TREE, GEOG)
        plain = b.BioGeoBEARS(*args).evaluate().lnL
        for times in ([], [10], [10, 20, 30], [2, 10]):
            self.assertAlmostEqual(b.BioGeoBEARS(*args, timeperiods=times).evaluate().lnL, plain)
        a = np.stack([np.ones((3, 3)) * (i + 1) for i in range(4)])
        m = b.BioGeoBEARS(*args, timeperiods=[0.5, 3, 6, 9], dispersal_multipliers=a)
        z = b.BioGeoBEARS(*args, timeperiods=[0.5, 1.5], dispersal_multipliers=a[:2])
        self.assertEqual(m.times, (0.5, 3.0))
        self.assertEqual(z.times, (0.5, 3.0))
        self.assertAlmostEqual(m.evaluate().lnL, z.evaluate().lnL, places=13)

    def _real(self, ident, model):
        return [r for r in _records() if r["screen"] == ident and r["model"] == model]

    def _audit(self, r):
        m = _make_model(r, "legacy_f32")
        x = m.evaluate(_frozen_params(r))
        self.assertTrue(np.isfinite(x.ancestral_probs_branch_top).all())
        self.assertTrue(np.isfinite(np.delete(x.ancestral_probs_branch_bottom, m.tree.root, axis=0)).all())
        self.assertLessEqual(abs(x.lnL - float(r["reference"]["LnL"])), 1e-12)
        self.assertLessEqual(_root_prob_max_error(m, x, r), 1e-14)
        return m

    def test_bsm169_9_states(self):
        for model in ("DIVALIKE", "DIVALIKE+J"):
            (r,) = self._real("BSM169", model)
            self.assertEqual(len(self._audit(r).states), 9)

    def test_bsm001_oldest_rules_extend_to_root(self):
        rows = self._real("BSM001", "DEC+J")
        self.assertGreaterEqual(len(rows), 2)
        for r in rows:
            with self.subTest(config=r["config"]):
                m = self._audit(r)
                self.assertEqual(m.raw_times, (1.0, 20.0, 60.0))
                self.assertGreater(m.times[-1], m.tree.ages[m.tree.root])

    def test_bsm085_masked_uppass(self):
        (r,) = self._real("BSM085", "BAYAREALIKE+J")
        m = self._audit(r)
        self.assertEqual(len(m.states), 99)


if __name__ == "__main__":
    unittest.main()
