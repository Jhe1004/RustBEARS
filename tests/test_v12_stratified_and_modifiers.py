"""Comprehensive mathematical, statistical and workflow tests for Tier-1/Tier-2
spatial modifiers and time-stratified monotone island-emergence probability contracts."""
from __future__ import annotations

import itertools
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.linalg import expm

import rustbears as b
from rustbears._ffi import kernel
from rustbears.operators import _states


class IndependentStratifiedReference:
    """First-principles multi-stratum continuous-time Markov + cladogenesis reference.

    Uses plain Python rate assembly and SciPy ``expm`` along exact branch slices,
    without touching ``CompiledContext``, ``CompiledStratifiedContext``, ``uppass``
    or ``PreparedSampler``.  Marginal ancestral probabilities at branch tops and
    bottoms are computed directly from the definition of conditional probability
    ``P(X_node = s, Y = y) / P(Y = y)`` via indicator projection.
    """

    def __init__(self, model: b.BioGeoBEARS, params: dict, *, complete_empty: bool = True):
        self.model = model
        self.states = tuple(tuple(s) for s in model.states)
        self.n = len(self.states)
        self.nt = model.tree.tip_count
        self.nn = len(model.tree.labels)
        self.root = int(model.tree.root)
        self.ages = np.asarray(model.tree.ages, dtype=float)
        self.parent = np.asarray(model.tree.parent, dtype=int)
        self.times = tuple(float(t) for t in model.times)
        self.masks = tuple(np.asarray(m, dtype=bool) for m in model.masks)
        self.root_weights = np.asarray(model.root_likelihood_weights, dtype=float).copy()

        pair = model.tree.edges.reshape(-1, 2, 2)
        self.steps = tuple((int(p), int(l), int(r)) for p, l, r in zip(pair[:, 0, 0], pair[:, 0, 1], pair[:, 1, 1]))
        self.layer_at_node = (
            np.searchsorted(self.times, self.ages, side="right")
            if self.times else np.zeros(self.nn, dtype=int)
        )

        # Build independent Q^(l) and C^(l) per stratum from raw modifiers.
        self.q_layers = []
        self.clado_layers = []
        for lid, mask in enumerate(self.masks):
            dist = (
                np.asarray(model.dist[lid], float) ** params["x"]
                * np.asarray(model.env[lid], float) ** params["n"]
                * np.asarray(model.manual[lid], float) ** params["w"]
            )
            dmat = dist * params["d"]
            elist = (np.asarray(model.area_size[lid], float) ** params["u"]) * params["e"]
            allowed_idx = np.flatnonzero(mask)
            lookup = {frozenset(self.states[idx]): idx for idx in allowed_idx}

            q = np.zeros((self.n, self.n), dtype=float)
            for idx in allowed_idx:
                src = self.states[idx]
                if not src:
                    continue
                occ = frozenset(src)
                for area in range(model.na):
                    if area in occ:
                        dst = lookup.get(occ - {area})
                        rate = float(elist[area])
                    else:
                        dst = lookup.get(occ | {area})
                        rate = float(dmat[list(src), area].sum())
                    if dst is not None:
                        q[idx, dst] += rate
                q[idx, idx] = -q[idx].sum()
            self.q_layers.append(q)

            nonnull_allowed = [int(i) for i in allowed_idx if self.states[i]]
            local_states = tuple(self.states[i] for i in nonnull_allowed)
            _, _, _, weights, jump, *priors = model._layer_arrays(lid, params)
            events = kernel(4, *_states(local_states), weights, jump, *priors,
                            model.na, [], precision="strict_f64")
            clado = np.zeros((self.n, self.n, self.n), dtype=float)
            for anc, lc, rc, w in zip(*events):
                clado[nonnull_allowed[int(anc)], nonnull_allowed[int(lc)], nonnull_allowed[int(rc)]] += float(w)
            for idx in nonnull_allowed:
                row_total = clado[idx].sum()
                if row_total <= 0:
                    raise AssertionError(f"zero cladogenesis row at state {idx} in stratum {lid}")
                clado[idx] /= row_total
            if complete_empty:
                clado[0, 0, 0] = 1.0
            self.clado_layers.append(clado)

        # Exact branch transition matrix P_child = P(X_bot -> X_top) as product of slice exponentials.
        self.branch_p = {}
        for child in range(self.nn):
            par = int(self.parent[child])
            if par < 0:
                continue
            old, young = float(self.ages[par]), float(self.ages[child])
            if self.times:
                cuts = [old, *reversed([x for x in self.times if young < x < old]), young]
            else:
                cuts = [old, young]
            p_mat = np.eye(self.n, dtype=float)
            for start, end in zip(cuts[:-1], cuts[1:]):
                lid = int(np.searchsorted(self.times, 0.5 * (start + end), side="right")) if self.times else 0
                allowed = self.masks[lid].astype(float)
                step_p = expm(self.q_layers[lid] * (start - end))
                step_p[step_p < 0] = 0.0
                # Enforce allowed subspace projection at both ends of the slice.
                step_p = (allowed[:, None] * step_p) * allowed[None, :]
                p_mat = p_mat @ step_p
            self.branch_p[child] = p_mat

    def joint_probability(self, tip_likes: np.ndarray, *,
                          pin_top: tuple[int, int] | None = None,
                          pin_bottom: tuple[int, int] | None = None) -> float:
        """Exact joint probability with optional indicator at a node top or branch bottom."""
        top_vec = {}
        bot_vec = {}
        for i in range(self.nt):
            v = tip_likes[i].astype(float).copy()
            v[~self.masks[self.layer_at_node[i]]] = 0.0
            if pin_top is not None and pin_top[0] == i:
                pinned = np.zeros_like(v)
                pinned[pin_top[1]] = v[pin_top[1]]
                v = pinned
            top_vec[i] = v
            b_vec = self.branch_p[i] @ v
            if pin_bottom is not None and pin_bottom[0] == i:
                pinned = np.zeros_like(b_vec)
                pinned[pin_bottom[1]] = b_vec[pin_bottom[1]]
                b_vec = pinned
            bot_vec[i] = b_vec

        for node, left, right in self.steps:
            lid = int(self.layer_at_node[node])
            clado = self.clado_layers[lid]
            v = np.einsum("aij,i,j->a", clado, bot_vec[left], bot_vec[right])
            v[~self.masks[lid]] = 0.0
            if pin_top is not None and pin_top[0] == node:
                pinned = np.zeros_like(v)
                pinned[pin_top[1]] = v[pin_top[1]]
                v = pinned
            top_vec[node] = v
            if node != self.root:
                b_vec = self.branch_p[node] @ v
                if pin_bottom is not None and pin_bottom[0] == node:
                    pinned = np.zeros_like(b_vec)
                    pinned[pin_bottom[1]] = b_vec[pin_bottom[1]]
                    b_vec = pinned
                bot_vec[node] = b_vec

        return float(np.dot(self.root_weights, top_vec[self.root]))

    def pattern_probabilities(self, patterns: np.ndarray, *, conditioning: str) -> tuple[np.ndarray, float]:
        survival_tips = np.zeros((self.nt, self.n), dtype=float)
        present_nonnull = self.masks[0].copy()
        present_nonnull[0] = False
        survival_tips[:, present_nonnull] = 1.0
        surv_prob = self.joint_probability(survival_tips)

        probs = np.empty(len(patterns), dtype=float)
        eye = np.eye(self.n, dtype=float)
        for k, pat in enumerate(patterns):
            probs[k] = self.joint_probability(eye[pat])
        if conditioning == "tip_survival":
            probs /= surv_prob
        return probs, surv_prob

    def marginal_posteriors(self, tip_likes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        denom = self.joint_probability(tip_likes)
        top = np.zeros((self.nn, self.n), dtype=float)
        bottom = np.full((self.nn, self.n), np.nan, dtype=float)
        for node in range(self.nn):
            for s in range(self.n):
                top[node, s] = self.joint_probability(tip_likes, pin_top=(node, s)) / denom
            if node != self.root:
                for s in range(self.n):
                    bottom[node, s] = self.joint_probability(tip_likes, pin_bottom=(node, s)) / denom
        return top, bottom


TREE4_EXACT = "((T1:0.31234567,T2:0.31234567):0.48765433,(T3:0.52345678,T4:0.52345678):0.27654322);"
GEOG3 = {"T1": "100", "T2": "010", "T3": "001", "T4": "110"}
STRATA_TIMES = [0.25, 0.65, 2.0]
MONOTONE_ALLOWED = [
    [[1, 1, 1], [1, 1, 1], [1, 1, 1]],  # 0.00 - 0.25 Ma: A, B, C all emerged
    [[1, 1, 0], [1, 1, 0], [0, 0, 0]],  # 0.25 - 0.65 Ma: A, B emerged
    [[1, 0, 0], [0, 0, 0], [0, 0, 0]],  # 0.65 - 2.00 Ma: only oldest island A
]
ADJ_BLOCK_AC = [
    [[1, 1, 0], [1, 1, 1], [0, 1, 1]],
    [[1, 1, 0], [1, 1, 1], [0, 1, 1]],
    [[1, 1, 0], [1, 1, 1], [0, 1, 1]],
]
STRATA_DIST = [
    [[1.0, 2.0, 4.5], [2.0, 1.0, 2.2], [4.5, 2.2, 1.0]],
    [[1.0, 2.6, 5.0], [2.6, 1.0, 2.8], [5.0, 2.8, 1.0]],
    [[1.0, 3.2, 6.0], [3.2, 1.0, 3.5], [6.0, 3.5, 1.0]],
]
STRATA_ENV = [
    [[1.0, 1.3, 2.1], [1.3, 1.0, 1.5], [2.1, 1.5, 1.0]],
    [[1.0, 1.6, 2.4], [1.6, 1.0, 1.7], [2.4, 1.7, 1.0]],
    [[1.0, 1.9, 2.8], [1.9, 1.0, 2.0], [2.8, 2.0, 1.0]],
]
STRATA_MANUAL = [
    [[1.0, 0.85, 0.35], [0.90, 1.0, 0.75], [0.25, 0.65, 1.0]],
    [[1.0, 0.55, 0.05], [0.50, 1.0, 0.05], [0.05, 0.05, 1.0]],
    [[1.0, 0.05, 0.05], [0.05, 1.0, 0.05], [0.05, 0.05, 1.0]],
]
STRATA_AREA = [
    [1.6, 1.1, 0.7],
    [1.4, 0.8, 0.5],
    [1.2, 0.5, 0.5],
]


def build_full_stratified_model(family: str, contract: str, *, allow_null_tips: bool | None = None,
                                root_prior=None, d=0.18, e=0.11, j=0.14, x=-0.4, w=1.2, n=-0.3, u=-0.5):
    has_j = "+J" in family.upper()
    tpl = b.model_template(f"{family}+x+w+n+u").with_parameters(
        d=dict(type="fixed", init=d),
        e=dict(type="fixed", init=e),
        j=dict(type="fixed", init=j if has_j else 0.0),
        x=dict(type="fixed", init=x),
        w=dict(type="fixed", init=w),
        n=dict(type="fixed", init=n),
        u=dict(type="fixed", init=u),
    )
    if allow_null_tips is None:
        allow_null_tips = contract == "unconditional_complete"
    return b.BioGeoBEARS(
        TREE4_EXACT, GEOG3, tpl,
        timeperiods=STRATA_TIMES,
        areas_allowed=MONOTONE_ALLOWED,
        areas_adjacency=ADJ_BLOCK_AC,
        distances=STRATA_DIST,
        envdistances=STRATA_ENV,
        dispersal_multipliers=STRATA_MANUAL,
        area_of_areas=STRATA_AREA,
        likelihood_conditioning=contract,
        allow_null_tips=allow_null_tips,
        root_prior=root_prior,
    )


class IndependentStratifiedReferenceTests(unittest.TestCase):
    def test_single_layer_all_spatial_modifiers_exact_laws_and_marginals(self):
        # Single-layer model with areas_adjacency (excluding AC), distances (+x),
        # envdistances (+n), dispersal_multipliers (+w), and area_of_areas (+u).
        for family in ("DEC", "DEC+J", "DIVALIKE+J", "BAYAREALIKE+J"):
            with self.subTest(family=family):
                has_j = "+J" in family
                tpl = b.model_template(f"{family}+x+w+n+u").with_parameters(
                    d=dict(type="fixed", init=0.16),
                    e=dict(type="fixed", init=0.12),
                    j=dict(type="fixed", init=0.15 if has_j else 0.0),
                    x=dict(type="fixed", init=-0.45),
                    w=dict(type="fixed", init=1.15),
                    n=dict(type="fixed", init=-0.25),
                    u=dict(type="fixed", init=-0.4),
                )
                common = dict(
                    areas_adjacency=ADJ_BLOCK_AC[0],
                    distances=STRATA_DIST[0],
                    envdistances=STRATA_ENV[0],
                    dispersal_multipliers=STRATA_MANUAL[0],
                    area_of_areas=STRATA_AREA[0],
                )
                ma = b.BioGeoBEARS(TREE4_EXACT, GEOG3, tpl, likelihood_conditioning="tip_survival", **common)
                mb = b.BioGeoBEARS(TREE4_EXACT, GEOG3, tpl, likelihood_conditioning="unconditional_complete",
                                   allow_null_tips=True, **common)
                self.assertNotIn("AC", ma.ranges)

                ref_b = IndependentStratifiedReference(mb, mb.model.resolve(), complete_empty=True)
                pats_b = np.array(list(itertools.product(range(len(mb.states)), repeat=4)), dtype=np.int64)
                pats_a = np.array(list(itertools.product(range(1, len(ma.states)), repeat=4)), dtype=np.int64)

                exp_b, surv = ref_b.pattern_probabilities(pats_b, conditioning="unconditional_complete")
                exp_a, _ = ref_b.pattern_probabilities(pats_a, conditioning="tip_survival")
                with mb.prepared_likelihood() as pb:
                    obs_b = np.exp(pb.evaluate_tip_states(pats_b))
                with ma.prepared_likelihood() as pa:
                    obs_a = np.exp(pa.evaluate_tip_states(pats_a))

                self.assertAlmostEqual(float(obs_b.sum()), 1.0, places=12)
                self.assertAlmostEqual(float(obs_a.sum()), 1.0, places=12)
                np.testing.assert_allclose(obs_b, exp_b, atol=1e-12, rtol=0)
                np.testing.assert_allclose(obs_a, exp_a, atol=1e-12, rtol=0)

                ev_a = ma.evaluate()
                self.assertAlmostEqual(ev_a.log_tip_survival_prob, math.log(surv), places=12)
                ref_top, ref_bot = ref_b.marginal_posteriors(ma.tips)
                np.testing.assert_allclose(ev_a.ancestral_probs_branch_top, ref_top, atol=1e-12, rtol=0)
                non_root = [i for i in range(len(ma.tree.labels)) if i != ma.tree.root]
                np.testing.assert_allclose(ev_a.ancestral_probs_branch_bottom[non_root],
                                           ref_bot[non_root], atol=1e-12, rtol=0)

    def test_stratified_monotone_emergence_all_six_families_exact_reference(self):
        families = ("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
        for family in families:
            with self.subTest(family=family):
                ma = build_full_stratified_model(family, "tip_survival")
                mb = build_full_stratified_model(family, "unconditional_complete")
                # Oldest stratum only allows empty (0) and island A (1): uniform_nonnull puts mass 1.0 on A.
                expected_root = np.zeros(len(ma.states))
                expected_root[ma.ranges.index("A")] = 1.0
                np.testing.assert_array_equal(ma.root_sampling_probs, expected_root)

                valid_0 = np.flatnonzero(ma.masks[0])
                valid_nonnull_0 = valid_0[valid_0 != 0]
                pats_b = np.array(list(itertools.product(valid_0, repeat=4)), dtype=np.int64)
                pats_a = np.array(list(itertools.product(valid_nonnull_0, repeat=4)), dtype=np.int64)

                ref = IndependentStratifiedReference(mb, mb.model.resolve(), complete_empty=True)
                exp_b, surv = ref.pattern_probabilities(pats_b, conditioning="unconditional_complete")
                exp_a, _ = ref.pattern_probabilities(pats_a, conditioning="tip_survival")

                with mb.prepared_likelihood() as pb:
                    obs_b = np.exp(pb.evaluate_tip_states(pats_b))
                with ma.prepared_likelihood() as pa:
                    obs_a = np.exp(pa.evaluate_tip_states(pats_a))

                self.assertAlmostEqual(float(exp_b.sum()), 1.0, places=12)
                self.assertAlmostEqual(float(obs_b.sum()), 1.0, places=12)
                self.assertAlmostEqual(float(obs_a.sum()), 1.0, places=12)
                np.testing.assert_allclose(obs_b, exp_b, atol=1e-12, rtol=0)
                np.testing.assert_allclose(obs_a, exp_a, atol=1e-12, rtol=0)

                ev_a = ma.evaluate()
                ev_b = mb.evaluate()
                self.assertAlmostEqual(ev_a.log_tip_survival_prob, math.log(surv), places=12)
                self.assertAlmostEqual(ev_a.lnL, ev_b.lnL - ev_a.log_tip_survival_prob, places=13)

                # Verify exact marginal posteriors at every node (including empty-tip pattern under B).
                ref_top_a, ref_bot_a = ref.marginal_posteriors(ma.tips)
                np.testing.assert_allclose(ev_a.ancestral_probs_branch_top, ref_top_a, atol=1e-12, rtol=0)
                non_root = [i for i in range(len(ma.tree.labels)) if i != ma.tree.root]
                np.testing.assert_allclose(ev_a.ancestral_probs_branch_bottom[non_root],
                                           ref_bot_a[non_root], atol=1e-12, rtol=0)

                # Also verify Contract-B marginal posteriors when a clade has an empty tip.
                geog_with_null = b.TipRanges(("A", "B", "C"), {"T1": (), "T2": (0,), "T3": (1,), "T4": (0, 1)})
                mb_null = b.BioGeoBEARS(
                    TREE4_EXACT, geog_with_null, mb.model,
                    timeperiods=STRATA_TIMES, areas_allowed=MONOTONE_ALLOWED,
                    areas_adjacency=ADJ_BLOCK_AC, distances=STRATA_DIST,
                    envdistances=STRATA_ENV, dispersal_multipliers=STRATA_MANUAL,
                    area_of_areas=STRATA_AREA, likelihood_conditioning="unconditional_complete",
                    allow_null_tips=True,
                )
                ev_b_null = mb_null.evaluate()
                ref_top_b, ref_bot_b = ref.marginal_posteriors(mb_null.tips)
                np.testing.assert_allclose(ev_b_null.ancestral_probs_branch_top, ref_top_b, atol=1e-12, rtol=0)
                np.testing.assert_allclose(ev_b_null.ancestral_probs_branch_bottom[non_root],
                                           ref_bot_b[non_root], atol=1e-12, rtol=0)

    def test_structural_zero_patterns_under_stratified_dec_return_clean_neginf(self):
        # With d=0 and j=0 starting on island A in the oldest stratum, any tip on B or C has exact 0 mass.
        ma = build_full_stratified_model("DEC", "tip_survival", d=0.0, e=0.1, j=0.0)
        ev = ma.evaluate(calc_ancprobs=False)
        self.assertTrue(np.isneginf(ev.lnL))
        valid_nonnull_0 = np.flatnonzero(ma.masks[0])[1:]
        pats_a = np.array(list(itertools.product(valid_nonnull_0, repeat=4)), dtype=np.int64)
        with ma.prepared_likelihood() as pa:
            lnliks = pa.evaluate_tip_states(pats_a)
        probs = np.exp(lnliks)
        self.assertFalse(np.isnan(lnliks).any())
        # Only the all-A pattern (1, 1, 1, 1) is reachable, so its conditional probability must be 1.0!
        a_idx = ma.ranges.index("A")
        all_a = np.all(pats_a == a_idx, axis=1)
        self.assertAlmostEqual(float(probs[all_a][0]), 1.0, places=14)
        self.assertEqual(float(probs[~all_a].sum()), 0.0)


class StratifiedSimulationAndBSMDualityTests(unittest.TestCase):
    def test_forward_simulation_duality_and_batch_thread_determinism(self):
        for contract in ("tip_survival", "unconditional_complete"):
            with self.subTest(contract=contract):
                m = build_full_stratified_model("DEC+J", contract)
                valid_0 = np.flatnonzero(m.masks[0])
                if contract == "tip_survival":
                    valid_0 = valid_0[valid_0 != 0]
                pats = np.array(list(itertools.product(valid_0, repeat=4)), dtype=np.int64)
                with m.prepared_likelihood() as pl:
                    expected = np.exp(pl.evaluate_tip_states(pats))

                n_reps = 80_000
                with b.PreparedSimulation(m) as sim:
                    full = sim.draw(n_reps, seed=20261118, num_threads=4)
                    part1 = sim.draw(333, start_rep=1, seed=20261118, num_threads=1)
                    part2 = sim.draw(667, start_rep=334, seed=20261118, num_threads=3)
                np.testing.assert_array_equal(
                    full.tip_states[:1000],
                    np.concatenate([part1.tip_states, part2.tip_states]),
                )
                np.testing.assert_array_equal(
                    full.node_states[:1000],
                    np.concatenate([part1.node_states, part2.node_states]),
                )

                lookup = {tuple(p): i for i, p in enumerate(pats)}
                counts = np.bincount([lookup[tuple(row)] for row in full.tip_states], minlength=len(pats))
                # Simultaneous Bernstein confidence band across all patterns (alpha = 1e-6).
                logterm = math.log(2.0 * len(pats) / 1e-6)
                limit = np.sqrt(2.0 * expected * (1.0 - expected) * logterm / n_reps) + 2.0 * logterm / (3.0 * n_reps)
                self.assertEqual(int(counts[expected == 0].sum()), 0)
                self.assertTrue(np.all(np.abs(counts / n_reps - expected) <= limit))

    def test_empirical_posterior_root_simulation_under_partial_stratum_emergence(self):
        # Oldest stratum allows A and B (and AB); verify empirical_posterior root draws match root marginal.
        allowed_two = [
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
        ]
        tpl = b.model_template("DEC+J").with_parameters(
            d=dict(type="fixed", init=0.2), e=dict(type="fixed", init=0.08), j=dict(type="fixed", init=0.12)
        )
        fit = b.BioGeoBEARS(
            TREE4_EXACT, GEOG3, tpl, timeperiods=STRATA_TIMES, areas_allowed=allowed_two,
            likelihood_conditioning="tip_survival",
        ).fit()
        expected_root = fit.evaluation.ancestral_probs_branch_top[fit.model.tree.root]
        self.assertEqual(float(expected_root[fit.model.ranges.index("C")]), 0.0)
        with b.PreparedSimulation(fit, root_distribution="empirical_posterior") as sim:
            draws = sim.draw(40_000, seed=77, num_threads=4)
        emp = np.bincount(draws.node_states[:, fit.model.tree.root], minlength=len(expected_root)) / 40_000
        np.testing.assert_allclose(emp, expected_root, atol=0.01, rtol=0)

    def test_stratified_bsm_respects_geological_emergence_and_node_posteriors(self):
        allowed_two = [
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],  # stratum 1 (youngest): A, B, C
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],  # stratum 2: A, B only
            [[1, 0, 0], [0, 0, 0], [0, 0, 0]],  # stratum 3 (oldest): A only
        ]
        tpl = b.model_template("DEC+J+x").with_parameters(
            d=dict(type="fixed", init=0.22), e=dict(type="fixed", init=0.09),
            j=dict(type="fixed", init=0.16), x=dict(type="fixed", init=-0.35),
        )
        fit = b.BioGeoBEARS(
            TREE4_EXACT, GEOG3, tpl,
            timeperiods=STRATA_TIMES, areas_allowed=allowed_two,
            areas_adjacency=ADJ_BLOCK_AC, distances=STRATA_DIST,
            likelihood_conditioning="tip_survival",
        ).fit()

        # 1) Under rare two-step colonization of newly emerged island C (0..0.25 Ma),
        # exact uniformization completes 8,000 maps effortlessly without rejection exhaustion!
        bsm = b.runBSM(fit, nummaps_goal=8000, bsm_mode="uniformization", num_threads=4, seedval=909)
        clado = bsm.clado_events_df
        ana = bsm.ana_events_df
        s3_clado = clado[clado["stratum"] == 3]
        for col in ("ancestral_range", "left_range", "right_range"):
            self.assertTrue(set(s3_clado[col]).issubset({"A"}), (col, set(s3_clado[col])))
        self.assertEqual(len(ana[ana["stratum"] == 3]), 0)

        s2_clado = clado[clado["stratum"] == 2]
        for col in ("ancestral_range", "left_range", "right_range"):
            self.assertFalse(any("C" in r for r in s2_clado[col]))
        s2_ana = ana[ana["stratum"] == 2]
        for col in ("from_range", "to_range"):
            self.assertFalse(any("C" in r for r in s2_ana[col]))
        self.assertNotIn("C", set(s2_ana["destination_area"].dropna()))

        for col in ("ancestral_range", "left_range", "right_range"):
            self.assertNotIn("AC", set(clado[col]))
        for col in ("from_range", "to_range"):
            self.assertNotIn("AC", set(ana[col]))

        for internal_node in (5, 6, 7):
            sub = clado[clado["node"] == internal_node]
            emp = np.bincount(sub["ancestral_state"].to_numpy(), minlength=len(fit.model.states)) / 8000
            exact = fit.evaluation.ancestral_probs_branch_top[internal_node - 1]
            np.testing.assert_allclose(emp, exact, atol=0.018, rtol=0)

        # 2) Verify legacy rejection explicitly detects 4-step stepping-stone exhaustion
        # while succeeding when direct colonization is non-rare:
        with self.assertRaisesRegex(ValueError, "use uniformization"):
            b.runBSM(fit, nummaps_goal=400, bsm_mode="rejection", num_threads=2, seedval=909)
        fit_mod = b.BioGeoBEARS(
            TREE4_EXACT, {"T1": "100", "T2": "010", "T3": "011", "T4": "110"},
            b.model_template("DEC+J").with_parameters(
                d=dict(type="fixed", init=0.45), e=dict(type="fixed", init=0.12), j=dict(type="fixed", init=0.2)
            ),
            timeperiods=STRATA_TIMES, areas_allowed=allowed_two,
            areas_adjacency=ADJ_BLOCK_AC, likelihood_conditioning="tip_survival",
        ).fit()
        bsm_rej = b.runBSM(fit_mod, nummaps_goal=200, bsm_mode="rejection", num_threads=4, seedval=909)
        self.assertNotIn("AC", set(bsm_rej.clado_events_df["ancestral_range"]))
        self.assertTrue(set(bsm_rej.clado_events_df[bsm_rej.clado_events_df["stratum"] == 3]["ancestral_range"]).issubset({"A"}))


class NegativeContractGateTests(unittest.TestCase):
    def test_non_monotone_area_sinking_is_rejected(self):
        # Island B exists in older stratum 2, then vanishes in younger stratum 1 -> violates forward conservation!
        sinking_allowed = [
            [[1, 0, 0], [0, 0, 0], [0, 0, 1]],
            [[1, 1, 0], [1, 1, 0], [0, 0, 1]],
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
        ]
        geog_ac = {"T1": "100", "T2": "100", "T3": "001", "T4": "001"}
        with self.assertRaisesRegex(b.UnsupportedModelContractError, "monotone"):
            b.BioGeoBEARS(
                TREE4_EXACT, geog_ac, "DEC+J",
                timeperiods=STRATA_TIMES, areas_allowed=sinking_allowed,
                likelihood_conditioning="tip_survival",
            )

    def test_non_binary_allowed_matrix_and_invalid_root_and_tip_constraints_rejected(self):
        bad_allowed = [
            [[1, 0.5, 1], [0.5, 1, 1], [1, 1, 1]],
            [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
        ]
        with self.assertRaisesRegex(b.UnsupportedModelContractError, "binary"):
            b.BioGeoBEARS(
                TREE4_EXACT, GEOG3, "DEC",
                timeperiods=STRATA_TIMES, areas_allowed=bad_allowed,
                likelihood_conditioning="tip_survival",
            )

        # Custom root prior assigning positive mass to unemerged island C in oldest stratum -> ValueError.
        with self.assertRaisesRegex(ValueError, "disallowed in the root stratum"):
            build_full_stratified_model(
                "DEC+J", "tip_survival",
                root_prior=b.RootPriorSpec("custom_weights", {"A": 0.5, "C": 0.5}),
            )

        # Tip T3 observed on C when C is disallowed even in the present stratum -> UnsupportedModelContractError.
        present_blocks_c = [
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
            [[1, 1, 0], [1, 1, 0], [0, 0, 0]],
            [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
        ]
        with self.assertRaisesRegex(b.UnsupportedModelContractError, "present stratum"):
            b.BioGeoBEARS(
                TREE4_EXACT, GEOG3, "DEC+J",
                timeperiods=STRATA_TIMES, areas_allowed=present_blocks_c,
                likelihood_conditioning="tip_survival",
            )

        # Free w with negative lower bound and zero manual multiplier -> ValueError; nonnegative min succeeds!
        zero_manual = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] * 3
        with self.assertRaisesRegex(ValueError, "strictly positive"):
            b.BioGeoBEARS(
                TREE4_EXACT, {"T1": "100", "T2": "100", "T3": "100", "T4": "100"},
                "DEC+w", timeperiods=STRATA_TIMES, dispersal_multipliers=zero_manual,
                likelihood_conditioning="tip_survival",
            )
        bounded_w = b.model_template("DEC+w").with_parameters(w=dict(min=0.0, init=1.0, max=5.0))
        ok_model = b.BioGeoBEARS(
            TREE4_EXACT, {"T1": "100", "T2": "100", "T3": "100", "T4": "100"},
            bounded_w, timeperiods=STRATA_TIMES, dispersal_multipliers=zero_manual,
            likelihood_conditioning="tip_survival",
        )
        self.assertTrue(np.isfinite(ok_model.evaluate().lnL))

    def test_tier3_features_remain_explicitly_blocked(self):
        # Fixed internal node and rootedge remain rejected under probability contracts.
        with self.assertRaises(b.UnsupportedModelContractError):
            b.BioGeoBEARS(
                TREE4_EXACT, GEOG3, "DEC",
                timeperiods=STRATA_TIMES, areas_allowed=MONOTONE_ALLOWED,
                fixnode=[5], fixlikes=[[0, 1, 0, 0, 0, 0, 0, 0]],
                likelihood_conditioning="tip_survival",
            )
        with self.assertRaises(b.UnsupportedModelContractError):
            b.BioGeoBEARS(
                TREE4_EXACT, GEOG3, "DEC",
                rootedge=True, likelihood_conditioning="tip_survival",
            )


class StratifiedUncertaintyWorkflowsEndToEndTests(unittest.TestCase):
    def test_profile_bootstrap_averaging_multitree_and_diagnostics_on_stratified_emergence(self):
        tree6 = "(((T1:0.22,T2:0.22):0.25,T3:0.47):0.43,((T4:0.18,T5:0.18):0.34,T6:0.52):0.38);"
        geog6 = {"T1": "100", "T2": "100", "T3": "010", "T4": "010", "T5": "001", "T6": "001"}
        times = [0.30, 0.65, 2.0]

        common = dict(
            timeperiods=times,
            areas_allowed=MONOTONE_ALLOWED,
            areas_adjacency=ADJ_BLOCK_AC,
            distances=STRATA_DIST,
            dispersal_multipliers=STRATA_MANUAL,
        )
        fit_dec_c = b.BioGeoBEARS(tree6, geog6, "DEC", likelihood_conditioning="legacy_unconditioned", **common).fit()
        fit_decj_c = b.BioGeoBEARS(tree6, geog6, "DEC+J", likelihood_conditioning="legacy_unconditioned", **common).fit()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            # 1. Cross-contract C -> A conversion + parametric bootstrap LRT + checkpoint resume.
            boot_ckpt = tmp_path / "boot.jsonl"
            boot1 = b.bootstrap_lrt(
                fit_dec_c, fit_decj_c, n_boot=6,
                target_conditioning="tip_survival", refit_observed=True,
                calc_surrogate_power=True, n_alt_sim=4,
                seed=101, num_threads=2, checkpoint_file=boot_ckpt,
            )
            self.assertTrue(boot1.recomputed_observed_under_new_contract)
            self.assertEqual(boot1.likelihood_conditioning, "tip_survival")
            self.assertEqual(boot1.n_converged + boot1.n_failed, 6)

            # Simulate interrupted checkpoint (keep header + first 4 committed records) and resume:
            ckpt_lines = boot_ckpt.read_bytes().splitlines(keepends=True)
            boot_ckpt.write_bytes(b"".join(ckpt_lines[:5]))
            boot2 = b.bootstrap_lrt(
                fit_dec_c, fit_decj_c, n_boot=6,
                target_conditioning="tip_survival", refit_observed=True,
                calc_surrogate_power=True, n_alt_sim=4,
                seed=101, num_threads=2, checkpoint_file=boot_ckpt, resume=True,
            )
            self.assertEqual(len(boot2.replicates_df), 6)
            np.testing.assert_allclose(
                boot1.replicates_df["d_stat"].to_numpy(),
                boot2.replicates_df["d_stat"].to_numpy(),
                atol=1e-12, rtol=0,
            )

            # 2. Profile CI (including +x spatial distance exponent) + 2D profile surface.
            fit_decjx = b.BioGeoBEARS(
                tree6, geog6, "DEC+J+x", likelihood_conditioning="tip_survival", **common
            ).fit()
            prof_ckpt = tmp_path / "prof.jsonl"
            prof = b.profile_ci(
                fit_decjx, params=["j", "x"], n_grid=6, tol=1e-3,
                num_threads=2, checkpoint_file=prof_ckpt,
            )
            df_prof = prof.to_dataframe()
            self.assertEqual(list(df_prof["param"]), ["j", "x"])
            self.assertTrue(np.all(df_prof["endpoint_deviance_residual_max"] < 1e-3))

            surf = b.profile_surface(
                fit_decjx, param_x="d", param_y="j", grid_x=3, grid_y=3, num_threads=2,
            )
            self.assertEqual(len(surf), 9)
            self.assertTrue(np.isfinite(surf["profile_obj"]).all())

            # 3. Model averaging across multiple candidate families + stratified BSM + diagnostics.
            fit_divj = b.BioGeoBEARS(
                tree6, geog6, "DIVALIKE+J", likelihood_conditioning="tip_survival", **common
            ).fit()
            ma_res = b.model_average([boot2.observed_null_fit, boot2.observed_alt_fit, fit_divj], criterion="AIC")
            self.assertAlmostEqual(float(ma_res.weights_df["weight"].sum()), 1.0, places=12)
            np.testing.assert_allclose(ma_res.ancestral_probs_branch_top.sum(axis=1), 1.0, atol=1e-12)
            ma_bsm = ma_res.run_bsm(num_maps=40, num_threads=2, seed=303)
            diag_single = b.bsm_diagnostics(ma_bsm)
            self.assertEqual(diag_single.num_maps_total, 40)
            self.assertTrue(len(diag_single.event_summary_df) > 0)
            # Forbidden range 'AC' and premature C in stratum 3 must remain absent in averaged BSM!
            self.assertNotIn("AC", set(ma_bsm.clado_events_df["ancestral_range"]))
            s3 = ma_bsm.clado_events_df[ma_bsm.clado_events_df["stratum"] == 3]
            self.assertTrue(set(s3["ancestral_range"]).issubset({"A"}))

            # 4. Multi-tree stratified fitting (including root-height variation across posterior trees) + posterior diagnostics.
            trees_sample = [
                tree6,
                "(((T1:0.20,T2:0.20):0.26,T3:0.46):0.41,((T4:0.19,T5:0.19):0.33,T6:0.52):0.35);",
                "(((T2:0.24,T1:0.24):0.24,T3:0.48):0.47,((T4:0.17,T5:0.17):0.36,T6:0.53):0.42);",
                "((T1:0.45,(T2:0.21,T3:0.21):0.24):0.44,((T4:0.18,T5:0.18):0.34,T6:0.52):0.37);",
                "(((T1:0.22,T2:0.22):0.25,T3:0.47):0.33,((T4:0.18,T5:0.18):0.32,T6:0.50):0.30);",
                "(((T1:0.23,T2:0.23):0.26,T3:0.49):0.45,((T4:0.20,T5:0.20):0.34,T6:0.54):0.40);",
            ]
            mt_ckpt = tmp_path / "mt.jsonl"
            mt = b.fit_over_trees(
                trees_sample, geog6, models=["DEC", "DEC+J"],
                likelihood_conditioning="tip_survival",
                criterion="AIC", bsm_maps_per_tree=10,
                num_threads=2, seed=505, checkpoint_file=mt_ckpt,
                **common,
            )
            self.assertEqual(mt.n_trees_completed, 6)
            diag_pop = b.bsm_diagnostics(
                mt.joint_bsm, tree_target="posterior_population", tree_dependence="independent",
            )
            self.assertEqual(diag_pop.num_maps_total, 60)
            self.assertTrue(len(diag_pop.event_summary_df) > 0)


if __name__ == "__main__":
    unittest.main()
