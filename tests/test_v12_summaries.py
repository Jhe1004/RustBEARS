"""Precision targets, model/tree alignment and lossless recovery contracts."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import rustbears as b
from rustbears.bsm import _ANA, _CLADO
from rustbears.uncertainty.averaging import model_average, align_node_array
from rustbears.uncertainty.checkpoint import (CheckpointStore, CheckpointMismatchError,
    task_fingerprint, frame_arrays, restore_frame)
from rustbears.uncertainty.diagnostics import (_make_bsm, bsm_diagnostics,
    run_bsm_adaptive, rank_status)


def fitted(tree="((T1:1,T2:1):1,T3:2);", d=.2, j=0):
    template = b.model_template("DEC+J" if j else "DEC").with_parameters(
        d=dict(type="fixed", init=d), e=dict(type="fixed", init=.1), j=dict(type="fixed", init=j))
    return b.BioGeoBEARS(tree, {"T1":"10", "T2":"10", "T3":"01"}, template,
                         likelihood_conditioning="tip_survival").fit()


def synthetic(forward, reverse=None, metadata=None):
    if reverse is None:
        reverse = np.zeros(len(forward), int)
    rows = []
    for i, (a, z) in enumerate(zip(forward, reverse), 1):
        for source, destination, count in (("A", "B", a), ("B", "A", z)):
            for _ in range(int(count)):
                rows.append([i, 1, 1, 1, 2, .5, .5, "d", source, destination])
    return _make_bsm(pd.DataFrame(columns=_CLADO), pd.DataFrame(rows, columns=_ANA),
                     ("A", "B"), len(forward), metadata)


class CheckpointTests(unittest.TestCase):
    def test_manifest_detects_float_array_and_parameter_changes(self):
        data = {"x":np.array([0., 1.]), "parameters":{"d":.1}, "seed":5}
        for other in ({**data, "x":np.array([-0., 1.])}, {**data, "seed":6},
                      {**data, "parameters":{"d":.2}}):
            self.assertNotEqual(task_fingerprint(data), task_fingerprint(other))

    def test_atomic_record_roundtrip_and_corruption(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"state.jsonl"
            store = CheckpointStore(path, {"seed":1})
            frame = pd.DataFrame({"n":[1,2], "source":[None,"A"], "time":[.1,.2]})
            meta, arrays = frame_arrays(frame, "events")
            store.append({"id":1, "status":"failed", "score":-np.inf, "events":meta}, arrays)
            with path.open("ab") as h:
                h.write(b'{"unfinished":')
            resumed = CheckpointStore(path, {"seed":1}, resume=True)
            self.assertEqual(resumed.completed_ids, (1,))
            self.assertEqual(resumed.records[0]["score"], -np.inf)
            pd.testing.assert_frame_equal(frame, restore_frame(meta, resumed.load_arrays(1)))
            with self.assertRaises(CheckpointMismatchError):
                CheckpointStore(path, {"seed":2}, resume=True)
            shard = next(resumed.shards.glob("*.npz"))
            shard.write_bytes(b"corrupted")
            with self.assertRaises(CheckpointMismatchError):
                CheckpointStore(path, {"seed":1}, resume=True)


class DiagnosticTests(unittest.TestCase):
    def test_single_tree_mcse_and_realization_interval_are_distinct(self):
        counts = [0,2,4,6]
        result = bsm_diagnostics(synthetic(counts), min_mean_count=.1)
        row = result.route_summary_df.iloc[0]
        self.assertAlmostEqual(row["mean"], 3)
        self.assertAlmostEqual(row.se_target, np.std(counts, ddof=1)/2)
        self.assertEqual((row.rpi_025, row.rpi_975), (0,6))
        self.assertEqual(result.stop_reason, "fixed_count")
        unseen = result.route_summary_df.iloc[1]
        self.assertEqual(unseen.assessment, "not_assessed")
        self.assertTrue(np.isnan(unseen.mean_mc_ci_upper))

    def test_weighted_strata_use_squared_weights_not_pooled_variance(self):
        metadata = pd.DataFrame({"map":range(1,7), "tree_id":1,
            "model_id":[1,1,2,2,2,2], "stratum_weight":[.8,.8,.2,.2,.2,.2]})
        result = bsm_diagnostics(synthetic([1,3,9,11,9,11], metadata=metadata), min_mean_count=.1)
        row = result.route_summary_df.iloc[0]
        expected = np.sqrt(.8**2*np.var([1,3],ddof=1)/2 + .2**2*np.var([9,11,9,11],ddof=1)/4)
        self.assertAlmostEqual(row["mean"], 3.6)
        self.assertAlmostEqual(row.se_target, expected)
        self.assertEqual(result.sampling_design, "model_stratified")

    def test_fixed_sample_and_population_have_different_targets(self):
        meta = pd.DataFrame({"map":range(1,7), "tree_id":[1,1,2,2,3,3],
                              "model_id":1, "stratum_weight":1/3})
        maps = synthetic([1,1,4,4,7,7], metadata=meta)
        fixed = bsm_diagnostics(maps, min_mean_count=.1)
        population = bsm_diagnostics(maps, min_mean_count=.1, tree_target="posterior_population", tree_dependence="independent")
        self.assertEqual(fixed.route_summary_df.iloc[0].se_target, 0)
        self.assertAlmostEqual(population.route_summary_df.iloc[0].se_target, np.std([1,4,7],ddof=1)/np.sqrt(3))
        with self.assertRaisesRegex(ValueError, "explicit"):
            bsm_diagnostics(maps, tree_target="posterior_population")

    def test_ess_is_scalar_specific_and_chain_order_is_preserved(self):
        n = 40
        forward = np.repeat(np.r_[np.full(20,2),np.full(20,8)],2)
        reverse = np.repeat(np.tile([2,8],20),2)
        meta = pd.DataFrame({"map":range(1,81), "tree_id":np.repeat(range(1,41),2),
                             "model_id":1, "stratum_weight":1/n})
        maps = synthetic(forward, reverse, meta)
        result = bsm_diagnostics(maps, min_mean_count=.1, tree_target="posterior_population", tree_dependence="mcmc_chain")
        ess = result.route_summary_df.effective_n_trees.to_numpy()
        self.assertLess(ess[0], ess[1])
        maps.tree_order = tuple(np.ravel(np.column_stack((np.arange(1,21),np.arange(21,41)))))
        altered = bsm_diagnostics(maps, min_mean_count=.1, tree_target="posterior_population", tree_dependence="mcmc_chain")
        self.assertNotEqual(ess[0], altered.route_summary_df.iloc[0].effective_n_trees)

    def test_paired_contrast_retains_covariance_and_equivalence_rules(self):
        maps = synthetic([1,2,3,4], [1,2,3,4])
        result = bsm_diagnostics(maps, min_mean_count=.1)
        row = result.route_rank_comparisons_df.iloc[0]
        self.assertEqual(row.se_diff, 0)
        self.assertEqual(row.status, "near_tie")
        for lower, upper, expected in ((-.2,.2,"near_tie"),(.3,.5,"resolved"),(-.5,-.3,"resolved"),(-.4,.4,"unresolved"),(.1,.4,"unresolved")):
            self.assertEqual(rank_status(lower,upper,.25),expected)


class ModelAverageTests(unittest.TestCase):
    def test_rotated_tree_nodes_are_aligned_before_averaging(self):
        left = fitted()
        right = fitted("(T3:2,(T2:1,T1:1):1);")
        result = model_average([left,right])
        np.testing.assert_allclose(result.ancestral_probs_branch_top,left.evaluation.ancestral_probs_branch_top,atol=1e-14)
        self.assertTrue(np.allclose(result.weights_df.weight,[.5,.5]))

    def test_aicc_unavailable_is_not_silently_dropped(self):
        f = fitted()
        with self.assertRaisesRegex(ValueError,"undefined AICc"):
            model_average([f,replace(f,aicc=np.inf)])

    def test_stratified_sampling_is_repeatable_and_carries_every_map(self):
        averaged = model_average([fitted(),fitted(j=.1)])
        a, z = averaged.run_bsm(20,seed=3), averaged.run_bsm(20,seed=3,num_threads=2)
        pd.testing.assert_frame_equal(a.clado_events_df,z.clado_events_df)
        pd.testing.assert_frame_equal(a.ana_events_df,z.ana_events_df)
        self.assertEqual(len(a.map_metadata),20)
        self.assertEqual(set(a.map_metadata.model_id),{1,2})
        self.assertEqual(bsm_diagnostics(a).sampling_design,"model_stratified")

    def test_adaptive_resume_restores_every_event_without_resampling(self):
        with tempfile.TemporaryDirectory() as temp:
            path = str(Path(temp)/"adaptive.jsonl")
            first = run_bsm_adaptive(fitted(),min_maps=6,batch_size=3,max_maps=9,
                                     target_rmcse=1e-12,checkpoint_file=path,seed=7)
            with patch("rustbears._ffi.PreparedSampler",side_effect=AssertionError("resampled")):
                second = run_bsm_adaptive(fitted(),min_maps=6,batch_size=3,max_maps=9,
                                     target_rmcse=1e-12,checkpoint_file=path,resume=True,seed=7)
            pd.testing.assert_frame_equal(first.clado_events_df,second.clado_events_df)
            pd.testing.assert_frame_equal(first.ana_events_df,second.ana_events_df)
            self.assertEqual(first.stop_reason,second.stop_reason)

    def test_batches_50_and_200_are_bit_identical_at_fixed_total(self):
        fit=fitted()
        small=run_bsm_adaptive(fit,min_maps=400,batch_size=50,max_maps=400,target_rmcse=1e-12,seed=17)
        large=run_bsm_adaptive(fit,min_maps=400,batch_size=200,max_maps=400,target_rmcse=1e-12,seed=17,num_threads=2)
        pd.testing.assert_frame_equal(small.clado_events_df,large.clado_events_df,check_exact=True)
        pd.testing.assert_frame_equal(small.ana_events_df,large.ana_events_df,check_exact=True)

    def test_legacy_empty_tip_score_is_not_silently_completed_for_new_bsm(self):
        template=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=.2),e=dict(type="fixed",init=.1))
        fit=b.BioGeoBEARS("((T1:1,T2:1):1,T3:2);",{"T1":"00","T2":"10","T3":"01"},template,
                         allow_null_tips=True).fit()
        with self.assertRaisesRegex(b.UnsupportedModelContractError,"empty"):
            model_average([fit]).run_bsm(4)
        with self.assertRaisesRegex(b.UnsupportedModelContractError,"empty"):
            run_bsm_adaptive(fit,min_maps=2,batch_size=2,max_maps=4)


if __name__ == "__main__":
    unittest.main()
