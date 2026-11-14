"""Actual small-tree fits, topology matching, injected failures and recovery."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import warnings

import numpy as np
import pandas as pd

import rustbears as b
from rustbears.uncertainty.averaging import node_clades, topology_hash
from rustbears.uncertainty.checkpoint import CheckpointStore, CheckpointMismatchError
from rustbears.uncertainty.multitree import fit_over_trees, read_tree_sample, TreeFailureError


TREES = ["((T1:1,T2:1):1,T3:2);", "((T1:1,T3:1):1,T2:2);",
         "(T3:3,(T2:1,T1:1):2);"]
GEOG = {"T1":"10", "T2":"10", "T3":"01"}


def models():
    return [b.model_template(name).with_parameters(d=dict(type="fixed",init=.2),
                e=dict(type="fixed",init=.1),j=dict(type="fixed",init=j))
            for name,j in (("DEC",0),("DEC+J",.1))]


class MultiTreeTests(unittest.TestCase):
    def test_nexus_translate_and_quoted_taxon_names_preserve_order(self):
        text = "#NEXUS\nBegin trees;\nTranslate 1 'T one', 2 T2, 3 T3;\n" \
               "tree STATE_0 = [&R] ((1:1,2:1):1,3:2);\n" \
               "tree STATE_10 = [&R] (3:2,(2:1,1:1):1);\nEnd;"
        trees=read_tree_sample(text)
        self.assertEqual(len(trees),2)
        self.assertEqual(trees[0].tip_labels,("T one","T2","T3"))
        self.assertEqual(trees[1].tip_labels,("T3","T2","T one"))

    def test_rotated_copied_trees_map_nodes_and_preserve_tree_ids(self):
        original="((T1:1,T2:1):1,T3:2);"
        rotated="(T3:2,(T2:1,T1:1):1);"
        result=fit_over_trees([original,rotated,original],GEOG,models(),
            init_strategy="full_cold_start",likelihood_conditioning="tip_survival",num_threads=2)
        fit=b.BioGeoBEARS(original,GEOG,models()[0],likelihood_conditioning="tip_survival").fit()
        self.assertEqual(result.tree_order,(1,2,3))
        np.testing.assert_array_equal(result.clade_support_all_trees,np.ones(5))
        np.testing.assert_allclose(result.prob_failure_bounds_top[0],result.complete_case_probs_top)
        np.testing.assert_allclose(result.cross_tree_prob_sd_top,0,atol=1e-14)
        self.assertEqual(result.n_trees_completed,3)
        self.assertTrue(np.isnan(result.complete_case_probs_bottom[fit.model.tree.root]).all())

    def test_failed_candidate_isolates_entire_tree_and_exact_missing_bounds(self):
        real=b.BioGeoBEARS.fit
        bad=topology_hash(b.read_newick(TREES[1]))
        def injected(model,*args,**kwargs):
            if topology_hash(model.tree)==bad and model.model.name=="DEC+J":
                raise RuntimeError("injected candidate failure")
            return real(model,*args,**kwargs)
        common=dict(init_strategy="full_cold_start",likelihood_conditioning="tip_survival",num_threads=1,
                    on_excess_failure="conservative_bounds")
        with patch.object(b.BioGeoBEARS,"fit",injected),warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            mono=fit_over_trees(TREES,GEOG,models(),**common)
            mrca=fit_over_trees(TREES,GEOG,models(),match_mode="mrca",**common)
        node=node_clades(mono.summary_tree).index(frozenset(("T1","T2")))
        self.assertAlmostEqual(mono.clade_support_all_trees[node],2/3)
        np.testing.assert_allclose(mono.prob_failure_bounds_top[1][node]-mono.prob_failure_bounds_top[0][node],0)
        np.testing.assert_allclose(mrca.prob_failure_bounds_top[1][node]-mrca.prob_failure_bounds_top[0][node],1/3)
        self.assertAlmostEqual(mrca.mean_mrca_excess_tips[node],1/3)
        self.assertTrue(np.isnan(mrca.complete_case_probs_bottom[node]).all())
        failed=mono.per_tree_df[mono.per_tree_df.tree_id==2]
        self.assertTrue((failed.tree_status=="incomplete_model_set").all())
        self.assertTrue(failed.weight.isna().all())
        self.assertEqual(set(failed.fit_status),{"failed","converged"})
        with patch.object(b.BioGeoBEARS,"fit",injected),self.assertRaises(TreeFailureError) as error:
            fit_over_trees(TREES,GEOG,models(),init_strategy="full_cold_start",num_threads=1)
        self.assertEqual(error.exception.result.n_trees_incomplete,1)

    def test_first_fit_exception_is_retried_without_dropping_the_tree(self):
        real=b.BioGeoBEARS.fit
        calls={}
        def transient(model,*args,**kwargs):
            name=model.model.name
            calls[name]=calls.get(name,0)+1
            if calls[name]==1:
                raise RuntimeError("transient optimizer failure")
            return real(model,*args,**kwargs)
        with patch.object(b.BioGeoBEARS,"fit",transient):
            result=fit_over_trees(TREES[:1],GEOG,models(),init_strategy="full_cold_start",num_threads=1)
        self.assertEqual(result.n_trees_completed,1)
        self.assertTrue(result.per_tree_df.retried.all())

    def test_no_complete_tree_returns_unknown_means_and_probability_bounds(self):
        with patch.object(b.BioGeoBEARS,"fit",side_effect=RuntimeError("injected")):
            result=fit_over_trees(TREES[:1],GEOG,models(),init_strategy="full_cold_start",num_threads=1,
                                  max_tree_failure_rate=1,on_excess_failure="conservative_bounds")
        self.assertTrue(np.isnan(result.complete_case_probs_top).all())
        np.testing.assert_array_equal(result.prob_failure_bounds_top[0],np.zeros((5,4)))
        np.testing.assert_array_equal(result.prob_failure_bounds_top[1],np.ones((5,4)))
        self.assertTrue(result.model_winner_freq.isna().all())

    def test_interrupted_checkpoint_restores_full_events_and_node_arrays(self):
        common=dict(init_strategy="full_cold_start",likelihood_conditioning="tip_survival",bsm_maps_per_tree=8,seed=91)
        expected=fit_over_trees(TREES,GEOG,models(),num_threads=1,**common)
        with tempfile.TemporaryDirectory() as temp:
            path=str(Path(temp)/"trees.jsonl")
            append=CheckpointStore.append
            def interrupted(store,record,arrays=None):
                append(store,record,arrays)
                if record["id"]==1:
                    raise KeyboardInterrupt("simulated process interruption")
            with patch.object(CheckpointStore,"append",interrupted),self.assertRaises(KeyboardInterrupt):
                fit_over_trees(TREES,GEOG,models(),num_threads=1,checkpoint_file=path,**common)
            resumed=fit_over_trees(TREES,GEOG,models(),num_threads=2,checkpoint_file=path,resume=True,**common)
            np.testing.assert_array_equal(expected.complete_case_probs_top,resumed.complete_case_probs_top)
            np.testing.assert_array_equal(expected.complete_case_probs_bottom,resumed.complete_case_probs_bottom)
            pd.testing.assert_frame_equal(expected.per_tree_df,resumed.per_tree_df)
            pd.testing.assert_frame_equal(expected.joint_bsm.clado_events_df,resumed.joint_bsm.clado_events_df)
            pd.testing.assert_frame_equal(expected.joint_bsm.ana_events_df,resumed.joint_bsm.ana_events_df)
            with patch.object(b.BioGeoBEARS,"fit",side_effect=AssertionError("unexpected refit")):
                again=fit_over_trees(TREES,GEOG,models(),num_threads=1,checkpoint_file=path,resume=True,**common)
            self.assertEqual(again.joint_bsm.num_maps,24)
            with self.assertRaises(CheckpointMismatchError):
                fit_over_trees(TREES,GEOG,models(),num_threads=1,checkpoint_file=path,resume=True,**{**common,"seed":92})


if __name__=="__main__":
    unittest.main()
