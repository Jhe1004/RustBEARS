"""Independent checks of probability/refit/recovery boundaries (small real fits)."""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import rustbears as b
from rustbears.uncertainty.bootstrap import _fit_replicate, _ReplicateTask
from rustbears.uncertainty.checkpoint import CheckpointStore, model_manifest
from rustbears.uncertainty.profile import _Curve
from rustbears.uncertainty.parallel import ParallelConfig


def extinction_model():
    obj=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=0),
        e=dict(init=.001,min=0,max=2))
    return b.BioGeoBEARS("((x:1,y:1):1,z:2);",{"x":"1","y":"0","z":"1"},obj,
        likelihood_conditioning="unconditional_complete",allow_null_tips=True)


def nested_fits():
    n=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=.2),e=dict(type="fixed",init=.1))
    a=b.model_template("DEC+J").with_parameters(d=dict(type="fixed",init=.2),e=dict(type="fixed",init=.1),
        j=dict(init=.05,min=0,max=.5))
    kw=dict(tree="((x:1,y:1):1,z:2);",geog={"x":"10","y":"10","z":"01"},
        likelihood_conditioning="tip_survival")
    return b.BioGeoBEARS(model=n,**kw).fit(),b.BioGeoBEARS(model=a,**kw).fit()


class IndependentInteractionTests(unittest.TestCase):
    def test_simulation_uses_actual_fitted_parameters_and_replicate_ids(self):
        model=extinction_model();fit=model.fit()
        self.assertGreater(fit.params["e"],model.model.params["e"].init*10)
        fitted=b.simulate_datasets(fit,n_reps=128,seed=723,num_threads=2)
        explicit=b.simulate_datasets(model,params=fit.params,n_reps=128,seed=723)
        initial=b.simulate_datasets(model,n_reps=128,seed=723)
        np.testing.assert_array_equal(fitted.node_states,explicit.node_states)
        self.assertFalse(np.array_equal(fitted.node_states,initial.node_states))
        with b.PreparedSimulation(fit) as sampler:
            first=sampler.draw(47,seed=723)
            last=sampler.draw(81,start_rep=48,seed=723,num_threads=2)
        np.testing.assert_array_equal(fitted.node_states,np.concatenate((first.node_states,last.node_states)))
        self.assertTrue((fitted.tip_states==0).any())
        with model.prepared_likelihood(fit.params) as prepared:
            self.assertTrue(np.isfinite(prepared.evaluate_tip_states(fitted.tip_states)).all())

    def test_actual_refit_produces_a_separate_profile_cache_generation(self):
        model=extinction_model();fit=model.fit()
        initial=replace(fit,params=model.model.resolve(),lnL=model.evaluate().lnL,
                        evaluation=model.evaluate())
        cache={};committed=[]
        first=_Curve(initial,"e",{},ParallelConfig(),map,cache,committed.append)
        after=_Curve(fit,"e",{},ParallelConfig(),map,cache,committed.append)
        self.assertNotEqual(first.generation,after.generation)
        first.point(.2);after.point(.2)
        self.assertEqual(len(committed),2)
        self.assertEqual({row["generation"] for row in committed},{first.generation,after.generation})

    def test_negative_nested_geometry_is_retried_and_kept_as_failure(self):
        nf,af=nested_fits();bad=replace(af,lnL=nf.lnL-1)
        # All optimizers claim convergence, but the alternative violates the
        # exact nesting inequality; it must not become a truncated D=0 success.
        def result(model,**kwargs):
            return nf if model.model.name=="DEC" else bad
        task=_ReplicateTask(2,"null",nf.model,af.model,nf.model.geog,nf.params,af.params,{},1,(41,42),
            {"test_identity":"must_survive"})
        with patch.object(b.BioGeoBEARS,"fit",result):
            row,arrays=_fit_replicate(task)
        self.assertEqual(row["id"],2)
        self.assertEqual(row["status"],"failed")
        self.assertTrue(row["retried"])
        self.assertEqual(row["error"].count("negative nested statistic"),2)
        self.assertEqual(row["rng"],task.rng)
        self.assertEqual(row["params_null"],nf.params)
        self.assertEqual(row["params_alt"],af.params)
        self.assertTrue(np.isnan(row["d_stat"]))
        self.assertEqual(arrays,{})

    def test_failed_replicates_and_probability_arrays_survive_full_resume(self):
        nf,af=nested_fits()
        def failed(task):
            row,arrays=_fit_replicate(task)
            if task.rep_id==2:
                row.update(status="failed",retried=True,d_stat=np.nan,error="independent injected failure")
            return row,arrays
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"bootstrap.jsonl"
            with patch("rustbears.uncertainty.bootstrap._fit_replicate",side_effect=failed):
                first=b.bootstrap_lrt(nf,af,n_boot=3,num_threads=1,checkpoint_file=path,seed=310)
            with patch("rustbears.uncertainty.bootstrap._fit_replicate",side_effect=AssertionError("completed replicate recomputed")):
                resumed=b.bootstrap_lrt(nf,af,n_boot=3,num_threads=2,checkpoint_file=path,seed=310,resume=True)
            pd.testing.assert_frame_equal(first.replicates_df,resumed.replicates_df)
            self.assertEqual(resumed.n_failed,1)
            self.assertEqual(resumed.pval_failure_bounds,first.pval_failure_bounds)
            self.assertAlmostEqual(resumed.pval_failure_bounds[1]-resumed.pval_failure_bounds[0],1/4)
            for rid,values in first.replicate_arrays.items():
                for key,value in values.items():np.testing.assert_array_equal(value,resumed.replicate_arrays[rid][key])

    def test_fitted_reference_and_raw_root_scale_are_fingerprinted(self):
        fit=extinction_model().fit()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"inputs.jsonl"
            CheckpointStore(path,model_manifest(fit))
            changed=replace(fit,params={**fit.params,"e":fit.params["e"]+.01})
            with self.assertRaises(b.CheckpointMismatchError):
                CheckpointStore(path,model_manifest(changed),resume=True)
        kw=dict(tree="(x:1,y:1);",geog={"x":"10","y":"01"},likelihood_conditioning="legacy_unconditioned")
        a=b.BioGeoBEARS(**kw,root_prior=b.RootPriorSpec("custom_weights",{"A":1,"B":2,"AB":3}))
        c=b.BioGeoBEARS(**kw,root_prior=b.RootPriorSpec("custom_weights",{"A":2,"B":4,"AB":6}))
        np.testing.assert_array_equal(a.root_sampling_probs,c.root_sampling_probs)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"root.jsonl";CheckpointStore(path,model_manifest(a))
            with self.assertRaises(b.CheckpointMismatchError):CheckpointStore(path,model_manifest(c),resume=True)

    def test_shared_prepared_sampler_serializes_concurrent_draws(self):
        fit=extinction_model().fit()
        with b.PreparedSimulation(fit) as sampler:
            serial=[sampler.draw(23,start_rep=s,seed=44).node_states for s in (1,24)]
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures=[pool.submit(sampler.draw,23,start_rep=s,seed=44) for s in (1,24)]
                concurrent=[f.result().node_states for f in futures]
        for expected,actual in zip(serial,concurrent):np.testing.assert_array_equal(expected,actual)
        with self.assertRaisesRegex(RuntimeError,"closed"):sampler.draw(1)


if __name__=="__main__":unittest.main()
