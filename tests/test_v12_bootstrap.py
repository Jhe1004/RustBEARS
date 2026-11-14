"""Actual paired fits, failure envelopes, nested boundaries and lossless resume."""
from dataclasses import replace
from pathlib import Path
import pickle
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import beta

import rustbears as b
from rustbears.uncertainty.bootstrap import failure_probability_intervals, _fit_replicate, _ReplicateTask


def fits(contract="tip_survival",j_lower=1e-5):
    n=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=.2),e=dict(type="fixed",init=.1))
    a=b.model_template("DEC+J").with_parameters(d=dict(type="fixed",init=.2),e=dict(type="fixed",init=.1),
                                             j=dict(init=.05,min=j_lower,max=.5))
    kw=dict(tree="((x:1,y:1):1,z:2);",geog={"x":"10","y":"10","z":"01"},
            likelihood_conditioning=contract,allow_null_tips=contract=="unconditional_complete")
    return b.BioGeoBEARS(model=n,**kw).fit(),b.BioGeoBEARS(model=a,**kw).fit()


class BootstrapTests(unittest.TestCase):
    def test_failed_negative_statistics_keep_both_attempts_without_becoming_observations(self):
        nf,af=fits();number=[0]
        def fake_fit(model,**kwargs):
            if model.model.name=="DEC":
                number[0]+=1
                return nf
            return replace(af,lnL=nf.lnL-(1. if number[0]==1 else .5))
        task=_ReplicateTask(1,"null",nf.model,af.model,nf.model.geog,nf.params,af.params,{},1,(41,42),{})
        with patch.object(b.BioGeoBEARS,"fit",fake_fit):
            row,arrays=_fit_replicate(task)
        self.assertEqual(row["status"],"failed")
        self.assertTrue(np.isnan(row["obj_null"]))
        self.assertTrue(np.isnan(row["obj_alt"]))
        self.assertTrue(np.isnan(row["d_stat"]))
        self.assertEqual(arrays,{})
        self.assertEqual(row["raw_obj_null"],nf.lnL)
        self.assertEqual(row["raw_obj_alt"],nf.lnL-.5)
        self.assertEqual(row["raw_d_stat"],-1.)
        np.testing.assert_allclose([a["raw_d_stat"] for a in row["attempts"]],[-2.,-1.],rtol=0,atol=1e-14)
        self.assertEqual([a["maxiter_multiplier"] for a in row["attempts"]],[1,2])
        self.assertTrue(all(len(a["alternative_fits"])==2 for a in row["attempts"]))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"audit.jsonl"
            b.CheckpointStore(path,{"test":"audit"}).append(row)
            restored=b.CheckpointStore(path,{"test":"audit"},resume=True).records[0]
            self.assertEqual(restored["attempts"],row["attempts"])
            self.assertTrue(np.isnan(restored["d_stat"]))

    def test_recovered_retry_retains_first_failure_raw_scores(self):
        nf,af=fits();number=[0]
        def fake_fit(model,**kwargs):
            if model.model.name=="DEC":
                number[0]+=1
                return nf
            return replace(af,lnL=nf.lnL+(-1. if number[0]==1 else .25))
        task=_ReplicateTask(1,"null",nf.model,af.model,nf.model.geog,nf.params,af.params,{},1,(41,42),{})
        with patch.object(b.BioGeoBEARS,"fit",fake_fit):
            row,_=_fit_replicate(task)
        self.assertEqual(row["status"],"converged")
        self.assertTrue(row["retried"])
        self.assertEqual(row["d_stat"],.5)
        self.assertEqual(row["raw_d_stat"],.5)
        self.assertEqual([a["status"] for a in row["attempts"]],["failed","converged"])
        np.testing.assert_allclose([a["raw_d_stat"] for a in row["attempts"]],[-2.,.5],rtol=0,atol=1e-14)

    def test_cp_failure_envelope_is_not_conservative_count_two_sided_interval(self):
        bounds,ci=failure_probability_intervals(3,2,10)
        self.assertEqual(bounds,(4/11,6/11))
        self.assertAlmostEqual(ci[0],beta.ppf(.025,3,8),places=14)
        self.assertAlmostEqual(ci[1],beta.ppf(.975,6,5),places=14)
        self.assertLess(ci[0],beta.ppf(.025,5,6))
        nofailure,mc=failure_probability_intervals(3,0,10)
        self.assertEqual(nofailure[0],nofailure[1])
        self.assertGreater(mc[1]-mc[0],0)

    def test_actual_refits_and_nested_physical_zero(self):
        null,alt=fits()
        result=b.bootstrap_lrt(null,alt,n_boot=5,num_threads=2,fit_options={"maxiter":30})
        self.assertEqual(result.n_converged+result.n_failed,5)
        self.assertEqual(result.n_failed,0)
        self.assertEqual(result.observed_alt_fit.model.model.params["j"].min,0)
        self.assertEqual(alt.model.model.params["j"].min,1e-5)
        self.assertEqual(result.pval_failure_bounds[0],result.pval_failure_bounds[1])
        self.assertTrue((result.replicates_df.n_evals_null>0).all())
        self.assertTrue((result.replicates_df.n_evals_alt>0).all())
        self.assertTrue((result.replicates_df.d_stat>=0).all())
        self.assertIn("null_ancestor_top",result.replicate_arrays[1])

    def test_custom_positive_bound_is_not_silently_expanded(self):
        null,alt=fits(j_lower=.01)
        with self.assertRaisesRegex(ValueError,"custom positive"):
            b.bootstrap_lrt(null,alt,n_boot=1,num_threads=1)

    def test_failure_is_retained_and_widens_both_intervals(self):
        null,alt=fits()
        def fail_one(task):
            record,arrays=_fit_replicate(task)
            if task.rep_id==2:
                record.update(status="failed",obj_null=np.nan,obj_alt=np.nan,d_stat=np.nan,retried=True,error="injected")
            return record,arrays
        with patch("rustbears.uncertainty.bootstrap._fit_replicate",side_effect=fail_one):
            result=b.bootstrap_lrt(null,alt,n_boot=4,num_threads=1)
        self.assertEqual(result.n_planned,4)
        self.assertEqual(result.n_failed,1)
        self.assertAlmostEqual(result.pval_failure_bounds[1]-result.pval_failure_bounds[0],1/5)
        expected=failure_probability_intervals(result.n_exceedances,1,4)
        self.assertEqual(result.pval_mc_ci95,expected[1])
        self.assertEqual(len(result.replicates_df),4)

    def test_explicit_contract_switch_and_legacy_generation_choice(self):
        null,alt=fits("legacy_unconditioned")
        with self.assertRaises(b.ContractMismatchError):
            b.bootstrap_lrt(null,alt,n_boot=1,num_threads=1)
        with self.assertRaises(b.ContractMismatchError):
            b.bootstrap_lrt(null,alt,n_boot=1,target_conditioning="tip_survival",num_threads=1)
        result=b.bootstrap_lrt(null,alt,n_boot=2,target_conditioning="tip_survival",refit_observed=True,num_threads=1)
        self.assertTrue(result.recomputed_observed_under_new_contract)
        self.assertEqual(result.likelihood_conditioning,"tip_survival")
        self.assertEqual(null.likelihood_conditioning,"legacy_unconditioned")
        raw=b.bootstrap_lrt(null,alt,n_boot=2,simulation_conditioning="tip_survival",num_threads=1)
        self.assertEqual(raw.likelihood_conditioning,"legacy_unconditioned")

    def test_checkpoint_resume_after_interruption_equals_uninterrupted_including_arrays(self):
        null,alt=fits()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"lrt.jsonl"
            def interrupt(task):
                if task.rep_id==3: raise KeyboardInterrupt("deliberate test stop")
                return _fit_replicate(task)
            with patch("rustbears.uncertainty.bootstrap._fit_replicate",side_effect=interrupt),self.assertRaises(KeyboardInterrupt):
                b.bootstrap_lrt(null,alt,n_boot=4,num_threads=1,checkpoint_file=path)
            resumed=b.bootstrap_lrt(null,alt,n_boot=4,num_threads=2,checkpoint_file=path,resume=True)
            full=b.bootstrap_lrt(null,alt,n_boot=4,num_threads=1)
            pd.testing.assert_frame_equal(full.replicates_df,resumed.replicates_df)
            for rep in full.replicate_arrays:
                for key,value in full.replicate_arrays[rep].items():
                    np.testing.assert_array_equal(value,resumed.replicate_arrays[rep][key])
            with self.assertRaises(b.CheckpointMismatchError):
                b.bootstrap_lrt(null,alt,n_boot=4,num_threads=1,checkpoint_file=path,resume=True,seed=999)

    def test_surrogate_power_and_task_serialization(self):
        null,alt=fits()
        task=_ReplicateTask(1,"null",null.model,alt.model,null.model.geog,null.params,alt.params,{},1,(1,2),{})
        self.assertEqual(pickle.loads(pickle.dumps(task)).rep_id,1)
        result=b.bootstrap_lrt(null,alt,n_boot=3,calc_surrogate_power=True,n_alt_sim=3,num_threads=2)
        self.assertIsNotNone(result.surrogate_power_fixed_threshold)
        self.assertEqual(len(result.surrogate_replicates_df),3)
        self.assertGreaterEqual(result.surrogate_power_failure_bounds[1],result.surrogate_power_fixed_threshold)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"lrt.png";fig=result.plot(path)
            self.assertGreater(path.stat().st_size,1000)
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_parallel_budget_rejects_nested_workers(self):
        with self.assertRaises(ValueError): b.ParallelConfig(workers=2,native_threads=2)
        with self.assertRaises(TypeError): b.ParallelConfig(workers=True)


if __name__=="__main__": unittest.main()
