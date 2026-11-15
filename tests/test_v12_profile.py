"""Analytical endpoints, physical zeros, checkpoint replay and true nuisance fits."""
from dataclasses import replace
import pickle
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import chi2

import rustbears as b
from rustbears.uncertainty.profile import _PointTask, _Curve
from rustbears.uncertainty.parallel import ParallelConfig


def extinction_fit(contract="unconditional_complete",upper=2.):
    model=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=0),
        e=dict(init=.1,min=0,max=upper))
    return b.BioGeoBEARS("((x:1,y:1):3,z:4);",{"x":"1","y":"1","z":"1"},model,
                        likelihood_conditioning=contract,allow_null_tips=contract=="unconditional_complete").fit()


class ProfileTests(unittest.TestCase):
    def test_rebase_changes_cache_generation_without_changing_parameter_box(self):
        fit=extinction_fit()
        other=replace(fit,params={**fit.params,"e":.1},lnL=-.9)
        curves=[_Curve(f,"e",{},ParallelConfig(),map,{},lambda _:None) for f in (fit,other)]
        self.assertNotEqual(curves[0].generation,curves[1].generation)

    def test_psychotria_bayarealike_j_profile_rebase_does_not_reuse_stale_basin(self):
        # This real 19-tip fixture formerly joined an e≈0 cached grid ordinate
        # to e=5 root-search ordinates after rebase, creating a false jump and
        # an endpoint residual of 0.0557264. The parameter box did not change.
        from rustbears.uncertainty.profile import _profile_point
        root=Path(__file__).resolve().parent/"fixtures/psychotria"
        model=b.BioGeoBEARS(root/"Psychotria_5.2_ultrametric.newick",root/"Psychotria_geog.data",
            "BAYAREALIKE+J",likelihood_conditioning="tip_survival",max_range_size=4)
        # Replay the recorded inferior D1 basin. An unconstrained initial fit
        # can already find the better basin on some platforms, in which case
        # asserting a later rebase would depend on the optimizer's trajectory.
        # The score and ancestral probabilities are recomputed by the native
        # likelihood; the profile and endpoint check below remain real fits.
        params=model.model.resolve(dict(d=1e-12,e=1e-12,j=.10811626692419736))
        evaluation=model.evaluate(params)
        self.assertAlmostEqual(evaluation.lnL,-24.260680130280434,places=8)
        k=len(model.model.free);n=model.tree.tip_count
        aic=2*k-2*evaluation.lnL
        fit=b.FitResult(lnL=evaluation.lnL,params=params,n_params=k,aic=aic,
            aicc=aic+2*k*(k+1)/(n-k-1),n_evals=1,converged=True,
            evaluation=evaluation,model=model,message="Replay recorded D1 local optimum")
        result=b.profile_ci(fit,params=["d"],n_grid=8,num_threads=2)
        self.assertTrue(result.rebased_mle)
        self.assertGreater(result.obj_mle,fit.lnL+.5)
        interval=result.intervals["d"]
        self.assertLess(interval.endpoint_deviance_residual_max,1e-4)
        self.assertIsNotNone(interval.upper)
        point=_profile_point(_PointTask(result.updated_fit.model,"d",interval.upper,
                                        result.updated_fit.params,{},1))
        self.assertLess(abs(2*(result.obj_mle-point["profile_obj"])-interval.critical_value),1e-4)
        self.assertEqual(interval.profile_curve_df.generation.nunique(),1)

    def test_analytical_extinction_endpoints_and_one_sided_reference(self):
        result=b.profile_ci(extinction_fit(),params=["e"],n_grid=5,include_onesided_zero_ref=True,num_threads=2)
        interval=result.intervals["e"]
        self.assertEqual(interval.lower,0.)
        self.assertEqual(interval.status_lower,"at_physical_zero")
        self.assertAlmostEqual(interval.upper,chi2.ppf(.95,1)/18,places=12)
        self.assertAlmostEqual(interval.onesided_upper_from_zero,chi2.ppf(.9,1)/18,places=12)
        self.assertLess(interval.endpoint_deviance_residual_max,1e-10)
        self.assertEqual(interval.interval_kind,"asymptotic_confidence_interval")

    def test_box_censoring_and_flat_conditional_likelihood_are_explicit(self):
        p=b.profile_ci(extinction_fit(upper=.11),params=["e"],n_grid=4).intervals["e"]
        self.assertIsNone(p.upper)
        self.assertEqual(p.status_upper,"censored_at_box_upper")
        flat=b.profile_ci(extinction_fit("tip_survival"),params=["e"],n_grid=4).intervals["e"]
        self.assertIsNone(flat.upper)
        self.assertEqual(flat.lower,0)

    def test_legacy_always_labels_support_interval(self):
        result=b.profile_ci(extinction_fit("legacy_unconditioned"),params=["e"],n_grid=4)
        self.assertEqual(result.interval_kind,"profile_support_interval")
        self.assertEqual(result.intervals["e"].interval_kind,"profile_support_interval")

    def test_custom_lower_box_is_respected_and_correctly_labelled(self):
        from rustbears.uncertainty.guards import clone_model
        original=extinction_fit()
        obj=original.model.model.with_parameters(e=dict(init=.1,min=.1,max=2))
        fit=clone_model(original.model,model_object=obj).fit()
        interval=b.profile_ci(fit,params=["e"],eval_exact_zero=False,n_grid=4).intervals["e"]
        self.assertIsNone(interval.lower)
        self.assertEqual(interval.status_lower,"censored_at_box_lower")
        with self.assertRaisesRegex(b.ProfileOptimizationError,"custom positive"):
            b.profile_ci(fit,params=["e"],n_grid=4)

    def test_structurally_impossible_zero_has_positive_certified_endpoint(self):
        model=b.model_template("DEC").with_parameters(d=dict(init=.2,min=0,max=5),e=dict(type="fixed",init=.1))
        fit=b.BioGeoBEARS("((x:1,y:1):1,z:2);",{"x":"10","y":"10","z":"01"},model,
            likelihood_conditioning="tip_survival",root_prior=b.RootPriorSpec("custom_weights",{"A":1})).fit()
        interval=b.profile_ci(fit,params=["d"],n_grid=5).intervals["d"]
        self.assertEqual(interval.status_lower,"infeasible_zero_approached")
        self.assertGreater(interval.lower,0)
        self.assertLess(interval.endpoint_deviance_residual_max,1e-4)
        self.assertTrue(np.isneginf(interval.profile_curve_df.loc[interval.profile_curve_df.value==0,"profile_obj"]).all())

    def test_discovered_better_profile_rebases_and_returns_updated_fit(self):
        fit=extinction_fit()
        bad=replace(fit,lnL=-.9,params={**fit.params,"e":.1},evaluation=fit.model.evaluate({"e":.1}))
        result=b.profile_ci(bad,params=["e"],n_grid=4)
        self.assertTrue(result.rebased_mle)
        self.assertAlmostEqual(result.obj_mle,0.)
        self.assertEqual(result.updated_fit.params["e"],0.)
        self.assertEqual(bad.params["e"],.1)
        with self.assertRaises(b.ProfileOptimizationError):
            b.profile_ci(bad,params=["e"],n_grid=4,multi_basin_guard=False)

    def test_plugin_uses_simulated_profile_drop_and_labels_local_approximation(self):
        result=b.profile_ci(extinction_fit(),params=["e"],method="bootstrap_plugin",n_boot_plugin=4,n_grid=4,num_threads=2)
        interval=result.intervals["e"]
        self.assertEqual(interval.interval_kind,"bootstrap_plugin_interval")
        self.assertEqual(len(interval.bootstrap_plugin_df),4)
        self.assertEqual(interval.critical_value,0.)
        self.assertEqual(interval.upper,0.)
        with self.assertRaises(b.ContractMismatchError):
            b.profile_ci(extinction_fit("legacy_unconditioned"),params=["e"],method="bootstrap_plugin",n_boot_plugin=2,n_grid=4)

    def test_profile_checkpoint_replays_all_points_without_fitting(self):
        fit=extinction_fit()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"profile.jsonl"
            first=b.profile_ci(fit,params=["e"],n_grid=5,checkpoint_file=path)
            with patch("rustbears.uncertainty.profile._profile_point",side_effect=AssertionError("unexpected recomputation")):
                replay=b.profile_ci(fit,params=["e"],n_grid=5,checkpoint_file=path,resume=True,num_threads=2)
            pd.testing.assert_frame_equal(first.intervals["e"].profile_curve_df,replay.intervals["e"].profile_curve_df)
            self.assertEqual(first.total_evals,replay.total_evals)
            with self.assertRaises(b.CheckpointMismatchError):
                b.profile_ci(fit,params=["e"],n_grid=6,checkpoint_file=path,resume=True)

    def test_task_is_serializable_and_plot_is_real(self):
        fit=extinction_fit()
        task=_PointTask(fit.model,"e",.1,fit.params,{},1)
        restored=pickle.loads(pickle.dumps(task))
        self.assertEqual(restored.param,"e")
        result=b.profile_ci(fit,params=["e"],n_grid=4)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"profile.png"
            fig=result.plot(path)
            self.assertGreater(path.stat().st_size,1000)
            import matplotlib.pyplot as plt
            plt.close(fig)

    def test_surface_refits_nuisance_and_preserves_grid(self):
        model=b.model_template("DEC+J").with_parameters(d=dict(init=.2,min=.01,max=.5),
            e=dict(type="fixed",init=.1),j=dict(init=.1,min=0,max=.5))
        fit=b.BioGeoBEARS("((x:1,y:1):1,z:2);",{"x":"10","y":"10","z":"01"},model,
                        likelihood_conditioning="tip_survival").fit()
        surface=b.profile_surface(fit,"d","j",grid_x=[.1,.2],grid_y=[0,.1],num_threads=2)
        self.assertEqual(len(surface),4)
        for row in surface.itertuples():
            self.assertAlmostEqual(row.profile_obj,fit.model.evaluate({"d":row.d,"j":row.j}).lnL,places=12)


if __name__=="__main__": unittest.main()
