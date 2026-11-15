"""Prospective-budget guards and retained truth/interval/reporting regressions."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import rustbears as b
from scripts import run_uncertainty_validation as validation
from scripts import v12_formal_protocol as protocol
from scripts.run_uncertainty_case_study import fit_posterior_case,write_bsm_outputs


class FormalReportingTests(unittest.TestCase):
    def interval(self,se):
        return SimpleNamespace(se_hessian_local=se,mle=.1,lower=0.,upper=.4,status_lower='physical_zero',
            status_upper='interior_bracketed',interval_kind='asymptotic_confidence_interval',
            critical_value=3.84145882,endpoint_deviance_residual_max=1e-8,onesided_upper_from_zero=.3)

    def test_wald_endpoints_persist_and_unavailable_boundary_is_not_imputed(self):
        valid=validation.interval_record(self.interval(.1),0.,2.)
        self.assertLess(valid['wald_lower'],0.)
        self.assertAlmostEqual(valid['wald_upper'],.2959963984540054)
        self.assertEqual(valid['se_hessian_local'],.1)
        missing=validation.interval_record(self.interval(None),0.,2.)
        self.assertIsNone(missing['wald_covers'])
        self.assertIsNone(missing['wald_upper'])
        self.assertFalse(missing['wald_available'])
        records=[dict(outputs={'asymptotic_chi2_e':valid,'bootstrap_plugin_e':valid}),
                 dict(outputs={'asymptotic_chi2_e':missing})]
        report=validation.summarize_v1(records,3)
        self.assertEqual(len(report),12)
        self.assertEqual(report['wald_local_e']['n_observed'],1)
        self.assertEqual(report['wald_local_e']['n_unavailable'],1)
        self.assertEqual(report['wald_local_e']['n_missing'],2)
        self.assertEqual(report['onesided_zero_reference_e']['n_observed'],2)
        self.assertEqual(report['bootstrap_plugin_e']['n_missing'],2)

    def test_failure_envelope_keeps_planned_denominator(self):
        result=validation.proportion_summary([True,False],4,unavailable=1)
        self.assertEqual(result['missing_outcome_bounds'],[.25,.75])
        self.assertEqual(result['n_observed'],2)
        self.assertLess(result['mc_cp95_envelope'][0],.25)
        self.assertGreater(result['mc_cp95_envelope'][1],.75)

    def test_v2c_uses_true_nodes_preserves_zero_probability_and_clusters_datasets(self):
        rows=[]
        for truth,avg,best in [([0,1],[[.8,.2],[.2,.8]],[[.6,.4],[.6,.4]]),
                               ([0,0],[[1.,0.],[0.,1.]],[[.5,.5],[.5,.5]])]:
            rows.append(dict(outputs={'model_average':validation.score_nodes(np.array(avg),np.array(truth)),
                'best_model':validation.score_nodes(np.array(best),np.array(truth)),
                'weights':[{'weight':.7},{'weight':.3}]}))
        report=validation.summarize_v2c(rows,3)
        self.assertEqual(report['n_missing'],1)
        avg=report['methods']['model_average']
        self.assertAlmostEqual(avg['mean_brier'],.54)
        self.assertEqual(avg['mean_log_probability'],-np.inf)
        self.assertEqual(avg['n_zero_true_probabilities'],1)
        self.assertEqual(avg['n_nonfinite_log_datasets'],1)
        self.assertIsNone(avg['mcse_log_probability'])
        self.assertAlmostEqual(avg['brier_missing_result_bounds'][0],.36)
        self.assertAlmostEqual(avg['brier_missing_result_bounds'][1],3.08/3)
        self.assertEqual(sum(x['n'] for x in avg['calibration']),8)
        self.assertEqual(avg['calibration'][-1]['observed_frequency'],.5)
        self.assertEqual(avg['calibration'][-1]['n_dataset_clusters'],2)
        self.assertAlmostEqual(report['paired_brier_difference_average_minus_best']['mean'],.03)

    def test_v2b_retains_inner_failures_and_labels_distinct_surrogate_target(self):
        row=dict(reject=False,p_failure_bounds=[.02,.08],n_failed=1,asymptotic_p_nominal=.01,
                 surrogate_power_fixed_threshold=.7,surrogate_power_failure_bounds=[.7,.9])
        result=validation.summarize_v2b([dict(outputs={'A':row})],2)['A']
        self.assertEqual(result['rejection_bounds_including_inner_and_outer_failures'],[0.,1.])
        self.assertEqual(result['nominal_chi2_rejection']['n_true'],1)
        self.assertEqual(result['surrogate']['missing_result_bounds'],[.35,.95])
        self.assertIn('Different generating-parameter targets',result['surrogate']['interpretation'])

    def make_budget(self,folder):
        gate=folder/'gate.json';gate.write_text(json.dumps({'gate2_passed':True,'gate3_passed':True}))
        job=dict(analysis='v2b',scenario='alt_j003',n_tips=30,areas=4,outer_replicates=20,
                 inner_bootstrap=199,plugin_bootstrap=99,surrogate_alt_replicates=0,finite_b_extrapolation=False,seed=20261115)
        payload=dict(schema='rustbears.formal_budget.v1',status='locked',locked_utc='2026-01-01T00:00:00+00:00',
            source_sha256={'scripts/run_uncertainty_validation.py':hashlib.sha256(Path(validation.__file__).read_bytes()).hexdigest(),
                           'scripts/v12_formal_protocol.py':hashlib.sha256(Path(protocol.__file__).read_bytes()).hexdigest()},
            gate_evidence={'path':'gate.json','sha256':hashlib.sha256(gate.read_bytes()).hexdigest()},jobs=[job])
        path=folder/'budget.json';path.write_text(json.dumps(payload));return path,payload,gate

    def test_budget_is_prospective_version_bound_and_gate_checked(self):
        with tempfile.TemporaryDirectory() as d:
            path,payload,gate=self.make_budget(Path(d))
            job,planned=validation.read_locked_budget(path,'v2b','alt_j003',30,4)
            self.assertEqual(job['outer_replicates'],20)
            _,executed=validation.read_locked_budget(path,'v2b','alt_j003',30,4,execute=True)
            self.assertEqual(planned,executed)
            gate.write_text('{}')
            with self.assertRaisesRegex(ValueError,'hashed Gate'):
                validation.read_locked_budget(path,'v2b','alt_j003',30,4,execute=True)
            payload['locked_utc']='2999-01-01T00:00:00+00:00';path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'past timezone'):
                validation.read_locked_budget(path,'v2b','alt_j003',30,4)

    def test_budget_distinguishes_paired_B_jobs_and_rejects_source_drift(self):
        with tempfile.TemporaryDirectory() as d:
            path,payload,_=self.make_budget(Path(d))
            payload['jobs'].append({**payload['jobs'][0],'inner_bootstrap':399,'finite_b_extrapolation':True})
            path.write_text(json.dumps(payload))
            job,_=validation.read_locked_budget(path,'v2b','alt_j003',30,4,finite_b_extrapolation=True)
            self.assertEqual(job['inner_bootstrap'],399)
            payload['source_sha256']['scripts/run_uncertainty_validation.py']='wrong';path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'source script'):
                validation.read_locked_budget(path,'v2b','alt_j003',30,4)

    def test_final_budget_cannot_change_frozen_design_or_accept_false_gate(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);path,payload,gate=self.make_budget(folder)
            design=folder/'design.json';design.write_text(json.dumps({'jobs':payload['jobs']}))
            payload['design_evidence']={'path':'design.json','sha256':hashlib.sha256(design.read_bytes()).hexdigest()}
            path.write_text(json.dumps(payload))
            validation.read_locked_budget(path,'v2b','alt_j003',30,4,execute=True)
            gate.write_text(json.dumps({'gate2_passed':True,'gate3_passed':False}))
            payload['gate_evidence']['sha256']=hashlib.sha256(gate.read_bytes()).hexdigest();path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'both pass'):
                validation.read_locked_budget(path,'v2b','alt_j003',30,4,execute=True)
            payload['jobs'][0]['outer_replicates']=30;path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'settings differ'):
                validation.read_locked_budget(path,'v2b','alt_j003',30,4)

    def test_posterior_receives_physical_zero_templates_for_all_families(self):
        with patch('scripts.run_uncertainty_case_study.b.fit_over_trees',return_value='sentinel') as run:
            result=fit_posterior_case('trees','geog',validation.FAMILIES,'summary','tip_survival',2,200,500,Path('/tmp'),2,9)
        self.assertEqual(result,'sentinel')
        models=run.call_args.kwargs['models']
        self.assertEqual([m.name for m in models],list(validation.FAMILIES))
        for model in models:
            for name in ('d','e','j'):
                if model.params[name].type=='free':self.assertEqual(model.params[name].min,0.)

    def test_real_mapping_export_keeps_events_design_and_joint_diagnostics(self):
        model=b.model_template('DEC').with_parameters(d=dict(type='fixed',init=.3),e=dict(type='fixed',init=.2))
        fit=b.BioGeoBEARS('((a:1,b:1):1,(c:1,d:1):1);',dict(a='10',b='01',c='11',d='10'),model,
                         likelihood_conditioning='tip_survival').fit()
        maps=b.model_average([fit],criterion='AIC').run_bsm(num_maps=20,seed=1234,num_threads=1)
        maps.summary_scope='complete_case';maps.n_trees_missing=2
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);result=write_bsm_outputs(maps,folder,'joint_')
            self.assertEqual(result['tree_target'],'fixed_sample')
            self.assertEqual(result['summary_scope'],'complete_case')
            self.assertEqual(result['n_trees_missing'],2)
            metadata=pd.read_csv(folder/'joint_map_metadata.csv')
            self.assertEqual(len(metadata),20)
            for field in ('tree_id','tree_topology_hash','model_id','map_id','stratum_weight','precision_mode','likelihood_conditioning'):
                self.assertIn(field,metadata.columns)
            self.assertEqual(len(pd.read_csv(folder/'joint_anagenetic_events.csv')),len(maps.ana_events_df))
            self.assertTrue((folder/'joint_route_precision.csv').is_file())
            self.assertTrue((folder/'joint_bsm_summary.json').is_file())


if __name__=='__main__':unittest.main()
