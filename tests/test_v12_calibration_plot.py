"""Nonvisual checks of Figure 2 provenance, missingness and fixed selection.

Fixtures are arithmetic examples only. These tests never render or create
manuscript figures, and their budgets are temporary test data, not authorization.
"""
import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import plot_v12_calibration as plot
from rustbears.uncertainty.checkpoint import canonical_manifest,task_fingerprint


class CalibrationPlotTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)

    def checkpoint(self,manifest,records,*,truncate=False):
        header=dict(schema='rustbears.checkpoint.v1',task_fingerprint=task_fingerprint(manifest),manifest=canonical_manifest(manifest))
        path=self.path/'records.jsonl'
        text='\n'.join(json.dumps(x) for x in [header,*[dict(record=r) for r in records]])+'\n'
        path.write_text(text[:-1] if truncate else text)
        return path,header

    def budget(self):
        gate=self.path/'gates.json';plot.write_json(gate,dict(gate2_passed=True,gate3_passed=True))
        jobs=[dict(analysis=a,scenario='test',n_tips=4,areas=2) for a in plot.ANALYSES]
        budget=dict(schema='rustbears.formal_budget.v1',status='locked',locked_utc='2020-01-01T00:00:00+00:00',
                    gate_evidence=dict(path=gate.name,sha256=plot.sha(gate)),jobs=jobs)
        path=self.path/'budget.json';plot.write_json(path,budget)
        inv=dict(schema='rustbears.figure2.inputs.v1',formal=True,budget=dict(path=path.name,sha256=plot.sha(path)),runs=[])
        return budget,path,inv

    def test_default_is_contract_only(self):
        output=self.path/'prepared'
        result=subprocess.run([sys.executable,str(ROOT/'scripts/plot_v12_calibration.py'),'--output',str(output)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual({p.suffix for p in output.iterdir()},{'.json'})
        self.assertFalse(plot.read_json(output/'input_audit.json')['rendered'])
        result=subprocess.run([sys.executable,str(ROOT/'scripts/plot_v12_calibration.py'),'--output',str(output),'--render'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(any(p.suffix in ('.pdf','.svg','.png','.tiff') for p in output.rglob('*')))

    def test_direct_render_requires_verified_eligibility(self):
        with self.assertRaisesRegex(ValueError,'eligibility'):
            plot.render({},self.path,self.path)
        with self.assertRaisesRegex(ValueError,'changed'):
            plot.render({},self.path,self.path,eligibility=dict(status='eligible_formal_inputs',formal=True,source_files={}))

    def test_project_glyph_floor_is_eight_points(self):
        self.assertGreaterEqual(plot.CONTRACT['minimum_glyph_pt'],8)
        self.assertGreaterEqual(plot.CONTRACT['primary_font_pt'],8)

    def test_recorded_workflow_digest_must_match_frozen_budget(self):
        budget=dict(source_sha256={'scripts/run_uncertainty_validation.py':'correct'})
        self.assertIn('matches',plot.verify_workflow_source(dict(script_sha256='correct'),'v1',budget))
        for digest in (None,'changed'):
            with self.assertRaisesRegex(ValueError,'workflow source'):
                plot.verify_workflow_source(dict(script_sha256=digest),'v2c',budget)
        self.assertIn('no separate',plot.verify_workflow_source({},'v3',budget))

    def test_outer_formal_flag_cannot_promote_underlying_smoke(self):
        manifest=dict(mode='smoke')
        plot.write_json(self.path/'preregistered_manifest.json',dict(manifest=canonical_manifest(manifest)))
        plot.write_json(self.path/'summary.json',dict(mode='smoke'))
        with self.assertRaisesRegex(ValueError,'smoke/pilot'):
            plot.load_run(dict(analysis='v1',path='.'),self.path,'arbitrary')

    def test_budget_requires_actual_passed_gates_and_hashes(self):
        budget,path,inv=self.budget()
        self.assertEqual(len(plot.budget_info(inv,self.path)[2]),5)
        inv['formal']=False
        with self.assertRaisesRegex(ValueError,'formal=true'):plot.budget_info(inv,self.path)
        inv['formal']=True
        (self.path/'gates.json').write_text('{"gate2_passed":true,"gate3_passed":false}')
        with self.assertRaisesRegex(ValueError,'SHA'):plot.budget_info(inv,self.path)
        budget['gate_evidence']['sha256']=plot.sha(self.path/'gates.json');plot.write_json(path,budget);inv['budget']['sha256']=plot.sha(path)
        with self.assertRaisesRegex(ValueError,'Gate 2 and Gate 3'):plot.budget_info(inv,self.path)

    def test_budget_rejects_changed_prospective_design(self):
        budget,path,inv=self.budget();design=self.path/'design.json'
        plot.write_json(design,dict(jobs=[]))
        budget['design_evidence']=dict(path=design.name,sha256=plot.sha(design));plot.write_json(path,budget)
        inv['budget']['sha256']=plot.sha(path)
        with self.assertRaisesRegex(ValueError,'frozen design'):plot.budget_info(inv,self.path)

    def test_every_planned_checkpoint_outcome_retained_including_failure(self):
        manifest=dict(analysis='v1',mode='formal',outer=2,precision='strict_f64')
        path,header=self.checkpoint(manifest,[dict(id=1,status='complete'),dict(id=2,status='failed')])
        rows,_=plot.checkpoint(path,2,expected_fingerprint=header['task_fingerprint'],expected_manifest=header['manifest'])
        self.assertEqual(rows[1]['status'],'failed')
        for records in ([dict(id=1)],[dict(id=1),dict(id=1)]):
            path,_=self.checkpoint(manifest,records)
            with self.assertRaisesRegex(ValueError,'every planned'):plot.checkpoint(path,2)
        path,_=self.checkpoint(manifest,[dict(id=1)],truncate=True)
        with self.assertRaisesRegex(ValueError,'truncated'):plot.checkpoint(path,1)

    def test_checkpoint_rejects_payload_and_configuration_tampering(self):
        path,header=self.checkpoint(dict(a=1),[dict(id=1)])
        changed=copy.deepcopy(header);changed['manifest']=canonical_manifest(dict(a=2))
        path.write_text(json.dumps(changed)+'\n'+json.dumps(dict(record=dict(id=1)))+'\n')
        with self.assertRaisesRegex(ValueError,'digest'):plot.checkpoint(path,1)
        path,header=self.checkpoint(dict(a=1),[dict(id=1)])
        with self.assertRaisesRegex(ValueError,'complete configuration'):
            plot.checkpoint(path,1,expected_manifest=canonical_manifest(dict(a=2)))

    def test_nonfinite_json_is_not_silently_rendered(self):
        path=self.path/'invalid.json';path.write_text('{"x":NaN}')
        with self.assertRaisesRegex(ValueError,'non-standard'):plot.read_json(path)

    def test_rate_missingness_and_binomial_interval(self):
        row=plot.rate_row([True]*19+[False],20)
        self.assertAlmostEqual(row['estimate'],.95)
        self.assertAlmostEqual(row['cp95_available'][0],.7512672372279722)
        self.assertAlmostEqual(row['cp95_available'][1],.9987349105020502)
        self.assertEqual(row['failure_bounds'],[.95,.95])
        absent=plot.rate_row([None,None],2)
        self.assertIsNone(absent['estimate']);self.assertIsNone(absent['cp95_available'])
        self.assertEqual(absent['failure_bounds'],[0,1])
        self.assertEqual(absent['cp95_planned_envelope'],[0,1])
        self.assertEqual(absent['status'],'unavailable')
        with self.assertRaisesRegex(ValueError,'booleans'):plot.rate_row(['False'],1)

    def fixture_runs(self):
        # Numerical unit fixtures; no files are eligible for the rendering gate.
        runs=[]
        def add(analysis,job,rows,summary,manifest=None):
            runs.append(dict(analysis=analysis,job=dict(analysis=analysis,**job),records=rows,
                             summary=summary,manifest=manifest or {},path=self.path/str(len(runs))))
        for tips in (40,150):
            for scenario in plot.V1_SCENARIOS:
                outputs={f'{method}_{parameter}':dict(covered=True,wald_covers=None,one_sided_reference_covers=False,
                    lower=0.,upper=1.,wald_lower=None,wald_upper=None) for method in ('asymptotic_chi2','bootstrap_plugin') for parameter in 'dej'}
                add('v1',dict(n_tips=tips,areas=4,scenario=scenario,outer_replicates=2),
                    [dict(id=1,outputs=outputs),dict(id=2,outputs={})],dict(negative_result=False),dict(truth=dict(j=.03)))
        for scenario in plot.V2_SCENARIOS:
            out={r:dict(p_failure_bounds=[.01,.2],n_failed=1,reject=False,finite_b_subsampling={'99':.4}) for r in plot.ROUTES}
            add('v2b',dict(scenario=scenario,finite_b_extrapolation=False,outer_replicates=2,inner_bootstrap=99),
                [dict(id=1,outputs=out),dict(id=2,outputs={})],dict(finite_b={'all_values':[0.,1.]}))
        support={k:dict(n_datasets_with_scored_clades=0,n_clades=0,paired_dataset_mean_brier_difference=None,
                  paired_dataset_bootstrap_ci95=None,dataset_scores=[]) for k in ('all','<0.7','0.7-0.95','>0.95')}
        add('v3',{},[dict(dataset=1,status='failed_retained')],dict(n_planned=1,n_missing=1,support_bins=support))
        rows=[dict(id=i+1,section=section,estimates={event:dict(mean=1.,se=.1) for event in 'de'}) for i,section in enumerate(('fixed','adaptive'))]
        add('v4',dict(settings={s:dict(replicates=2) for s in ('fixed','adaptive')}),rows,
            dict(summaries={s:{e:dict(truth=1.) for e in 'de'} for s in ('fixed','adaptive')}))
        return runs

    def test_predeclared_panels_keep_every_scenario_and_missing_method(self):
        tables=plot.main_tables(self.fixture_runs())
        self.assertEqual([len(tables[k]) for k in 'abcd'],[32,15,3,4])
        self.assertEqual({(r['n_tips'],r['scenario'],r['method']) for r in tables['a']},
            {(n,s,m) for n in (40,150) for s in plot.V1_SCENARIOS for m in plot.METHODS})
        self.assertTrue(all(r['estimate'] is None for r in tables['a'] if r['method']=='wald_local'))
        self.assertTrue(all(r['failure_bounds']==[0,1] for r in tables['b']))
        self.assertTrue(all(r['estimate'] is None for r in tables['c']))
        self.assertTrue(all(r['n_missing']==1 for r in tables['d']))

    def test_nonfinite_v4_estimate_stays_unavailable(self):
        runs=self.fixture_runs();runs[-1]['records'][0]['estimates']['d']['se']='nan'
        fixed=next(r for r in plot.main_tables(runs)['d'] if r['section']=='fixed' and r['event']=='d')
        self.assertIsNone(fixed['estimate']);self.assertEqual(fixed['n_missing'],2)

    def test_export_retains_secondary_results_and_flags(self):
        runs=self.fixture_runs();tables=plot.main_tables(runs)
        artifacts=plot.export_tables(runs,tables,self.path)
        folder=self.path/'source_data'
        with (folder/'v1_all_parameter_intervals.csv').open() as handle:rows=list(csv.DictReader(handle))
        self.assertEqual(len(rows),8*2*3)
        self.assertEqual({r['interval'].rsplit('_',1)[1] for r in rows},set('dej'))
        with (folder/'all_summary_values.csv').open() as handle:leaves=list(csv.DictReader(handle))
        self.assertTrue(any(r['value']=='False' for r in leaves))
        self.assertTrue(any(r['metric_path']=='/finite_b/all_values/0' and r['value']=='0.0' for r in leaves))
        self.assertIn('source_data/all_run_summaries.jsonl',artifacts)
        self.assertFalse(any(p.suffix in ('.pdf','.svg','.png','.tiff') for p in self.path.rglob('*')))


if __name__=='__main__':unittest.main()
