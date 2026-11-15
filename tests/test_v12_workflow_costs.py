"""Arithmetic/provenance fixtures only: never render synthetic manuscript data."""
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import v12_workflow_costs as costs
from scripts import plot_v12_performance as plot


class WorkflowCostTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.job=dict(id='test',analysis='v2c',reserved_cpus=2)
        self.manifest=dict(n_tips=30,areas=4,scenario='higher_jump',outer_replicates=2,
                           resolved_input=dict(states=list(range(16))))
        self.start=datetime(2020,1,1,tzinfo=timezone.utc)

    def row(self,**changes):
        row=dict(mode='formal',analysis='v2c',scenario='higher_jump',workers=2,native_threads=1,
            n_tips=30,areas=4,states=16,telemetry_source_sha256='test-helper',status='returned',phase='fit',
            wall_seconds=1.,process_cpu_seconds=.5,checkpoint_rows_before=0,checkpoint_rows_after=0,
            started_utc='2020-01-02T00:00:00+00:00',host='test-host',platform='test',pid=1,
            invocation='test-invocation',thread_environment={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},
            replicate=1,model='DEC',contract='tip_survival',call_workers=1,converged=True,resume_requested=True)
        row.update(changes);return row

    def read(self,rows,job=None,manifest=None):
        path=self.root/'timings.jsonl';path.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        return costs.read_timings(path,job or self.job,manifest or self.manifest,'test-helper',self.start)

    def test_default_writes_contract_without_images(self):
        result=subprocess.run([sys.executable,str(ROOT/'scripts/plot_v12_performance.py'),'--output',str(self.root)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual({p.suffix for p in self.root.iterdir()},{'.json'})
        self.assertGreaterEqual(plot.CONTRACT['minimum_font_pt'],8)
        self.assertIn('64/256',plot.CONTRACT['scope'])
        self.assertIn('supplementary only',plot.CONTRACT['selection'])

    def test_p0_cannot_replace_formal_workflow_inventory(self):
        result=subprocess.run([sys.executable,str(ROOT/'scripts/plot_v12_performance.py'),'--output',str(self.root),
            '--render','--audit-scripts',str(self.root),'--benchmark','a=unused','--benchmark','b=unused','--p0','a=unused'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('P0 cannot substitute',result.stderr)
        self.assertFalse(any(p.suffix in ('.pdf','.svg','.png') for p in self.root.iterdir()))

    def test_full_first_call_with_resume_flag_is_fresh(self):
        row=self.read([self.row()])[0]
        self.assertTrue(row['resume_requested']);self.assertTrue(row['api_complete_fresh'])

    def test_existing_checkpoint_and_raised_costs_are_retained(self):
        rows=self.read([self.row(checkpoint_rows_before=2,checkpoint_rows_after=4),
            self.row(invocation='second',status='raised',error='retained',checkpoint_rows_after=1)])
        self.assertEqual([r['classification'] for r in rows],['partial_resume','raised'])
        self.assertFalse(any(r['api_complete_fresh'] for r in rows))
        self.assertEqual(len(rows),2)

    def test_smoke_wrong_source_dimensions_workers_and_thread_budget_rejected(self):
        for change in [dict(mode='smoke'),dict(telemetry_source_sha256='wrong'),dict(n_tips=31),dict(states=32),
                       dict(workers=8),dict(native_threads=2),dict(thread_environment={})]:
            with self.subTest(change=change),self.assertRaises(ValueError):self.read([self.row(**change)])

    def test_negative_nonfinite_duplicate_and_truncated_timings_rejected(self):
        for change in [dict(wall_seconds=-1),dict(wall_seconds=float('nan')),dict(wall_seconds=0),
                       dict(checkpoint_rows_before=3,checkpoint_rows_after=2),dict(started_utc='2019-01-01T00:00:00Z')]:
            with self.subTest(change=change),self.assertRaises(ValueError):self.read([self.row(**change)])
        with self.assertRaisesRegex(ValueError,'duplicated'):self.read([self.row(),self.row()])
        path=self.root/'partial';path.write_text(json.dumps(self.row()))
        with self.assertRaisesRegex(ValueError,'truncated'):costs.read_timings(path,self.job,self.manifest,'test-helper',self.start)

    def six(self):
        rows=[self.row(model=model,started_utc=f'2020-01-02T00:00:{2*i:02d}+00:00') for i,model in enumerate(costs.MODELS)]
        rows.append(self.row(phase='model_average',model=None,n_models=6,started_utc='2020-01-02T00:00:12+00:00'))
        return rows

    def test_six_fit_cost_sums_only_same_invocation_replicate(self):
        rows=self.read(self.six());points=costs.cost_points(rows)
        self.assertEqual(len(points),1);self.assertEqual(points[0]['wall_seconds'],7)
        self.assertEqual(points[0]['component_count'],7)
        for field,value in [('replicate',2),('invocation','another')]:
            mixed=self.six();mixed[-1][field]=value
            self.assertEqual(costs.cost_points(self.read(mixed)),[])

    def test_partial_failed_fit_or_missing_model_does_not_become_six_fit_cost(self):
        for change in [dict(status='raised'),dict(converged=False),dict(checkpoint_rows_before=1,checkpoint_rows_after=1)]:
            rows=self.six();rows[2].update(change)
            self.assertEqual(costs.cost_points(self.read(rows)),[])
        self.assertEqual(costs.cost_points(self.read(self.six()[:-1])),[])

    def test_overlapping_api_durations_cannot_be_summed(self):
        rows=self.six();rows[0]['wall_seconds']=3
        with self.assertRaisesRegex(ValueError,'overlapping'):self.read(rows)

    def test_all_fresh_attempts_retained_not_fastest_selected(self):
        rows=self.read(self.six()+[dict(r,invocation='retry',wall_seconds=.5) for r in self.six()])
        points=costs.cost_points(rows)
        self.assertEqual([r['wall_seconds'] for r in points],[7,3.5])
        self.assertEqual(len({r['task_id'] for r in points}),1)

    def test_bootstrap_includes_surrogate_and_retained_inner_failures(self):
        job=dict(self.job,analysis='v2b');manifest=dict(self.manifest,inner_bootstrap=199,surrogate_alt_replicates=199)
        row=self.row(analysis='v2b',phase='bootstrap_lrt',n_boot=199,n_alt_sim=199,n_planned=199,n_failed=3,
                     route='C_to_A',checkpoint_rows_after=398)
        rows=self.read([row],job,manifest);point=costs.cost_points(rows)[0]
        self.assertEqual(point['n_failed'],3);self.assertEqual(point['workload'],'lrt_199_199')
        with self.assertRaisesRegex(ValueError,'incomplete'):self.read([dict(row,checkpoint_rows_after=199)],job,manifest)
        # A complete B399 call remains in raw supplementary rows, not B199 main points.
        other=dict(row,n_boot=399,n_planned=399,checkpoint_rows_after=598)
        self.assertEqual(costs.cost_points(self.read([other],job,dict(manifest,inner_bootstrap=399))),[])

    def test_planned_denominators_do_not_depend_on_observed_success(self):
        job=dict(self.job,analysis='v1')
        planned=costs.planned_cost_tasks(job,self.manifest)
        self.assertEqual(len(planned),4)
        self.assertEqual(len({x['task_id'] for x in planned}),4)
        self.assertEqual(costs.read_timings(self.root/'absent',job,self.manifest,'test-helper',self.start),[])

    def test_single_tree_native_workers_and_tree_failures_are_explicit(self):
        job=dict(self.job,analysis='d2')
        manifest=dict(self.manifest,case_id='case',bsm_maps=5000,max_trees=200,bsm_maps_per_tree=500,
                      likelihood_contracts=['tip_survival','legacy_unconditioned'])
        mapping=self.row(analysis='d2',phase='run_bsm',case_id='case',planned_maps=5000,num_maps=5000,native_threads=2)
        trees=self.row(analysis='d2',phase='fit_over_trees',case_id='case',planned_trees=200,maps_per_tree=500,
            n_trees_planned=200,n_trees_completed=198,n_trees_incomplete=2,checkpoint_rows_after=200,
            started_utc='2020-01-02T00:00:02+00:00')
        points=costs.cost_points(self.read([mapping,trees],job,manifest))
        self.assertEqual([p['workload'] for p in points],['maps_5000','trees_200_maps_500'])
        self.assertEqual(points[1]['n_trees_incomplete'],2)
        with self.assertRaisesRegex(ValueError,'native thread'):self.read([dict(mapping,native_threads=1)],job,manifest)
        with self.assertRaisesRegex(ValueError,'incomplete'):self.read([dict(trees,checkpoint_rows_after=199)],job,manifest)

    def test_suite_counts_require_exact_ids_settings_workers_and_source(self):
        jobs=[dict(analysis='v1',scenario='x')]
        budget=dict(jobs=jobs,source_sha256={'scripts/run_v12_formal_suite.py':'source'})
        job=dict(id='one',analysis='v1',settings=jobs[0],reserved_cpus=2,argv=['--workers','2'])
        plan=dict(schema='rustbears.formal_suite.v1',budget_status='locked',budget_sha256='budget',
                  script_sha256='source',jobs=[job],cpu_budget=4)
        ledger=dict(plan_sha256='plan',source_unchanged=True,status='all_workflows_processed',outcomes=[dict(id='one',status='processed')])
        self.assertIn('one',costs.verify_plan(plan,ledger,budget,'budget','plan'))
        for update in [dict(outcomes=[]),dict(outcomes=[ledger['outcomes'][0]]*2),dict(source_unchanged=False),
                       dict(outcomes=[dict(id='one',status='failed_retained')]),dict(plan_sha256='changed')]:
            with self.subTest(update=update),self.assertRaises(ValueError):costs.verify_plan(plan,dict(ledger,**update),budget,'budget','plan')
        for update in [dict(reserved_cpus=8),dict(argv=['--workers','3']),dict(settings=dict(analysis='v1',scenario='changed'))]:
            with self.subTest(update=update),self.assertRaises(ValueError):costs.verify_plan(dict(plan,jobs=[dict(job,**update)]),ledger,budget,'budget','plan')

    def test_summary_only_processed_flag_does_not_replace_checkpoint(self):
        path=self.root/'summary.json';path.write_text('{"mode":"formal","n_records":2,"n_planned":2}')
        job=dict(self.job,output=str(self.root),result=str(path),settings=dict(analysis='v2c'))
        outcome=dict(result=dict(path=str(path),sha256=costs.sha(path),evidence_files={str(path):costs.sha(path)}))
        with self.assertRaises(FileNotFoundError) as error:costs.verify_processed_job(job,{},'budget',outcome)
        self.assertIn('preregistered_manifest',str(error.exception))

    def test_no_fcntl_import_on_portable_reader_path(self):
        # Exercise the actual import graph in a fresh process, including on Windows.
        result=subprocess.run([sys.executable,'-c','import sys; from scripts import v12_workflow_costs; assert "scripts.run_v12_formal_suite" not in sys.modules'],cwd=ROOT,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_dedicated_reader_rejects_smoke_even_with_a_valid_budget_envelope(self):
        from scripts import benchmark_v12_workflows as bench
        gate=self.root/'gates.json';gate.write_text('{"gate2_passed":true,"gate3_passed":true}')
        budget=dict(schema='rustbears.formal_budget.v1',status='locked',locked_utc='2020-01-01T00:00:00Z',
            gate_evidence=dict(path=gate.name,sha256=costs.sha(gate)),
            jobs=[dict(analysis='workflow_benchmark',**bench.DEFAULT)],
            source_sha256={name:costs.sha(ROOT/name) for name in ('scripts/benchmark_v12_workflows.py','scripts/v12_formal_protocol.py')})
        budget_path=self.root/'budget.json';budget_path.write_text(json.dumps(budget))
        folder=self.root/'smoke';folder.mkdir();(folder/'plan.json').write_text('{"mode":"smoke","formal":false}');(folder/'results.json').write_text('{"formal":false}')
        inventory=self.root/'inventory.json';inventory.write_text(json.dumps(dict(schema='rustbears.workflow_benchmark.inputs.v1',formal=True,
            budget=dict(path=budget_path.name,sha256=costs.sha(budget_path)),benchmark=folder.name)))
        with self.assertRaisesRegex(ValueError,'smoke or pilot'):costs.inspect_dedicated_benchmark(inventory)

    def test_all_failed_tree_maps_require_explicit_unavailable_marker(self):
        with self.assertRaises(FileNotFoundError):costs._verify_cost_maps(self.root,'joint_',19,0)
        (self.root/'no_complete_tree_maps.json').write_text('{"status":"no_complete_tree_maps","n_complete":0,"n_incomplete":20}')
        costs._verify_cost_maps(self.root,'joint_',19,0)
        (self.root/'joint_map_metadata.csv').write_text('map\n1\n')
        with self.assertRaisesRegex(ValueError,'zero complete trees'):costs._verify_cost_maps(self.root,'joint_',19,0)

    def test_plugin_cost_requires_all_B_rows_for_each_parameter(self):
        (self.root/'intervals.csv').write_text('param\nd\ne\nj\n')
        for name in ('d','e','j'):
            (self.root/f'curve_{name}.csv').write_text('value,profile_obj\n0,-1\n')
            (self.root/f'plugin_{name}.csv').write_text('rep_id,d_stat\n1,0\n2,0\n3,0\n')
        costs.verify_dedicated_phase(self.root,self.root,'profile_plugin',dict(n_boot_plugin=3),{})
        (self.root/'plugin_j.csv').write_text('rep_id,d_stat\n1,0\n2,0\n')
        with self.assertRaisesRegex(ValueError,'calibration replicate'):
            costs.verify_dedicated_phase(self.root,self.root,'profile_plugin',dict(n_boot_plugin=3),{})


if __name__=='__main__':unittest.main()
