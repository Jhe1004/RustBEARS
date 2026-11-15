"""Plan/cache/process-control checks; no formal cost experiment is executed."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts import benchmark_v12_workflows as bench


class WorkflowBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.path=Path(self.temp.name)

    def test_default_grid_is_complete_and_plan_only(self):
        result=subprocess.run([sys.executable,str(ROOT/'scripts/benchmark_v12_workflows.py'),'--output',str(self.path)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        plan=json.loads((self.path/'plan.json').read_text())
        self.assertEqual(len(plan['cells']),18);self.assertEqual(plan['settings']['states'],[16,64,256])
        self.assertEqual(plan['settings']['phase_timeout_seconds'],14400)
        self.assertEqual(plan['phases'],list(bench.PHASES));self.assertIsNone(plan['budget'])
        self.assertFalse((self.path/'results.json').exists())
        self.assertFalse(any(p.name=='worker_result.json' for p in self.path.rglob('*')))

    def test_formal_execute_requires_locked_budget_before_output_creation(self):
        output=self.path/'forbidden'
        result=subprocess.run([sys.executable,str(ROOT/'scripts/benchmark_v12_workflows.py'),'--output',str(output),'--execute'],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0);self.assertIn('own locked budget',result.stderr);self.assertFalse(output.exists())

    def test_paired_workers_share_input_seed_without_replacing_replicates(self):
        a=dict(states=16,replicate=1,workers=1);b=dict(a,workers=8)
        self.assertEqual(bench.seeds_for(bench.DEFAULT,a),bench.seeds_for(bench.DEFAULT,b))
        self.assertNotEqual(bench.seeds_for(bench.DEFAULT,a),bench.seeds_for(bench.DEFAULT,dict(a,replicate=2)))
        self.assertEqual(len(set(bench.seeds_for(bench.DEFAULT,a).values())),12)

    def test_invalid_or_duplicate_dimensions_rejected(self):
        for change in (dict(states=[16,16]),dict(states=[57]),dict(workers=[0]),dict(workers=[48]),dict(repeats=0),dict(n_tips=4),dict(phase_timeout_seconds=0),dict(maps_per_tree=10)):
            with self.subTest(change=change),self.assertRaises(ValueError):bench.validate_settings(dict(bench.DEFAULT,**change))

    def test_formal_core_must_belong_to_complete_execution_checkout(self):
        with patch.object(bench,'ROOT',self.path):
            with self.assertRaisesRegex(ValueError,'full source checkout'):bench.require_checkout_core()

    def test_all_initial_and_tree_templates_share_physical_parameter_domain(self):
        for name in bench.MODELS:
            definition=bench.template(name)
            self.assertEqual((definition.params['d'].min,definition.params['d'].max,definition.params['d'].init),(0,2,.015))
            self.assertEqual((definition.params['e'].min,definition.params['e'].max,definition.params['e'].init),(0,2,.008))
            if '+J' in name:self.assertEqual((definition.params['j'].min,definition.params['j'].init),(0,.03))

    def test_scaled_trees_are_deterministic_computational_inputs(self):
        import rustbears as b
        cell=bench.cells_for(bench.SMOKE)[0]
        a=bench.prepare_inputs(self.path/'a',bench.SMOKE,cell);z=bench.prepare_inputs(self.path/'z',bench.SMOKE,dict(cell,workers=8))
        self.assertEqual(a,z);self.assertEqual(a['tree_scales'],[.9,1.1])
        self.assertIn('not posterior',a['interpretation'])
        trees=[b.read_newick(line) for line in (self.path/'a/scaled_trees.newick').read_text().splitlines()]
        self.assertEqual([t.tip_count for t in trees],[8,8])

    def test_fit_cache_tampering_detected_before_object_load(self):
        (self.path/'fits.json').write_text('{"array_sha256":"wrong","records":[]}')
        (self.path/'fits.npz').write_bytes(b'arbitrary test payload')
        with self.assertRaisesRegex(ValueError,'metadata changed'):bench.load_fits(self.path,self.path,'wrong')
        with self.assertRaisesRegex(ValueError,'arrays changed'):bench.load_fits(self.path,self.path,bench.sha(self.path/'fits.json'))

    def test_after_suite_requires_exact_terminal_ids_and_no_reserved_cpu(self):
        bench.write(self.path/'suite_plan.json',dict(jobs=[dict(id='a'),dict(id='b')]))
        baseline=dict(plan_sha256=bench.sha(self.path/'suite_plan.json'),status='completed_with_retained_failures',source_unchanged=True,
            outcomes=[dict(id='a',status='processed'),dict(id='b',status='failed_retained')],cpu_reservations=dict(active_reserved_cpus=0))
        bench.write(self.path/'suite_ledger.json',baseline)
        self.assertEqual(bench.completed_statistical_suite(self.path)['status'],'completed_with_retained_failures')
        for change in (dict(outcomes=baseline['outcomes'][:1]),dict(status='running'),dict(cpu_reservations=dict(active_reserved_cpus=8)),dict(plan_sha256='bad')):
            bench.write(self.path/'suite_ledger.json',dict(baseline,**change))
            with self.subTest(change=change),self.assertRaises(ValueError):bench.completed_statistical_suite(self.path)

    def test_predecessor_binding_rejects_unrelated_terminal_suite(self):
        bench.write(self.path/'suite_plan.json',dict(jobs=[dict(id='a',settings={})],budget_sha256='wrong'))
        bench.write(self.path/'suite_ledger.json',dict(outcomes=[]))
        expected=dict(jobs=[dict(analysis='v1')],source_sha256={})
        with self.assertRaisesRegex(ValueError,'complete locked statistical budget'):bench.completed_statistical_suite(self.path,expected,'authorized')
        original=self.path/'statistical_budget.json'
        bench.write(original,dict(schema='rustbears.formal_budget.v1',status='locked',gate_evidence=dict(sha256='gates')))
        budget=dict(predecessor_budget=dict(path=original.name,sha256=bench.sha(original)),
            predecessor_statistical_budget_sha256=bench.sha(original),gate_evidence=dict(sha256='gates'))
        self.assertEqual(bench.predecessor_budget(budget,self.path/'cost_budget.json')[1],bench.sha(original))
        budget['predecessor_statistical_budget_sha256']='another'
        with self.assertRaisesRegex(ValueError,'original statistical budget'):bench.predecessor_budget(budget,self.path/'cost_budget.json')

    @unittest.skipUnless(hasattr(os,'wait4'),'bounded phase execution requires POSIX wait4/process groups')
    def test_timeout_preserves_logs_and_censored_elapsed(self):
        script=self.path/'sleep.py'
        script.write_text('import json,time,sys\nfrom pathlib import Path\np=Path(sys.argv[1]);(p/"api_started.json").write_text(json.dumps({"monotonic":time.monotonic()}));print("retained stdout",flush=True);print("retained stderr",file=sys.stderr,flush=True);time.sleep(5)\n')
        result=bench.measured([sys.executable,str(script),str(self.path)],self.path,.2,threading.Event())
        self.assertTrue(result['timed_out']);self.assertFalse(result['interrupted']);self.assertGreater(result['process_wall_seconds'],.19)
        self.assertGreater(result['api_right_censored_elapsed_seconds'],0);self.assertLess(result['process_wall_seconds'],3)
        self.assertIn('retained stdout',(self.path/'stdout.log').read_text());self.assertIn('retained stderr',(self.path/'stderr.log').read_text())

    @unittest.skipUnless(hasattr(os,'wait4'),'bounded phase execution requires POSIX wait4/process groups')
    def test_stop_interrupts_owned_process_without_waiting_for_timeout(self):
        stop=threading.Event();thread=threading.Timer(.15,stop.set);thread.start()
        try:result=bench.measured([sys.executable,'-c','import time; time.sleep(5)'],self.path,30,stop)
        finally:thread.cancel()
        self.assertTrue(result['interrupted']);self.assertFalse(result['timed_out']);self.assertLess(result['process_wall_seconds'],3)

    @unittest.skipUnless(hasattr(os,'wait4'),'bounded phase execution requires POSIX wait4/process groups')
    def test_export_timeout_does_not_censor_already_completed_api(self):
        script=self.path/'slow_export.py'
        script.write_text('import json,time,sys\nfrom pathlib import Path\np=Path(sys.argv[1]);(p/"api_started.json").write_text(json.dumps({"monotonic":time.monotonic()}));(p/"workflow_timings.jsonl").write_text(json.dumps({"status":"returned","wall_seconds":.01})+"\\n");time.sleep(5)\n')
        result=bench.measured([sys.executable,str(script),str(self.path)],self.path,.2,threading.Event())
        self.assertTrue(result['timed_out']);self.assertIsNone(result['api_right_censored_elapsed_seconds'])
        self.assertEqual(result['terminal_api_status'],'returned');self.assertEqual(result['terminal_api_wall_seconds'],.01)

    @unittest.skipUnless(hasattr(os,'wait4'),'phase launcher execution requires POSIX')
    def test_terminal_timeout_is_not_rerun_and_its_evidence_is_immutable(self):
        # A private driver replaces only process execution; it never calls an
        # inference API or uses formal mode. Production timeout policy is intact.
        driver=self.path/'driver.py'
        driver.write_text(f'''import json,sys
from pathlib import Path
sys.path.insert(0,{str(ROOT)!r})
from scripts import benchmark_v12_workflows as b
b.os.getpriority=lambda *_:0
counter=Path({str(self.path/'calls.json')!r})
def dummy(command,folder,timeout,stop,lock_fds=()):
 n=json.loads(counter.read_text()) if counter.exists() else 0
 counter.write_text(json.dumps(n+1))
 (folder/'stdout.log').write_text('retained dummy timeout')
 (folder/'stderr.log').write_text('retained dummy stderr')
 return dict(returncode=-9,timed_out=True,interrupted=False,process_wall_seconds=.2,api_right_censored_elapsed_seconds=.1)
b.measured=dummy
b.main()
''')
        output=self.path/'dummy_costs';command=[sys.executable,str(driver),'--mode','smoke','--output',str(output),'--execute']
        first=subprocess.run(command,capture_output=True,text=True);self.assertEqual(first.returncode,0,first.stderr)
        records=json.loads((output/'results.json').read_text())['records']
        self.assertEqual(sum(r['status']=='timeout_retained' for r in records),2)
        self.assertEqual(sum(r['status']=='blocked_initial_fit_failure' for r in records),4)
        second=subprocess.run(command+['--resume'],capture_output=True,text=True);self.assertEqual(second.returncode,0,second.stderr)
        self.assertEqual(json.loads((self.path/'calls.json').read_text()),2)
        # Reproduce the last durable coordinator record after parent SIGKILL:
        # the child lock has now closed, but no exit/RSS/clock was observed.
        ledger=json.loads((output/'results.json').read_text());orphan=ledger['records'][0]
        orphan['status']='running';orphan['attempts'][-1]=dict(directory=orphan['attempts'][-1]['directory'],status='running')
        bench.write(output/'results.json',ledger)
        resumed=subprocess.run(command+['--resume'],capture_output=True,text=True);self.assertEqual(resumed.returncode,0,resumed.stderr)
        repaired=json.loads((output/'results.json').read_text())['records'][0]['attempts']
        self.assertEqual(len(repaired),2);self.assertEqual(repaired[0]['status'],'unobserved_process_exit')
        for key in ('returncode','process_wall_seconds','process_cpu_seconds','process_peak_rss_bytes'):
            self.assertIsNone(repaired[0][key])
        self.assertEqual(json.loads((self.path/'calls.json').read_text()),3)
        path=Path(records[0]['attempts'][0]['directory'])/'stdout.log';path.write_text('changed')
        changed=subprocess.run(command+['--resume'],capture_output=True,text=True)
        self.assertNotEqual(changed.returncode,0);self.assertIn('evidence changed',changed.stderr)


if __name__=='__main__':unittest.main()
