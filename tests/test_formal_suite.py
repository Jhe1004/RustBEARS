"""Test the scheduler with short local dummy processes; never run science jobs.

The driver replaces command construction and result checking only inside its
private test process. Temporary fixture budgets do not authorize formal studies.
"""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

try:
    import fcntl
except ImportError:  # The server launcher owns POSIX process groups and flock.
    fcntl=None

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
if fcntl is not None:
    import run_v12_formal_suite as suite
else:
    suite=None


@unittest.skipIf(fcntl is None,'formal-suite process-group launcher requires POSIX')
class FormalSuiteTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)

    def test_summary_counters_alone_cannot_certify_any_workflow(self):
        for analysis in ('v1','v2b','v2c','v3','v4','d1','d2'):
            with self.subTest(analysis=analysis):
                directory=self.path/analysis;directory.mkdir()
                settings=dict(analysis=analysis,scenario='unit',n_tips=4,areas=2,seed=1,
                    finite_b_extrapolation=False,outer_replicates=2,replicates=2,settings={},case_id='unit')
                data=dict(mode='formal',formal=True,locked_budget_sha256='unit',budget_sha256='unit',
                    n_records=2,n_planned=2,outcomes=[{},{}],settings={},source_unchanged=True,
                    status='all_planned_datasets_processed' if analysis=='v3' else 'completed' if analysis=='v4' else 'complete')
                result=directory/('summary.json' if analysis in ('v1','v2b','v2c') else 'results.json')
                budget=dict(status='locked',jobs=[settings],source_sha256={'scripts/run_uncertainty_case_study.py':'unit-source'})
                budget_path=directory/'budget.json';budget_path.write_text(json.dumps(budget));digest=suite.sha(budget_path)
                data.update(locked_budget_sha256=digest,budget_sha256=digest)
                result.write_text(json.dumps(data))
                if analysis in ('d1','d2'):
                    # Supply the formal provenance envelope so rejection reaches
                    # the missing scientific artifacts, not absent call fields.
                    from v12_case_evidence import MODELS,CONTRACTS
                    manifest=dict(settings,mode='formal',formal=True,original_posterior_required=analysis=='d2',
                        models=MODELS,likelihood_contracts=CONTRACTS,workflow_script_sha256='unit-source',
                        locked_budget=dict(sha256=digest,payload=budget,locked_utc='2020-01-01T00:00:00+00:00'))
                    (directory/'manifest.json').write_text(json.dumps(manifest))
                    (directory/'execution_started.json').write_text(json.dumps(dict(budget_sha256=digest,
                        gate2_gate3_verified=True,first_started_utc='2021-01-01T00:00:00+00:00')))
                    pattern='planned contract/model family'
                else:pattern='preregistered_manifest|execution_manifest|manifest.json'
                with self.assertRaisesRegex((ValueError,OSError),pattern):
                    suite.check_result(dict(analysis=analysis,result=str(result),output=str(directory),settings=settings,
                        argv=['unit-command','--budget-json',str(budget_path)]),digest)

    def make_driver(self,n=1,seconds=.1,*,jobs=None,cpu_budget=2,workers=1,large_workers=None,cpu_count_override=None,priority_override=0):
        script=self.path/'dummy_worker.py'
        script.write_text('''import fcntl,json,os,sys,time
from pathlib import Path
base=Path(sys.argv[1]); state=base/'dummy_state.json'
def change(delta):
 with (base/'dummy_state.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  v=json.loads(state.read_text()) if state.exists() else dict(active=0,peak=0,starts=0,pids=[],thread_env=[],active_reserved=0,peak_reserved=0,allocations=[])
  v['active']+=delta;v['peak']=max(v['peak'],v['active'])
  v['active_reserved']+=delta*int(sys.argv[3]);v['peak_reserved']=max(v['peak_reserved'],v['active_reserved'])
  if delta>0:
   v['allocations'].append(int(sys.argv[3]))
   v['starts']+=1;v['pids'].append(os.getpid());v['thread_env'].append({k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
  state.write_text(json.dumps(v))
change(1)
time.sleep(float(sys.argv[2]))
change(-1)
''')
        driver=self.path/'driver.py'
        driver.write_text(f'''import sys
from pathlib import Path
sys.path.insert(0,{str(ROOT/'scripts')!r})
import run_v12_formal_suite as suite
# Host CI runners may inherit nice != 0. Isolate this environmental input for
# sleeping dummy workers only; the production launcher's gate stays unchanged.
suite.os.getpriority=lambda which,who: {priority_override!r}
if {cpu_count_override!r} is not None:
 suite.os.cpu_count=lambda: {cpu_count_override!r}
def dummy(job,index,output,budget_path,revbayes,d2,workers):
 name=f'{{index:02d}}_unit';directory=output/name
 return dict(id=name,analysis=job['analysis'],settings=job,argv=[sys.executable,{str(script)!r},{str(self.path)!r},{str(seconds)!r},str(workers)],output=str(directory),result=str(directory/'unused.json'),reserved_cpus=workers,cwd={str(ROOT)!r},environment=suite.THREAD_ENV)
suite.command_for=dummy
suite.check_result=lambda job,digest: dict(path=job['result'],sha256='unit-test-only',meaning='dummy process completed')
suite.main()
''')
        gate=self.path/'unit_gates.json';gate.write_text(json.dumps(dict(gate2_passed=True,gate3_passed=True)))
        budget=dict(schema='rustbears.formal_budget.v1',status='locked',locked_utc='2020-01-01T00:00:00+00:00',
            gate_evidence=dict(path=str(gate),sha256=suite.sha(gate)),
            source_sha256={name:suite.sha(ROOT/name) for name in suite.VALIDATOR_SOURCES},
            jobs=jobs if jobs is not None else [dict(analysis='v1',unit_test=i) for i in range(n)])
        budget_path=self.path/'unit_budget.json';budget_path.write_text(json.dumps(budget))
        output=self.path/'suite'
        self.argv=[sys.executable,str(driver),'--budget-json',str(budget_path),'--output',str(output),
                   '--revbayes',sys.executable,'--d2-input-dir',str(self.path),'--cpu-budget',str(cpu_budget),'--workers-per-job',str(workers)]
        if large_workers is not None:self.argv+=['--large-v1-workers',str(large_workers)]
        self.output=output
        return output

    def run_driver(self,*flags):
        return subprocess.run([*self.argv,*flags],capture_output=True,text=True,timeout=15,
                              env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))

    def test_default_plan_runs_no_worker(self):
        self.make_driver(3)
        result=self.run_driver();self.assertEqual(result.returncode,0,result.stderr)
        plan=json.loads((self.output/'suite_plan.json').read_text())
        self.assertEqual(plan['max_concurrent_jobs'],2)
        self.assertFalse((self.path/'dummy_state.json').exists())
        self.assertFalse((self.output/'suite_ledger.json').exists())
        self.assertLessEqual({p.name for p in self.output.iterdir()},{'suite_plan.json','suite.lock'})

    def test_execution_respects_cpu_and_native_thread_caps(self):
        self.make_driver(4,.2)
        result=self.run_driver('--execute');self.assertEqual(result.returncode,0,result.stderr)
        state=json.loads((self.path/'dummy_state.json').read_text())
        self.assertEqual(state['starts'],4);self.assertLessEqual(state['peak'],2)
        self.assertEqual(state['active'],0)
        self.assertTrue(all(set(row.values())=={'1'} for row in state['thread_env']))
        ledger=json.loads((self.output/'suite_ledger.json').read_text())
        self.assertEqual(len(ledger['outcomes']),4)

    def test_nonzero_priority_is_rejected_before_any_execution_or_plan(self):
        self.make_driver(priority_override=5)
        result=self.run_driver('--execute')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('normal process priority (nice=0)',result.stderr)
        self.assertFalse((self.path/'dummy_state.json').exists())
        self.assertFalse(self.output.exists())

    def test_fixed_priority_and_large_allocation_leave_scientific_settings_unchanged(self):
        jobs=[dict(analysis='v1',n_tips=40,areas=4,scenario='interior',seed=1),
              dict(analysis='v2b',scenario='null_e0',seed=2),
              dict(analysis='v1',n_tips=150,areas=5,scenario='j_zero',seed=3),
              dict(analysis='d1',seed=4),dict(analysis='v3',seed=5),dict(analysis='d2',seed=6),
              dict(analysis='v4',seed=7),dict(analysis='v1',n_tips=150,areas=5,scenario='interior',seed=8)]
        self.make_driver(jobs=jobs,cpu_budget=40,workers=2,large_workers=8,cpu_count_override=64)
        result=self.run_driver();self.assertEqual(result.returncode,0,result.stderr)
        plan=json.loads((self.output/'suite_plan.json').read_text())
        self.assertEqual([r['settings']['seed'] for r in plan['jobs']],[4,6,5,7,3,8,2,1])
        self.assertEqual({r['settings']['seed']:r['settings'] for r in plan['jobs']},{r['seed']:r for r in jobs})
        for job in plan['jobs']:
            expected=8 if suite.is_large_v1(job['settings']) else 2
            self.assertEqual(job['reserved_cpus'],expected)
            self.assertEqual(int(job['argv'][-1]),expected)
        # Check the real workflow command, not only the dummy argument builder.
        task=dict(analysis='v1',n_tips=150,areas=5,scenario='interior',seed=20261115)
        command=suite.command_for(task,1,self.path,self.path/'budget.json',Path(sys.executable),self.path,
                                  suite.allocated_workers(task,2,8))
        self.assertEqual(command['argv'][command['argv'].index('--workers')+1],'8')
        self.assertEqual(command['reserved_cpus'],8)

    def test_mixed_eight_and_two_cpu_jobs_respect_reserved_budget(self):
        # Only short sleeping processes run. The mocked inventory permits this
        # accounting test on two-core CI hosts without consuming 40 CPU cores.
        jobs=[dict(analysis=a,unit_test=a) for a in ('d1','d2','v3','v4')]
        jobs += [dict(analysis='v1',n_tips=150,areas=5,unit_test=i) for i in range(4)]
        jobs += [dict(analysis='v2c',unit_test=i) for i in range(2)]
        self.make_driver(jobs=jobs,seconds=.3,cpu_budget=40,workers=2,large_workers=8,cpu_count_override=64)
        result=self.run_driver('--execute');self.assertEqual(result.returncode,0,result.stderr)
        state=json.loads((self.path/'dummy_state.json').read_text())
        self.assertEqual(state['starts'],10);self.assertEqual(state['allocations'].count(8),4)
        self.assertEqual(state['allocations'].count(2),6)
        self.assertLessEqual(state['peak_reserved'],40);self.assertEqual(state['active_reserved'],0)
        ledger=json.loads((self.output/'suite_ledger.json').read_text())
        self.assertLessEqual(ledger['cpu_reservations']['peak_reserved_cpus_this_launch'],40)
        self.assertEqual(ledger['cpu_reservations']['active_reserved_cpus'],0)
        self.assertTrue(all(r['total_reserved_cpus_at_admission']<=40 for r in ledger['outcomes']))

    def test_allocation_cannot_change_on_resume_and_default_matches_regular(self):
        jobs=[dict(analysis='v1',n_tips=150,areas=5,unit_test=1)]
        self.make_driver(jobs=jobs,cpu_budget=8,workers=2,cpu_count_override=64)
        result=self.run_driver();self.assertEqual(result.returncode,0,result.stderr)
        path=self.output/'suite_plan.json';before=path.read_bytes()
        self.assertEqual(json.loads(before)['large_v1_workers'],2)
        same=self.run_driver('--large-v1-workers','2');self.assertEqual(same.returncode,0,same.stderr)
        changed=self.run_driver('--large-v1-workers','8','--resume')
        self.assertNotEqual(changed.returncode,0);self.assertIn('different immutable plan',changed.stderr)
        self.assertEqual(path.read_bytes(),before)

    def test_signal_releases_waiting_weighted_jobs_without_starting_them(self):
        jobs=[dict(analysis='v1',n_tips=150,areas=5,unit_test=1),dict(analysis='v1',n_tips=40,areas=4,unit_test=2)]
        self.make_driver(jobs=jobs,seconds=5,cpu_budget=8,workers=2,large_workers=8,cpu_count_override=64)
        process=subprocess.Popen([*self.argv,'--execute'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,
                                 start_new_session=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        try:
            deadline=time.monotonic()+5
            while not (self.path/'dummy_state.json').exists() and time.monotonic()<deadline:time.sleep(.02)
            self.assertTrue((self.path/'dummy_state.json').exists())
            os.kill(process.pid,signal.SIGTERM);process.wait(timeout=4)
            state=json.loads((self.path/'dummy_state.json').read_text());self.assertEqual(state['starts'],1)
            ledger=json.loads((self.output/'suite_ledger.json').read_text())
            self.assertEqual(ledger['status'],'interrupted')
            self.assertEqual({r['id']:r['status'] for r in ledger['outcomes']},{'01_unit':'interrupted','02_unit':'not_started'})
            self.assertEqual(ledger['cpu_reservations']['active_reserved_cpus'],0)
            self.assertEqual(ledger['cpu_reservations']['waiting_jobs'],[])
            self.assertFalse((self.output/'02_unit').exists())
        finally:
            if process.poll() is None:process.kill();process.wait(timeout=3)

    def seed_ledger(self,rows):
        result=self.run_driver();self.assertEqual(result.returncode,0,result.stderr)
        path=self.output/'suite_ledger.json'
        path.write_text(json.dumps(dict(plan_sha256=suite.sha(self.output/'suite_plan.json'),outcomes=rows)))
        return path

    def test_terminal_failure_never_automatically_reruns(self):
        self.make_driver()
        failed=dict(id='01_unit',status='failed_retained',error='retained original numerical failure')
        self.seed_ledger([failed]);result=self.run_driver('--execute','--resume')
        self.assertNotEqual(result.returncode,0)
        self.assertFalse((self.path/'dummy_state.json').exists())
        ledger=json.loads((self.output/'suite_ledger.json').read_text())
        self.assertEqual(ledger['outcomes'],[failed])

    def test_interrupted_owned_task_resumes_without_erasing_its_directory(self):
        self.make_driver()
        self.seed_ledger([dict(id='01_unit',status='interrupted')])
        directory=self.output/'01_unit';directory.mkdir();marker=directory/'prior_checkpoint.txt';marker.write_text('retained')
        result=self.run_driver('--execute','--resume');self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(marker.read_text(),'retained')
        self.assertEqual(json.loads((self.path/'dummy_state.json').read_text())['starts'],1)
        self.assertEqual(json.loads((self.output/'suite_ledger.json').read_text())['outcomes'][0]['status'],'processed')

    def test_existing_ledger_does_not_authorize_unowned_job_directory(self):
        self.make_driver();path=self.seed_ledger([]);before=path.read_bytes()
        (self.output/'01_unit').mkdir()
        result=self.run_driver('--execute','--resume')
        self.assertNotEqual(result.returncode,0);self.assertIn('unowned',result.stderr)
        self.assertEqual(path.read_bytes(),before)
        self.assertFalse((self.path/'dummy_state.json').exists())

    def test_duplicate_or_unknown_ledger_ids_are_rejected_without_rewriting(self):
        self.make_driver()
        for rows in ([dict(id='01_unit',status='failed_retained')]*2,[dict(id='not_in_plan',status='failed_retained')],
                     [dict(id='01_unit',status='invented')]):
            with self.subTest(rows=rows):
                path=self.seed_ledger(rows);before=path.read_bytes()
                result=self.run_driver('--execute','--resume')
                self.assertNotEqual(result.returncode,0)
                self.assertEqual(path.read_bytes(),before)
                self.assertFalse((self.path/'dummy_state.json').exists())

    def test_child_retains_lock_after_abrupt_launcher_death(self):
        self.make_driver(1,2.)
        process=subprocess.Popen([*self.argv,'--execute'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,
                                 start_new_session=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        try:
            deadline=time.monotonic()+5
            while not (self.path/'dummy_state.json').exists() and time.monotonic()<deadline:time.sleep(.02)
            self.assertTrue((self.path/'dummy_state.json').exists(),'dummy child failed to start')
            os.kill(process.pid,signal.SIGKILL);process.wait(timeout=3)
            attempt=self.run_driver('--execute','--resume')
            self.assertNotEqual(attempt.returncode,0)
            self.assertIn('owns this output',attempt.stderr)
            self.assertEqual(json.loads((self.path/'dummy_state.json').read_text())['starts'],1)
            with (self.output/'suite.lock').open('a+') as lock:
                deadline=time.monotonic()+5
                while True:
                    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);break
                    except BlockingIOError:
                        if time.monotonic()>=deadline:self.fail('child did not release inherited lock')
                        time.sleep(.02)
            self.assertEqual(json.loads((self.path/'dummy_state.json').read_text())['active'],0)
        finally:
            if process.poll() is None:process.kill();process.wait(timeout=3)


if __name__=='__main__':unittest.main()
