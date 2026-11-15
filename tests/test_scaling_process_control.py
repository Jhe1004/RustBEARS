"""Engineering process-control regressions; no scientific model is evaluated."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock

from scripts import benchmark_v12_scaling as bench


@unittest.skipUnless(os.name == "posix" and hasattr(os, "wait4"), "POSIX process-group benchmark")
class ScalingProcessControlTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def wait_file(self, path, process, timeout=10):
        end = time.monotonic()+timeout
        while not path.exists() and time.monotonic() < end:
            if process.poll() is not None:
                self.fail("dummy controller exited before its child became ready")
            time.sleep(.01)
        self.assertTrue(path.exists(), "dummy child did not become ready")

    def test_normal_measurement_and_stop_before_launch(self):
        result = bench.measured([sys.executable, "-c", "print('dummy')"], self.folder, "normal", 10)
        self.assertEqual(result['returncode'], 0)
        self.assertFalse(result['timed_out'])
        self.assertFalse(result['interrupted'])
        control = bench.ExecutionControl(); control.request_stop(signal.SIGTERM)
        with mock.patch.object(bench.subprocess, 'Popen') as popen:
            with self.assertRaises(bench.BenchmarkInterrupted):
                bench.measured([sys.executable], self.folder, "never", 10, control)
            popen.assert_not_called()

    def signal_case(self, ignore=False):
        child = self.folder/'child.py'; ready = self.folder/'ready'; result = self.folder/'result.json'
        child.write_text("import os,signal,time\nfrom pathlib import Path\n"+
            ("signal.signal(signal.SIGTERM,signal.SIG_IGN)\n" if ignore else "")+
            f"Path({str(ready)!r}).write_text(str(os.getpid()))\ntime.sleep(60)\n")
        driver = self.folder/'driver.py'
        driver.write_text("import json,sys\nfrom pathlib import Path\n"+
            f"sys.path.insert(0,{str(bench.ROOT)!r})\n"+
            "from scripts import benchmark_v12_scaling as b\n"+
            "with b.ExecutionControl() as c:\n"+
            f" r=b.measured([sys.executable,{str(child)!r}],Path({str(self.folder)!r}),'dummy',30,c)\n"+
            f"Path({str(result)!r}).write_text(json.dumps(r))\n")
        other = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'], start_new_session=True)
        process = subprocess.Popen([sys.executable, str(driver)], start_new_session=True)
        try:
            self.wait_file(ready, process)
            child_pid = int(ready.read_text()); process.send_signal(signal.SIGTERM)
            self.assertEqual(process.wait(timeout=10), 0)
            value = json.loads(result.read_text())
            self.assertTrue(value['interrupted']); self.assertFalse(value['timed_out'])
            self.assertEqual(value['returncode'], -signal.SIGKILL if ignore else -signal.SIGTERM)
            with self.assertRaises(ProcessLookupError): os.kill(child_pid, 0)
            self.assertIsNone(other.poll(), "an unrelated process was affected")
        finally:
            for p in (process, other):
                if p.poll() is None:
                    os.killpg(p.pid, signal.SIGKILL); p.wait()

    def test_signal_stops_and_reaps_only_the_owned_process(self):
        self.signal_case()

    def test_signal_escalates_when_owned_process_ignores_term(self):
        self.signal_case(ignore=True)

    def test_reaped_leader_does_not_leave_ignoring_descendant_running(self):
        import fcntl
        ready = self.folder/'grandchild.ready'; leader_pid = self.folder/'leader.pid'
        lock = self.folder/'grandchild.lock'; result = self.folder/'result.json'
        child = self.folder/'grandchild.py'; leader = self.folder/'leader.py'; driver = self.folder/'driver.py'
        child.write_text("import fcntl,os,signal,time\nfrom pathlib import Path\n"+
            "signal.signal(signal.SIGTERM,signal.SIG_IGN)\n"+
            f"lock=open({str(lock)!r},'a+');fcntl.flock(lock,fcntl.LOCK_EX)\n"+
            f"Path({str(ready)!r}).write_text(str(os.getpid()))\ntime.sleep(60)\n")
        leader.write_text("import os,subprocess,sys,time\nfrom pathlib import Path\n"+
            f"Path({str(leader_pid)!r}).write_text(str(os.getpid()))\n"+
            f"subprocess.Popen([sys.executable,{str(child)!r}])\ntime.sleep(60)\n")
        driver.write_text("import json,sys\nfrom pathlib import Path\n"+
            f"sys.path.insert(0,{str(bench.ROOT)!r})\n"+
            "from scripts import benchmark_v12_scaling as b\nwith b.ExecutionControl() as c:\n"+
            f" r=b.measured([sys.executable,{str(leader)!r}],Path({str(self.folder)!r}),'family',30,c)\n"+
            f"Path({str(result)!r}).write_text(json.dumps(r))\n")
        process = subprocess.Popen([sys.executable,str(driver)],start_new_session=True)
        try:
            self.wait_file(ready,process); process.send_signal(signal.SIGTERM)
            self.assertEqual(process.wait(timeout=10),0)
            value = json.loads(result.read_text())
            self.assertTrue(value['interrupted']); self.assertEqual(value['returncode'],-signal.SIGTERM)
            with lock.open('a+') as check:
                deadline=time.monotonic()+5
                while True:
                    try:
                        fcntl.flock(check,fcntl.LOCK_EX|fcntl.LOCK_NB);break
                    except BlockingIOError:
                        if time.monotonic()>=deadline:raise
                        time.sleep(.01)
        finally:
            if process.poll() is None:process.kill();process.wait()
            if leader_pid.exists():
                try:os.killpg(int(leader_pid.read_text()),signal.SIGKILL)
                except ProcessLookupError:pass

    def test_signal_between_popen_and_assignment_is_not_lost(self):
        control = bench.ExecutionControl(); real = subprocess.Popen; spawned = []
        def launch(*args, **kwargs):
            process = real(*args, **kwargs); spawned.append(process.pid)
            control.request_stop(signal.SIGTERM)
            return process
        with mock.patch.object(bench.subprocess, 'Popen', side_effect=launch):
            result = bench.measured([sys.executable, '-c', 'import time; time.sleep(60)'],
                                    self.folder, 'window', 30, control)
        self.assertTrue(result['interrupted'])
        with self.assertRaises(ProcessLookupError): os.kill(spawned[0], 0)

    def test_timeout_still_reaps_owned_process(self):
        result = bench.measured([sys.executable, '-c', 'import time; time.sleep(60)'],
                                self.folder, 'timeout', .05)
        self.assertTrue(result['timed_out']); self.assertFalse(result['interrupted'])
        self.assertEqual(result['returncode'], -signal.SIGKILL)

    def test_inherited_lock_is_validated_preserved_and_forwarded(self):
        import fcntl
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(bench.inherited_pipeline_lock(), ())
        with (self.folder/'pipeline.lock').open('a+') as owner:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fd = owner.fileno()
            with mock.patch.dict(os.environ, {'RUSTBEARS_PIPELINE_LOCK_FD': str(fd)}):
                self.assertEqual(bench.inherited_pipeline_lock(), (fd,))
                code = f"import os; assert os.fstat({fd}).st_ino == {os.fstat(fd).st_ino}"
                result = bench.measured([sys.executable, '-c', code], self.folder, 'lock', 10,
                                        lock_fds=bench.inherited_pipeline_lock())
                self.assertEqual(result['returncode'], 0)
                os.fstat(fd)  # The benchmark must not close its caller's descriptor.
            closed_fd = os.dup(fd); os.close(closed_fd)
            with mock.patch.dict(os.environ, {'RUSTBEARS_PIPELINE_LOCK_FD': str(closed_fd)}):
                with self.assertRaises(OSError): bench.inherited_pipeline_lock()
        for value in ('0', '-1', 'abc'):
            with mock.patch.dict(os.environ, {'RUSTBEARS_PIPELINE_LOCK_FD': value}):
                with self.assertRaises(ValueError): bench.inherited_pipeline_lock()

    def test_child_keeps_pipeline_lock_after_controller_is_killed(self):
        import fcntl
        lock_path = self.folder/'pipeline.lock'; ready = self.folder/'ready'
        child = self.folder/'child.py'; driver = self.folder/'driver.py'
        child.write_text("import os,time\nfrom pathlib import Path\n"+
            f"Path({str(ready)!r}).write_text(str(os.getpid()))\ntime.sleep(60)\n")
        driver.write_text("import sys\nfrom pathlib import Path\n"+
            f"sys.path.insert(0,{str(bench.ROOT)!r})\n"+
            "from scripts import benchmark_v12_scaling as b\n"+
            f"b.measured([sys.executable,{str(child)!r}],Path({str(self.folder)!r}),'locked',30,lock_fds=b.inherited_pipeline_lock())\n")
        child_pid = None; controller = None
        with lock_path.open('a+') as owner:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fd = owner.fileno()
            env = dict(os.environ, RUSTBEARS_PIPELINE_LOCK_FD=str(fd))
            controller = subprocess.Popen([sys.executable, str(driver)], env=env,
                                          pass_fds=(fd,), start_new_session=True)
            try:
                self.wait_file(ready, controller); child_pid = int(ready.read_text())
            finally:
                controller.kill(); controller.wait()
        try:
            with lock_path.open('a+') as challenger:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(challenger, fcntl.LOCK_EX | fcntl.LOCK_NB)
        finally:
            if child_pid is not None:
                try: os.killpg(child_pid, signal.SIGKILL)
                except ProcessLookupError: pass
        # Release is observable without treating an orphan's temporary zombie PID
        # as a still-running workload.
        with lock_path.open('a+') as challenger:
            deadline = time.monotonic()+5
            while True:
                try:
                    fcntl.flock(challenger, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    if time.monotonic() >= deadline: raise
                    time.sleep(.01)

    def test_interruption_retains_partial_measurement_without_launching_next(self):
        control = bench.ExecutionControl(); calls = []
        def inputs(folder, *_):
            (folder/'tree.newick').write_text('dummy')
            (folder/'geography.data').write_text('dummy')
        def measured(*args):
            calls.append(args[2]); control.request_stop(signal.SIGTERM)
            return dict(returncode=-signal.SIGTERM, interrupted=True, timed_out=False)
        args = SimpleNamespace(repeats=1, oracle_library=self.folder, timeout=10, mode='smoke')
        manifest = dict(tips=[3], states=[4], params={})
        with mock.patch.object(bench, 'make_inputs', side_effect=inputs), \
             mock.patch.object(bench, 'commands_for', return_value={'a':['dummy'], 'b':['dummy']}), \
             mock.patch.object(bench, 'measured', side_effect=measured):
            with self.assertRaises(bench.BenchmarkInterrupted):
                bench.execute_workloads(args,self.folder,manifest,{},self.folder/'unused',{},control,())
        self.assertEqual(len(calls), 1)
        for filename in ('progress.json', 'results.json'):
            value = json.loads((self.folder/filename).read_text())
            self.assertEqual(value['status'], 'interrupted')
            self.assertFalse(value['formal'])
            self.assertEqual(len(value['current_attempt']['execution']), 1)
            self.assertEqual(value['records'], [])

    def test_spawn_failure_is_retained_as_failure(self):
        args = SimpleNamespace(repeats=1, oracle_library=self.folder, timeout=10, mode='smoke')
        def inputs(folder, *_):
            (folder/'tree.newick').write_text('dummy'); (folder/'geography.data').write_text('dummy')
        with mock.patch.object(bench, 'make_inputs', side_effect=inputs), \
             mock.patch.object(bench, 'commands_for', return_value={'missing':[str(self.folder/'no-executable')]}):
            with self.assertRaises(FileNotFoundError):
                bench.execute_workloads(args,self.folder,dict(tips=[3],states=[4],params={}),{},
                    self.folder/'unused',{},bench.ExecutionControl(),())
        value = json.loads((self.folder/'results.json').read_text())
        self.assertEqual(value['status'], 'failed_retained')
        self.assertIn('FileNotFoundError', value['error'])


if __name__ == '__main__':
    unittest.main()
