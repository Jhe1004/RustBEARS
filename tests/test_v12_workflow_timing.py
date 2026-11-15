"""Timing evidence must retain failed calls and identify partial checkpoint resumes."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from scripts.v12_workflow_timing import WorkflowTimer,checkpoint_rows


class WorkflowTimingTests(unittest.TestCase):
    def test_actual_call_result_and_checkpoint_progress_are_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'timings.jsonl';checkpoint=Path(d)/'checkpoint.jsonl'
            checkpoint.write_text('{}\n{}\n')
            result=SimpleNamespace(n_planned=4,n_failed=1)
            def work():
                with checkpoint.open('a') as f:f.write('{}\n')
                return result
            actual=WorkflowTimer(path,resume_requested=True).run('bootstrap_lrt',work,checkpoint=checkpoint)
            self.assertIs(actual,result)
            row=json.loads(path.read_text())
            self.assertEqual((row['checkpoint_rows_before'],row['checkpoint_rows_after']),(1,2))
            self.assertEqual(row['n_failed'],1)
            self.assertEqual(row['status'],'returned')
            self.assertGreater(row['wall_seconds'],0)
            self.assertIn('excludes outer result-table export',row['timing_scope'])

    def test_exception_is_rethrown_and_retained_without_success_record(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'timings.jsonl'
            def fail():raise ValueError('failed fit')
            with self.assertRaisesRegex(ValueError,'failed fit'):
                WorkflowTimer(path).run('profile_bootstrap_plugin',fail)
            row=json.loads(path.read_text())
            self.assertEqual(row['status'],'raised')
            self.assertEqual(row['error'],'ValueError: failed fit')
            self.assertGreater(row['wall_seconds'],0)

    def test_first_launch_resume_flag_is_not_cached_work_and_partial_tail_is_not_complete(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'timings.jsonl';checkpoint=Path(d)/'checkpoint.jsonl'
            checkpoint.write_bytes(b'{}\n{"partial":')
            self.assertEqual(checkpoint_rows(checkpoint),0)
            WorkflowTimer(path,resume_requested=True).run('profile_asymptotic_chi2',lambda:None,checkpoint=checkpoint)
            row=json.loads(path.read_text())
            self.assertTrue(row['resume_requested'])
            self.assertEqual(row['checkpoint_rows_before'],0)


if __name__=='__main__':unittest.main()
