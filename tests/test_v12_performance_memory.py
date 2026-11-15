"""Synthetic arithmetic/reader fixtures only; never render quantitative figures."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.plot_v12_performance import benchmark_rows, benchmark_counts, CONTRACT


class PerformanceMemoryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)/'arithmetic_fixture.json'
        execution=dict(returncode=0,timed_out=False,wall_seconds=8.,peak_rss_bytes=40*1024**2)
        record=dict(preflight=False,tips=3,states=4,replicate=0,input_sha256={'tree.newick':'fixture'},
            execution={},comparisons={})
        for name,reference in [('original_biogeobears','rust_original'),('revbayes','rust_revbayes'),('lagrange_ng','rust_lagrange')]:
            record['execution'][name]=dict(execution)
            record['execution'][reference]=dict(execution,wall_seconds=2.,peak_rss_bytes=10*1024**2)
            record['comparisons'][name]=dict(rust_reference=reference,alignment_passed=True,maximum_lnL_error=0.)
        self.data=dict(formal=True,source_stable=True,original_oracle_stable=True,
            manifest=dict(mode='formal',tips=[3],states=[4],repeats=2),
            records=[record,copy.deepcopy(record)])
        self.data['records'][1]['replicate']=1

    def read(self):
        self.path.write_text(json.dumps(self.data))
        return benchmark_rows([('arithmetic-host',self.path)])[0]

    def test_all_matched_repeats_keep_both_ratios_and_raw_bytes(self):
        before=copy.deepcopy(self.data)
        rows=self.read()
        self.assertEqual(len(rows),6)
        self.assertEqual({r['replicate'] for r in rows},{0,1})
        self.assertTrue(all(r['time_ratio']==4 and r['memory_ratio']==4 for r in rows))
        self.assertTrue(all(r['external_peak_rss_mib']==40 and r['rust_peak_rss_mib']==10 for r in rows))
        self.assertEqual(self.data,before)
        self.assertEqual(benchmark_counts(rows),dict(planned=6,time_available=6,memory_ratio_available=6,
            timeouts=0,other_failed=0,missing_rss_pairs=0))

    def test_missing_and_zero_memory_are_unavailable_not_zero_ratios(self):
        self.data['records'][0]['execution']['original_biogeobears']['peak_rss_bytes']=None
        self.data['records'][1]['execution']['rust_original']['peak_rss_bytes']=0
        rows=self.read();subset=[r for r in rows if r['software']=='original_biogeobears']
        self.assertTrue(all(r['time_ratio']==4 and r['memory_ratio'] is None for r in subset))
        self.assertEqual(subset[1]['rust_peak_rss_bytes'],0)
        self.assertIsNone(subset[1]['rust_peak_rss_mib'])
        self.assertEqual(benchmark_counts(rows)['missing_rss_pairs'],2)
        self.assertEqual(benchmark_counts(rows)['memory_ratio_available'],4)

    def test_failure_and_timeout_keep_observed_partial_memory_without_ratio(self):
        self.data['records'][0]['execution']['revbayes'].update(returncode=1)
        self.data['records'][1]['execution']['rust_revbayes'].update(returncode=-9,timed_out=True)
        rows=self.read();subset=[r for r in rows if r['software']=='revbayes']
        self.assertEqual([r['status'] for r in subset],['failed_or_unsupported','timeout'])
        self.assertTrue(all(r['time_ratio'] is None and r['memory_ratio'] is None for r in subset))
        self.assertTrue(all(r['memory_available'] for r in subset))
        self.assertTrue(all(r['external_peak_rss_mib']==40 for r in subset))
        self.assertEqual(benchmark_counts(rows)['planned'],6)
        self.assertEqual(benchmark_counts(rows)['time_available'],4)
        self.assertEqual(benchmark_counts(rows)['timeouts'],1)

    def test_unaligned_pair_is_not_memory_comparison_even_with_positive_rss(self):
        self.data['records'][0]['comparisons']['lagrange_ng']['alignment_passed']=False
        row=[r for r in self.read() if r['software']=='lagrange_ng'][0]
        self.assertIsNone(row['memory_ratio']);self.assertTrue(row['memory_available'])

    def test_invalid_memory_values_rejected_even_for_failed_process(self):
        for value in [-1,float('nan'),float('inf'),True,'4096']:
            with self.subTest(value=value):
                self.data['records'][0]['execution']['revbayes'].update(returncode=1,peak_rss_bytes=value)
                with self.assertRaisesRegex(ValueError,'peak_rss_bytes'):self.read()

    def test_invalid_time_rejected_and_missing_failed_time_retained(self):
        for value in [0,None,-1,float('nan'),True]:
            with self.subTest(value=value):
                self.data['records'][0]['execution']['revbayes']['wall_seconds']=value
                with self.assertRaises(ValueError):self.read()
        self.data['records'][0]['execution']['revbayes'].update(wall_seconds=None,returncode=1)
        self.assertIsNone([r for r in self.read() if r['software']=='revbayes'][0]['time_ratio'])

    def test_missing_repeat_and_smoke_cannot_enter_formal_table(self):
        saved=copy.deepcopy(self.data)
        self.data['records'].pop()
        with self.assertRaisesRegex(ValueError,'incomplete'):self.read()
        self.data=saved;self.data['formal']=False;self.data['manifest']['mode']='smoke'
        with self.assertRaisesRegex(ValueError,'never relabelled smoke'):self.read()

    def test_contract_keeps_process_scope_font_floor_and_single_figure(self):
        self.assertEqual(CONTRACT['final_size_mm'],[180,225])
        self.assertGreaterEqual(CONTRACT['minimum_font_pt'],8)
        self.assertIn('whole fresh process',CONTRACT['scope'])
        self.assertIn('wait4',CONTRACT['scope'])
        self.assertIn('shared',CONTRACT['legends']['b'])


if __name__=='__main__':unittest.main()
