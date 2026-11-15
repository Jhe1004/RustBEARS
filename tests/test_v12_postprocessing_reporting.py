"""Postprocessing semantics only: fixtures never render manuscript data."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from scripts import plot_v12_performance as plot
from scripts import v12_workflow_costs as costs


class ProfileReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.records=[dict(param=p,interval_kind='profile_confidence_interval',mle=.1,lower=0,upper=.3,status_lower='at_physical_zero',
            status_upper='interior_bracketed',disconnected_support_detected='False',
            endpoint_deviance_residual_max=1e-8,se_hessian_local=.2,onesided_upper_from_zero=.25,
            critical_value=3.84) for p in ('d','e','j')]
        self.records[1].update(upper='',status_upper='censored_at_box_upper',se_hessian_local='')
        self.records[2].update(disconnected_support_detected='True')
        self.write_intervals()
        (self.folder/'profile_result.json').write_text(json.dumps(dict(rebased_mle=True)))
        self.evidence={str(p):costs.sha(p) for p in self.folder.iterdir()}

    def write_intervals(self):
        with (self.folder/'intervals.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(self.records[0]))
            writer.writeheader();writer.writerows(self.records)

    def read(self):return costs.profile_cost_diagnostics(self.folder,self.evidence,'returned')

    def test_returned_cost_keeps_distinct_interval_diagnostics(self):
        result=self.read()
        self.assertEqual(result['profile_box_censored_endpoints'],1)
        self.assertEqual(result['profile_disconnected_parameters'],1)
        self.assertEqual(result['profile_missing_se_parameters'],1)
        self.assertTrue(result['profile_rebased_mle'])
        e=result['profile_endpoint_diagnostics'][1]
        self.assertIsNone(e['upper']);self.assertIsNone(e['se_hessian_local'])
        self.assertTrue(e['se_hessian_missing'])
        self.assertNotIn('n_failed',result)
        self.assertEqual(result['profile_artifacts'],self.evidence)

    def test_changed_or_unhashed_diagnostics_are_rejected(self):
        for name in ('intervals.csv','profile_result.json'):
            missing=dict(self.evidence);del missing[str(self.folder/name)]
            with self.subTest(name=name),self.assertRaisesRegex(ValueError,'hash-bound'):
                costs.profile_cost_diagnostics(self.folder,missing,'returned')
        self.records[0]['status_upper']='numerical_failure';self.write_intervals()
        with self.assertRaisesRegex(ValueError,'hash-bound'):self.read()

    def test_unknown_or_uncertified_status_cannot_be_silently_returned(self):
        for status in ('numerical_failure','numerically_uncertified','new_unreviewed_status'):
            self.records[0]['status_upper']=status;self.write_intervals()
            self.evidence[str(self.folder/'intervals.csv')]=costs.sha(self.folder/'intervals.csv')
            with self.subTest(status=status),self.assertRaisesRegex(ValueError,'uncertified returned-profile'):
                self.read()

    def test_raised_call_cannot_acquire_returned_diagnostics(self):
        with self.assertRaisesRegex(ValueError,'returned API'):
            costs.profile_cost_diagnostics(self.folder,self.evidence,'raised')

    def test_missing_endpoint_must_be_explicitly_censored(self):
        self.records[0]['upper']='';self.write_intervals()
        self.evidence[str(self.folder/'intervals.csv')]=costs.sha(self.folder/'intervals.csv')
        with self.assertRaisesRegex(ValueError,'endpoint/status mismatch'):self.read()

    def test_diagnostics_export_preserves_na_and_api_status(self):
        fields=self.read();out=self.folder/'source';out.mkdir()
        row=dict(fields,workload='profile_chi2',job_id='example',task_id='example:profile_chi2',
            row_id='example:profile_chi2:1',states=16,workers=1,replicate=1,
            status='returned',process_status='complete',api_complete_fresh=True,wall_seconds=2.)
        counts=[dict(workload='profile_chi2',states=256,workers=8,unique_tasks_observed=0,
                     planned_tasks=3,full_attempts=0,unavailable_tasks=3)]
        audit=dict(status='eligible_formal_workflow_costs',formal=True,rows=[row],points=[row],
                   jobs=[],cost_tasks=[],stratum_counts=counts)
        costs.export_costs(audit,out)
        with (out/'workflow_profile_diagnostics.csv').open() as stream:records=list(csv.DictReader(stream))
        self.assertEqual(len(records),3)
        self.assertEqual(records[1]['status'],'returned')
        self.assertEqual(records[1]['upper'],'');self.assertEqual(records[1]['se_hessian_local'],'')
        self.assertEqual(records[1]['status_upper'],'censored_at_box_upper')
        with (out/'workflow_complete_costs.csv').open() as stream:record=next(csv.DictReader(stream))
        self.assertEqual(json.loads(record['profile_endpoint_diagnostics'])[1]['upper'],None)
        self.assertEqual(record['n_failed'],'')
        with (out/'workflow_stratum_counts.csv').open() as stream:record=next(csv.DictReader(stream))
        self.assertEqual((record['states'],record['workers'],record['unique_tasks_observed'],record['planned_tasks']),('256','8','0','3'))
        self.assertFalse(any(p.suffix in ('.pdf','.svg','.png','.tiff') for p in out.iterdir()))


class PlannedLegendTests(unittest.TestCase):
    def setUp(self):
        self.settings=dict(states=[16,64,256],workers=[1,8],repeats=3)
        self.workloads=('profile_chi2','profile_plugin')
        self.tasks=[dict(task_id=f'{phase}:{state}:{workers}:{rep}',workload=phase,states=state,
            workers=workers,replicate=rep,full_cost_attempts=0)
            for phase in self.workloads for state in self.settings['states']
            for workers in self.settings['workers'] for rep in range(1,4)]

    def audit(self,points=()):
        tasks=copy.deepcopy(self.tasks)
        for task in tasks:task['full_cost_attempts']=sum(p['task_id']==task['task_id'] for p in points)
        return dict(settings=self.settings,points=list(points),
            stratum_counts=costs.dedicated_stratum_counts(self.settings,tasks,points,self.workloads))

    def point(self,index=0):return dict(self.tasks[index],wall_seconds=1.)

    def test_wholly_missing_states_and_workers_keep_encodings_and_denominators(self):
        audit=self.audit([self.point()]);legend=plot.workflow_legend_spec(audit)
        full=plot.workflow_legend_spec(self.audit([dict(t,wall_seconds=1.) for t in self.tasks]))
        self.assertEqual(legend['state_colors'],full['state_colors'])
        self.assertEqual(legend['cpu_markers'],full['cpu_markers'])
        self.assertEqual(list(legend['state_colors']),[16,64,256])
        self.assertEqual(legend['state_labels'][256],'256 states (0/12)')
        self.assertEqual(legend['cpu_labels'][8],'8 CPUs reserved (0/18)')
        missing=[r for r in audit['stratum_counts'] if r['states']==256]
        self.assertEqual(len(missing),4)
        self.assertTrue(all(r['unique_tasks_observed']==0 and r['planned_tasks']==3 for r in missing))
        self.assertEqual(len(audit['points']),1)

    def test_no_full_observations_do_not_create_zero_duration_points(self):
        audit=self.audit();legend=plot.workflow_legend_spec(audit)
        self.assertEqual(audit['points'],[])
        self.assertEqual(len(audit['stratum_counts']),12)
        self.assertTrue(all(r['planned_tasks']==3 and r['unique_tasks_observed']==0 for r in audit['stratum_counts']))
        self.assertTrue(all('(0/' in label for label in legend['state_labels'].values()))

    def test_repeated_attempts_are_not_multiple_observed_tasks(self):
        audit=self.audit([self.point(),dict(self.point(),wall_seconds=2.)])
        first=audit['stratum_counts'][0]
        self.assertEqual((first['full_attempts'],first['unique_tasks_observed'],first['planned_tasks']),(2,1,3))

    def test_planned_order_does_not_follow_observation_order(self):
        self.settings['states']=[256,16,64];self.settings['workers']=[8,1]
        legend=plot.workflow_legend_spec(self.audit([self.point()]))
        self.assertEqual(list(legend['state_colors']),[256,16,64])
        self.assertEqual(list(legend['cpu_markers']),[8,1])

    def test_missing_or_duplicate_tasks_are_not_zero_success_strata(self):
        for tasks in (self.tasks[:-1],self.tasks+[self.tasks[0]]):
            with self.assertRaisesRegex(ValueError,'planned cost strata'):
                costs.dedicated_stratum_counts(self.settings,tasks,[],self.workloads)

    def test_unplanned_or_mislabelled_point_is_rejected(self):
        for change in (dict(states=256),dict(task_id='unknown')):
            with self.assertRaisesRegex(ValueError,'planned stratum'):
                self.audit([dict(self.point(),**change)])


if __name__=='__main__':unittest.main()
