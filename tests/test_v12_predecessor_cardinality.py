"""Synthetic predecessor ledgers only; no statistical or benchmark execution."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from scripts import benchmark_v12_workflows as bench
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_v12_formal_suite import command_for, launch_priority


def jobs_for(core):
    jobs=[]
    for tips in ([40] if core else [40,150]):
        for areas in ([4] if core else [4,5]):
            for scenario in ('interior','e_zero','j_zero','double_zero'):
                jobs.append(dict(analysis='v1',n_tips=tips,areas=areas,scenario=scenario,
                                 finite_b_extrapolation=False,outer_replicates=200))
    for scenario in (('null_e02','alt_j003') if core else
                     ('null_e0','null_e02','alt_j001','alt_j003','alt_j008')):
        jobs.append(dict(analysis='v2b',n_tips=30,areas=4,scenario=scenario,
                         finite_b_extrapolation=False,outer_replicates=200))
    if not core:
        jobs.append(dict(analysis='v2b',n_tips=30,areas=4,scenario='null_e02',
                         finite_b_extrapolation=True,outer_replicates=200))
    for scenario in ('higher_jump','lower_jump'):
        jobs.append(dict(analysis='v2c',n_tips=30,areas=4,scenario=scenario,
                         finite_b_extrapolation=False,outer_replicates=300))
    jobs.append(dict(analysis='v3',replicates=10 if core else 100))
    if not core:jobs.append(dict(analysis='v4',settings=dict(engineering_fixture=True)))
    jobs.extend([dict(analysis='d1',case_id='Psychotria'),dict(analysis='d2',case_id='Trimezieae')])
    for job in jobs:job['seed']=20261115
    assert len(jobs)==(11 if core else 28)
    return jobs


class PredecessorCardinalityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)

    def fixture(self,core):
        gate=self.folder/'gate.json';bench.write(gate,dict(gate2_passed=True,gate3_passed=True))
        reference=dict(path=gate.name,sha256=bench.sha(gate))
        budget=dict(schema='rustbears.formal_budget.v1',status='locked',jobs=jobs_for(core),
                    gate_evidence=reference,
                    source_sha256={'scripts/run_v12_formal_suite.py':bench.sha(bench.ROOT/'scripts/run_v12_formal_suite.py')})
        path=self.folder/'statistical_budget.json';bench.write(path,budget);digest=bench.sha(path)
        cost=dict(predecessor_budget=dict(path=path.name,sha256=digest),
                  predecessor_statistical_budget_sha256=digest,gate_evidence=reference)
        plan_jobs=[command_for(job,i,self.folder,path,Path('unused'),Path('unused'),
                              4 if core else 8 if job['analysis']=='v1' and job.get('n_tips')==150 and job.get('areas')==5 else 2)
                   for i,job in enumerate(budget['jobs'],1)]
        plan_jobs.sort(key=launch_priority)
        plan=dict(jobs=plan_jobs,budget_sha256=digest,
                  script_sha256=budget['source_sha256']['scripts/run_v12_formal_suite.py'])
        bench.write(self.folder/'suite_plan.json',plan)
        ledger=dict(plan_sha256=bench.sha(self.folder/'suite_plan.json'),
                    status='completed_with_retained_failures',source_unchanged=True,
                    outcomes=[dict(id=j['id'],status='failed_retained' if i==0 else 'processed')
                              for i,j in enumerate(plan_jobs)],
                    cpu_reservations=dict(active_reserved_cpus=0))
        bench.write(self.folder/'suite_ledger.json',ledger)
        return budget,cost,plan,ledger

    def check(self,cost):
        budget,digest=bench.predecessor_budget(cost,self.folder/'cost_budget.json')
        return bench.completed_statistical_suite(self.folder,budget,digest)

    def write_plan(self,plan,ledger):
        bench.write(self.folder/'suite_plan.json',plan)
        ledger=dict(ledger,plan_sha256=bench.sha(self.folder/'suite_plan.json'))
        bench.write(self.folder/'suite_ledger.json',ledger)

    def test_original_28_and_core_11_exact_hashed_predecessors_pass(self):
        for core in (False,True):
            with self.subTest(core=core):
                budget,cost,plan,ledger=self.fixture(core)
                self.assertEqual(len(plan['jobs']),11 if core else 28)
                self.assertEqual(self.check(cost)['status'],'completed_with_retained_failures')
                self.assertEqual({j['reserved_cpus'] for j in plan['jobs']},{4} if core else {2,8})

    def test_missing_extra_or_duplicate_jobs_are_rejected_for_both_sizes(self):
        for core in (False,True):
            for mutation in ('missing','extra','duplicate'):
                with self.subTest(core=core,mutation=mutation):
                    budget,cost,plan,ledger=self.fixture(core)
                    if mutation=='missing':plan['jobs'].pop()
                    else:
                        row=copy.deepcopy(plan['jobs'][0])
                        if mutation=='extra':row['id']='unregistered'
                        plan['jobs'].append(row)
                    self.write_plan(plan,ledger)
                    with self.assertRaises(ValueError):self.check(cost)

    def test_ids_cannot_be_swapped_or_relabelled_even_with_complete_same_settings(self):
        for mutation in ('renamed','swapped','analysis'):
            with self.subTest(mutation=mutation):
                budget,cost,plan,ledger=self.fixture(True)
                if mutation=='renamed':
                    old=plan['jobs'][0]['id'];plan['jobs'][0]['id']='foreign_id'
                    for row in ledger['outcomes']:
                        if row['id']==old:row['id']='foreign_id'
                elif mutation=='swapped':
                    plan['jobs'][0]['settings'],plan['jobs'][1]['settings']=plan['jobs'][1]['settings'],plan['jobs'][0]['settings']
                else:plan['jobs'][0]['analysis']='v4'
                self.write_plan(plan,ledger)
                with self.assertRaisesRegex(ValueError,'IDs/configurations'):self.check(cost)

    def test_configuration_budget_and_recorded_launcher_drift_rejected(self):
        for mutation in ('settings','budget_hash','source_hash'):
            with self.subTest(mutation=mutation):
                budget,cost,plan,ledger=self.fixture(True)
                if mutation=='settings':plan['jobs'][0]['settings']['seed']+=1
                elif mutation=='budget_hash':plan['budget_sha256']='wrong'
                else:plan['script_sha256']='wrong'
                self.write_plan(plan,ledger)
                with self.assertRaises(ValueError):self.check(cost)

    def test_actual_launcher_source_drift_rejected_before_using_constructor(self):
        budget,cost,plan,ledger=self.fixture(True)
        source=self.folder/'scripts/run_v12_formal_suite.py';source.parent.mkdir()
        source.write_text('# changed source; never imported\n')
        with patch.object(bench,'ROOT',self.folder),self.assertRaisesRegex(ValueError,'launcher source'):
            self.check(cost)

    def test_missing_duplicate_nonterminal_or_busy_ledger_is_rejected(self):
        for core in (False,True):
            for mutation in ('missing','duplicate','running','busy','unstable','plan_hash'):
                with self.subTest(core=core,mutation=mutation):
                    budget,cost,plan,ledger=self.fixture(core)
                    if mutation=='missing':ledger['outcomes'].pop()
                    elif mutation=='duplicate':ledger['outcomes'][-1]=copy.deepcopy(ledger['outcomes'][0])
                    elif mutation=='running':ledger['outcomes'][0]['status']='running'
                    elif mutation=='busy':ledger['cpu_reservations']['active_reserved_cpus']=4
                    elif mutation=='unstable':ledger['source_unchanged']=False
                    else:ledger['plan_sha256']='wrong'
                    bench.write(self.folder/'suite_ledger.json',ledger)
                    with self.assertRaises(ValueError):self.check(cost)

    def test_upstream_budget_bytes_and_same_gate_binding_remain_required(self):
        budget,cost,plan,ledger=self.fixture(True)
        changed=copy.deepcopy(cost);changed['gate_evidence']['sha256']='different gate'
        with self.assertRaisesRegex(ValueError,'same passed gates'):self.check(changed)
        budget['jobs'][0]['seed']+=1;bench.write(self.folder/'statistical_budget.json',budget)
        with self.assertRaisesRegex(ValueError,'original statistical budget'):self.check(cost)


if __name__=='__main__':unittest.main()
