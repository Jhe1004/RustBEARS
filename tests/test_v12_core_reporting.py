"""Reporting/provenance fixtures only; never execute or render science."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import plot_v12_calibration as plot
from scripts import plot_v12_performance as performance
from scripts import v12_core_reporting as core
from scripts import v12_workflow_costs as costs
from scripts import benchmark_v12_workflows as bench


def design_fixture():
    jobs=[dict(analysis='v1',scenario=s,n_tips=40,areas=4,outer_replicates=200,
               plugin_bootstrap=199,seed=20261115) for s in plot.V1_SCENARIOS]
    jobs += [dict(analysis='v2b',scenario=s,n_tips=30,areas=4,outer_replicates=200,
                  inner_bootstrap=99,surrogate_alt_replicates=0,finite_b_extrapolation=False,seed=20261115)
             for s in ('null_e02','alt_j003')]
    jobs += [dict(analysis='v2c',scenario=s,n_tips=30,areas=4,outer_replicates=300,seed=20261115)
             for s in ('higher_jump','lower_jump')]
    jobs += [dict(analysis='v3',replicates=10,seed=20261115)]
    jobs += [dict(analysis=a,case_id=c,seed=20261115) for a,c in [('d1','D1_Psychotria'),('d2','D2_Trimezieae')]]
    registry=[dict(id=f'{i:02d}_{j["analysis"]}',settings_index=i-1,analysis=j['analysis'],
                   expected_task_ids=core.planned_ids(j)) for i,j in enumerate(jobs,1)]
    old=dict(analysis='v4',seed=20261115,settings={s:dict(replicates=2) for s in
             ('fixed','fixed_trees','independent','mcmc','routes','adaptive')})
    return dict(schema='rustbears.core_statistical_design.v1',scope='core_v1',source_commit=None,
        specified_utc='2019-01-01T00:00:00+00:00',jobs=jobs,job_registry=registry,
        sampling_policy=dict(v3_dataset_ids=list(range(1,11))),
        inherited_evidence=[dict(id='26_v4',analysis='v4',settings=old,fresh_execution=False)])


class CoreReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name).resolve();self.design=design_fixture()

    def write(self,name,value):
        path=self.root/name;path.parent.mkdir(parents=True,exist_ok=True)
        plot.write_json(path,value);return path

    def ref(self,path):return dict(path=str(path),sha256=plot.sha(path))

    def test_registry_retains_twelve_jobs_and_exact_fixed_ids(self):
        registry=core.verify_registry(self.design)
        self.assertEqual(len(registry),12)
        self.assertEqual(sum(v['origin']=='fresh' for v in registry.values()),11)
        self.assertEqual(registry['09_v3']['expected_task_ids']['ids'],list(range(1,11)))
        for change in ('missing','duplicate','favorable_ids','wrong_scenario'):
            value=copy.deepcopy(self.design)
            if change=='missing':value['job_registry'].pop()
            elif change=='duplicate':value['job_registry'][-1]['id']=value['job_registry'][0]['id']
            elif change=='favorable_ids':value['sampling_policy']['v3_dataset_ids']=list(range(2,12))
            else:value['jobs'][4]['scenario']='null_e05'
            with self.subTest(change=change),self.assertRaises(ValueError):core.verify_registry(value)

    def test_expected_ids_include_failures_and_reject_count_only_evidence(self):
        expected=dict(kind='dataset',ids=[1,2],count=2)
        core.actual_ids(dict(analysis='v3',records=[dict(dataset=1,status='failed_retained'),dict(dataset=2,status='complete')]),expected)
        for ids in ([1],[1,1],[2,3]):
            with self.assertRaisesRegex(ValueError,'missing or duplicated'):
                core.actual_ids(dict(analysis='v3',records=[dict(dataset=i) for i in ids]),expected)
        records=[dict(section=s,replicate=i) for s in ('fixed','adaptive') for i in (1,2)]
        core.actual_ids(dict(analysis='v4',records=records),dict(fixed=[1,2],adaptive=[1,2]))
        records[-1]['replicate']=1
        with self.assertRaisesRegex(ValueError,'V4 task IDs'):core.actual_ids(dict(analysis='v4',records=records),dict(fixed=[1,2],adaptive=[1,2]))

    def budget(self,name,design,commit):
        d=self.write(name+'/design.json',design)
        g=self.write(name+'/gate.json',dict(gate2_passed=True,gate3_passed=True))
        source=self.root/name/'source';source.mkdir()
        names=['rustbears/__init__.py',*bench.SOURCE_FILES,'scripts/run_v12_formal_suite.py','scripts/validate_v12_v4_protocol.py','native.so']
        for filename in names:
            p=source/filename;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(name+filename)
        hashes={n:plot.sha(source/n) for n in names}
        plot.write_json(g,dict(gate2_passed=True,gate3_passed=True,source_commit=commit,native_sha256=hashes['native.so']))
        payload=dict(schema='rustbears.formal_budget.v1',status='locked',locked_utc='2020-01-01T00:00:00+00:00',
            gate_evidence=self.ref(g),design_evidence=self.ref(d),jobs=design['jobs'],source_commit=commit,
            source_sha256=hashes,native_sha256=hashes['native.so'])
        path=self.write(name+'/budget.json',payload)
        return payload,path,d,dict(root=str(source),commit=commit,sha256=hashes)

    def collection(self):
        old=self.design['inherited_evidence'][0]
        ob,op,od,osource=self.budget('old',dict(jobs=[old['settings']]),'b'*40)
        old.update(original_budget=self.ref(op),original_source_commit='b'*40,original_native_sha256=ob['native_sha256'])
        budget,bpath,dpath,source=self.budget('new',self.design,'a'*40)
        inv=dict(schema='rustbears.figure2.inputs.v2',formal=True,scope='core_v1',core_design=self.ref(dpath),
                 panel_selection=core.core_panels(),runs=[])
        results={}
        for key,value in core.verify_registry(self.design).items():
            job=value['job'];analysis=job['analysis'];isold=value['origin']=='inherited'
            b,bp,dp,src=(ob,op,od,osource) if isold else (budget,bpath,dpath,source)
            folder=self.root/key;folder.mkdir()
            item=dict(id=key,analysis=analysis,origin=value['origin'],path=str(folder),budget=self.ref(bp),design=self.ref(dp),
                      source=src,expected_task_ids=value['expected_task_ids'])
            metadata=dict(payload=b,sha256=plot.sha(bp),locked_utc=b['locked_utc'])
            manifest=dict(job,mode='formal',source_sha256={'rustbears/__init__.py':src['sha256']['rustbears/__init__.py']},locked_budget=metadata)
            if analysis=='v3':
                manifest=dict(settings=job,budget=metadata,production_source_sha256={'__init__.py':src['sha256']['rustbears/__init__.py']},native_sha256=b['native_sha256'])
                plan=self.write(key+'/suite_plan.json',dict(budget_sha256=plot.sha(bp),jobs=[dict(id=key,settings=job,output=str(folder))]))
                ledger=self.write(key+'/suite_ledger.json',dict(plan_sha256=plot.sha(plan),outcomes=[dict(id=key,status='processed',started_utc='2021-01-01T00:00:00+00:00')]))
                item['coordinator_evidence']=dict(plan=self.ref(plan),ledger=self.ref(ledger))
                rows=[dict(dataset=i,status='failed_retained') for i in range(1,11)]
            elif analysis=='v4':
                manifest.update(execution_settings=job['settings'],source_sha256={
                    'scripts/validate_v12_v4_protocol.py':src['sha256']['scripts/validate_v12_v4_protocol.py']})
                rows=[dict(section=s,replicate=i) for s in job['settings'] for i in (1,2)]
            else:
                rows=[dict(id=i,status='failed_retained') for i in range(1,job.get('outer_replicates',0)+1)]
                manifest['script_sha256']=src['sha256']['scripts/run_uncertainty_validation.py']
            start=self.write(key+'/execution_started.json',dict(first_started_utc='2021-01-01T00:00:00+00:00'))
            mp=self.write(key+'/manifest.json',manifest)
            results[str(folder)]=dict(analysis=analysis,path=folder,manifest=manifest,records=rows,summary={},files={str(mp):plot.sha(mp),str(start):plot.sha(start)})
            inv['runs'].append(item)
        path=self.write('collection.json',inv)
        return inv,path,results

    def inspect_fixture(self,path,results):
        seen=[]
        def read_run(item,base,digest):
            result=copy.deepcopy(results[item['path']]);seen.append((result['analysis'],digest));return result
        def case(folder,*args):
            return dict(status='incomplete_profiles',files=results[str(folder)]['files'])
        with patch.object(plot,'load_run',side_effect=read_run),patch('scripts.v12_case_evidence.verify_case',side_effect=case):
            result=core.inspect_collection(path)
        return result,seen

    def test_collection_passes_each_original_budget_to_unchanged_scientific_reader(self):
        inv,path,results=self.collection();(runs,audit),seen=self.inspect_fixture(path,results)
        self.assertEqual((len(runs),audit['jobs'],len(audit['cases'])),(10,12,2))
        old=next(r for r in inv['runs'] if r['analysis']=='v4')
        self.assertEqual(next(h for a,h in seen if a=='v4'),old['budget']['sha256'])
        self.assertNotEqual(seen[0][1],old['budget']['sha256'])
        self.assertEqual(runs[-1]['origin'],'inherited')
        self.assertEqual(audit['scope'],'core_v1')

    def test_collection_rejects_old_budget_relabel_missing_job_and_changed_source(self):
        inv,path,results=self.collection()
        for change in ('old_budget','missing','duplicate','source','scope'):
            value=copy.deepcopy(inv)
            if change=='old_budget':value['runs'][-1]['budget']=value['runs'][0]['budget']
            elif change=='missing':value['runs'].pop()
            elif change=='duplicate':value['runs'][-1]=value['runs'][0]
            elif change=='source':value['runs'][0]['source']['commit']='f'*40
            else:value['scope']='full_v1'
            self.write('collection.json',value)
            with self.subTest(change=change),self.assertRaises(ValueError):self.inspect_fixture(path,results)

    def test_actual_source_drift_and_design_after_lock_fail_closed(self):
        inv,path,results=self.collection();item=inv['runs'][0]
        root=Path(item['source']['root']);(root/'rustbears/__init__.py').write_text('drift')
        with self.assertRaisesRegex(ValueError,'execution source changed'):self.inspect_fixture(path,results)
        b=plot.read_json(item['budget']['path']);design_path=Path(item['design']['path']);d=plot.read_json(design_path)
        d['specified_utc']='2021-01-01T00:00:00+00:00';plot.write_json(design_path,d)
        b['design_evidence']=self.ref(design_path);plot.write_json(Path(item['budget']['path']),b)
        with self.assertRaisesRegex(ValueError,'design postdates'):
            core.bound_budget(self.ref(Path(item['budget']['path'])),self.root,{})

    def test_full_default_rejects_core_budget_without_explicit_v2_scope(self):
        inv,path,_=self.collection()
        legacy=dict(schema='rustbears.figure2.inputs.v1',formal=True,budget=inv['runs'][0]['budget'],runs=[])
        with self.assertRaisesRegex(ValueError,'all V1/V2b/V2c/V3/V4'):plot.budget_info(legacy,self.root)

    def test_gate_failure_cannot_be_hidden_by_updated_reference_hash(self):
        budget,bpath,_,_=self.budget('badgate',self.design,'a'*40)
        gate=Path(budget['gate_evidence']['path']);plot.write_json(gate,dict(gate2_passed=True,gate3_passed=False))
        budget['gate_evidence']=self.ref(gate);plot.write_json(bpath,budget)
        with self.assertRaisesRegex(ValueError,'Gate 2 and Gate 3'):core.bound_budget(self.ref(bpath),self.root,{})

    def test_passed_gate_for_another_source_or_native_is_rejected(self):
        budget,bpath,_,_=self.budget('wrongsourcegate',self.design,'a'*40)
        gate=Path(budget['gate_evidence']['path']);original=plot.read_json(gate)
        for change in (dict(source_commit='b'*40),dict(native_sha256='different'),dict(source_sha256={'rustbears/__init__.py':'different'})):
            plot.write_json(gate,dict(original,**change));budget['gate_evidence']=self.ref(gate);plot.write_json(bpath,budget)
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,'gate source'):
                core.bound_budget(self.ref(bpath),self.root,{})

    def test_v3_coordinator_cannot_claim_a_start_before_the_budget_lock(self):
        inv,path,results=self.collection();item=next(r for r in inv['runs'] if r['analysis']=='v3')
        ledger_path=Path(item['coordinator_evidence']['ledger']['path']);ledger=plot.read_json(ledger_path)
        ledger['outcomes'][0]['started_utc']='2019-06-01T00:00:00+00:00';plot.write_json(ledger_path,ledger)
        item['coordinator_evidence']['ledger']=self.ref(ledger_path);self.write('collection.json',inv)
        with self.assertRaisesRegex(ValueError,'before its budget lock'):self.inspect_fixture(path,results)

    def predecessor(self):
        registry=core.verify_registry(self.design);jobs=[dict(id=k,settings=v['job']) for k,v in registry.items() if v['origin']=='fresh']
        budget=dict(jobs=self.design['jobs'],source_sha256={'scripts/run_v12_formal_suite.py':'launcher'})
        plan=self.write('suite_plan.json',dict(budget_sha256='budget',script_sha256='launcher',jobs=jobs))
        ledger=dict(plan_sha256=plot.sha(plan),source_unchanged=True,status='completed_with_retained_failures',
            outcomes=[dict(id=j['id'],status='failed_retained') for j in jobs],cpu_reservations=dict(active_reserved_cpus=0))
        self.write('suite_ledger.json',ledger);return registry,budget,ledger

    def test_predecessor_accepts_eleven_terminal_retained_failures_without_success_claim(self):
        registry,budget,ledger=self.predecessor()
        out=core.core_predecessor(self.root,budget,'budget',registry,{})
        self.assertEqual(out['status'],'completed_with_retained_failures')
        for change in ('duplicate','missing','active','running','false_processed'):
            value=copy.deepcopy(ledger)
            if change=='duplicate':value['outcomes'][-1]=value['outcomes'][0]
            elif change=='missing':value['outcomes'].pop()
            elif change=='active':value['cpu_reservations']['active_reserved_cpus']=4
            elif change=='running':value['outcomes'][0]['status']='running'
            else:value['outcomes'][0].update(status='processed',result={})
            self.write('suite_ledger.json',value)
            with self.subTest(change=change),self.assertRaises(ValueError):core.core_predecessor(self.root,budget,'budget',registry,{})

    def test_predecessor_processed_evidence_drift_rejected(self):
        registry,budget,ledger=self.predecessor();artifact=self.write('result.json',dict(status='retained'))
        ledger['outcomes'][0].update(status='processed',result=dict(path=str(artifact),sha256=plot.sha(artifact),evidence_files={str(artifact):plot.sha(artifact)}))
        self.write('suite_ledger.json',ledger);core.core_predecessor(self.root,budget,'budget',registry,{})
        artifact.write_text('{}')
        with self.assertRaisesRegex(ValueError,'evidence changed'):core.core_predecessor(self.root,budget,'budget',registry,{})

    def cost_inventory(self):
        collection,_,results=self.collection();fresh=collection['runs'][0]
        predecessor_path=Path(fresh['budget']['path']);predecessor=plot.read_json(predecessor_path)
        settings=dict(bench.DEFAULT,repeats=1);cells=bench.cells_for(settings)
        workflow=dict(schema='rustbears.core_workflow_cost_design.v1',scope='core_v1',
            source_commit=None,specified_utc=self.design['specified_utc'],jobs=[dict(analysis='workflow_benchmark',**settings)],
            expected_cell_ids=[c['id'] for c in cells],expected_phases=list(bench.PHASES),planned_tasks=36)
        workflow_path=self.write('workflow_design.json',workflow)
        budget=dict(predecessor,jobs=workflow['jobs'],design_evidence=self.ref(workflow_path),
            predecessor_budget=self.ref(predecessor_path),predecessor_statistical_budget_sha256=plot.sha(predecessor_path))
        budget_path=self.write('cost_budget.json',budget)
        registry=core.verify_registry(self.design);jobs=[dict(id=k,settings=v['job']) for k,v in registry.items() if v['origin']=='fresh']
        plan=self.write('coordinator/suite_plan.json',dict(budget_sha256=plot.sha(predecessor_path),
            script_sha256=predecessor['source_sha256']['scripts/run_v12_formal_suite.py'],jobs=jobs))
        self.write('coordinator/suite_ledger.json',dict(plan_sha256=plot.sha(plan),source_unchanged=True,
            status='completed_with_retained_failures',cpu_reservations=dict(active_reserved_cpus=0),
            outcomes=[dict(id=j['id'],status='failed_retained') for j in jobs]))
        self.write('costs/plan.json',dict(after_suite=str(plan.parent)))
        inv=dict(schema='rustbears.workflow_benchmark.inputs.v2',formal=True,scope='core_v1',
            budget=self.ref(budget_path),workflow_design=self.ref(workflow_path),core_design=collection['core_design'],
            source=fresh['source'],benchmark=str(self.root/'costs'),
            expected_task_ids=[c['id']+':'+p for c in cells for p in bench.PHASES],
            workload_labels=['χ² profiles: d, e, j','Plug-in profiles: B=199','LRT: B=199',
                             'Six models + average','5000 single-tree maps','20 trees × 100 maps'],
            inherited_v4=collection['runs'][-1])
        path=self.write('cost_inventory.json',inv)
        return inv,path,results

    def test_core_cost_context_binds_separate_budget_payload_and_inherited_v4(self):
        inv,path,results=self.cost_inventory()
        with patch.object(plot,'load_run',side_effect=lambda item,*_:copy.deepcopy(results[item['path']])):
            context=core.core_cost_context(path,inv,bench)
        self.assertEqual(context['scope'],'core_v1')
        self.assertIn('One planned repetition',context['repetition_interpretation'])
        self.assertEqual(context['labels'][2],'LRT: B=199')
        self.assertIn(inv['inherited_v4']['budget']['path'],context['files'])
        self.assertEqual(context['predecessor_checked']['status'],'completed_with_retained_failures')

    def test_core_cost_rejects_missing_task_labels_and_inherited_identity_replacement(self):
        inv,path,results=self.cost_inventory()
        for change in ('task','duplicate','label','inherited_budget'):
            value=copy.deepcopy(inv)
            if change=='task':value['expected_task_ids'].pop()
            elif change=='duplicate':value['expected_task_ids'][-1]=value['expected_task_ids'][0]
            elif change=='label':value['workload_labels'][2]='LRT: 199 + 199'
            else:value['inherited_v4']['budget']=value['budget']
            self.write('cost_inventory.json',value)
            with patch.object(plot,'load_run',side_effect=lambda item,*_:copy.deepcopy(results[item['path']])),self.subTest(change=change),self.assertRaises(ValueError):
                core.core_cost_context(path,value,bench)

    def test_complete_core_cost_reader_retains_all_failed_attempts_without_full_cost_points(self):
        inv,path,old_runs=self.cost_inventory();budget=plot.read_json(inv['budget']['path'])
        settings=dict(bench.DEFAULT,repeats=1);root=Path(inv['source']['root']);folder=Path(inv['benchmark'])
        source=dict(files={n:plot.sha(root/n) for n in [*bench.SOURCE_FILES,'rustbears/__init__.py']},
            native_path=str(root/'native.so'),native_sha256=plot.sha(root/'native.so'))
        plan=dict(mode='formal',formal=True,budget=dict(sha256=inv['budget']['sha256'],payload=budget),
            settings=settings,cells=bench.cells_for(settings),phases=list(bench.PHASES),sources=source,
            after_suite=str(self.root/'coordinator'),host='engineering-fixture',full_proposal_grid=False)
        plan_path=self.write('costs/plan.json',plan);records=[]
        for cell in plan['cells']:
            self.write('costs/'+cell['id']+'/input_manifest.json',dict(files={}))
            for phase in bench.PHASES:
                directory=folder/cell['id']/phase/'attempt_001';directory.mkdir(parents=True)
                plot.write_json(directory/'spec.json',dict(mode='formal',cell=cell,settings=settings,sources=source,
                    input_sha256={},started_utc='2021-01-01T00:00:00+00:00'))
                (directory/'stderr.log').write_text('Retained engineering failure. No scientific execution occurred.\n')
                files={str(p):plot.sha(p) for p in directory.iterdir()}
                attempt=dict(directory=str(directory),status='failed_retained',returncode=1,timed_out=False,
                    interrupted=False,process_wall_seconds=.1,process_peak_rss_bytes=1,api_right_censored_elapsed_seconds=None)
                records.append(dict(cell_id=cell['id'],phase=phase,status='failed_retained',attempts=[attempt],files=files))
        with patch.object(plot,'load_run',side_effect=lambda item,*_:copy.deepcopy(old_runs[item['path']])):
            context=core.core_cost_context(path,inv,bench)
            self.write('costs/results.json',dict(formal=True,status='all_planned_phases_processed',source_unchanged=True,
                plan_sha256=plot.sha(plan_path),first_started_utc='2021-01-01T00:00:00+00:00',
                after_suite=context['predecessor_checked'],records=records))
            audit=costs.inspect_dedicated_benchmark(path)
        self.assertEqual(audit['scope'],'core_v1');self.assertEqual(audit['independent_repetitions'],1)
        self.assertFalse(audit['between_repeat_variability_estimated'])
        self.assertEqual(audit['points'],[]);self.assertEqual(len(audit['rows']),36)
        self.assertEqual(len(audit['stratum_counts']),36)
        self.assertTrue(all(r['wall_seconds'] is None for r in audit['rows']))
        self.assertTrue(all(r['planned_tasks']==1 and r['unique_tasks_observed']==0 for r in audit['stratum_counts']))


class CoreTableTests(unittest.TestCase):
    def fixture(self):
        spec=importlib.util.spec_from_file_location('calibration_fixtures',Path(__file__).with_name('test_v12_calibration_plot.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        fixture=module.CalibrationPlotTests();fixture.path=Path('/nonexistent/arithmetic_fixture')
        return [r for r in fixture.fixture_runs() if not (r['analysis']=='v1' and r['job']['n_tips']!=40)
                and not (r['analysis']=='v2b' and r['job']['scenario'] not in ('null_e02','alt_j003'))]

    def test_core_fixed_panel_ids_preserve_missingness_and_failure_bounds(self):
        runs=self.fixture();tables=plot.main_tables(runs,selection=core.core_panels())
        self.assertEqual([len(tables[k]) for k in 'abcd'],[16,6,3,4])
        self.assertTrue(all(r['estimate'] is None for r in tables['a'] if r['method']=='wald_local'))
        self.assertTrue(all(r['failure_bounds']==[0,1] for r in tables['b']))
        self.assertTrue(all(r['estimate'] is None for r in tables['c']))
        with self.assertRaisesRegex(ValueError,'incomplete or duplicated'):plot.main_tables(runs)
        runs[1]=copy.deepcopy(runs[0])
        with self.assertRaisesRegex(ValueError,'incomplete or duplicated'):plot.main_tables(runs,selection=core.core_panels())

    def test_one_repetition_preserves_all_thirtysix_planned_cells_and_zero_of_n(self):
        settings=dict(states=[16,64,256],workers=[1,8],repeats=1)
        tasks=[dict(task_id=f'{phase}:{state}:{worker}:1',workload=phase,states=state,workers=worker,replicate=1,full_cost_attempts=0)
               for phase in bench.PHASES for state in settings['states'] for worker in settings['workers']]
        counts=costs.dedicated_stratum_counts(settings,tasks,[],bench.PHASES)
        self.assertEqual(len(counts),36);self.assertTrue(all(r['planned_tasks']==1 and r['unique_tasks_observed']==0 for r in counts))
        legend=performance.workflow_legend_spec(dict(settings=settings,points=[],stratum_counts=counts))
        self.assertEqual(legend['state_labels'][256],'256 states (0/12)')
        self.assertEqual(legend['cpu_labels'][8],'8 CPUs reserved (0/18)')


if __name__=='__main__':unittest.main()
