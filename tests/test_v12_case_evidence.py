"""Tiny artifact fixtures verify retained failures and missing-row rejection.

These files test the evidence parser; they never represent experimental data.
"""
import copy
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from v12_case_evidence import verify_case_content,MODELS,CONTRACTS,E4
from plot_v12_calibration import decode
from rustbears.uncertainty.checkpoint import canonical_manifest,task_fingerprint


class CaseEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name)
        self.manifest=dict(models=MODELS,likelihood_contracts=CONTRACTS,tip_count=3,bsm_maps=2,
            bsm_maps_per_tree=2,n_boot=2,max_trees=2,seed=1,original_posterior_required=True,
            model_parameter_bounds={m:{'d':[0,5],'e':[0,5],**({'j':[0,3]} if m.endswith('+J') else {})} for m in MODELS})
        self.report=dict(status='incomplete_profiles',contracts={})
        self.make_fixture()

    def dump(self,path,value):path.write_text(json.dumps(value))

    def csv(self,path,rows,fields=None):
        with path.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields or list(rows[0]));writer.writeheader();writer.writerows(rows)

    def store(self,path,manifest,records):
        header=dict(schema='rustbears.checkpoint.v1',task_fingerprint=task_fingerprint(manifest),manifest=canonical_manifest(manifest))
        path.write_text('\n'.join(json.dumps(r) for r in [header,*[dict(record=r) for r in records]])+'\n')

    def maps(self,directory,prefix):
        rows=[dict(map=i,map_id=i,tree_id=1,model_id=1,tree_topology_hash='unit',stratum_weight=1.,
                   precision_mode='strict_f64',likelihood_conditioning=directory.name) for i in (1,2)]
        self.csv(directory/(prefix+'map_metadata.csv'),rows)
        events=[dict(row,node=node) for row in rows for node in (4,5)]
        self.csv(directory/(prefix+'cladogenetic_events.csv'),events)
        self.csv(directory/(prefix+'anagenetic_events.csv'),[],fields=list(events[0]))
        for name in ('event_precision','route_precision','route_comparisons','event_through_time'):
            self.csv(directory/(prefix+name+'.csv'),[dict(status='unit_fixture')])
        summary=dict(status='complete',maps=2,precision_target_met=False)
        self.dump(directory/(prefix+'bsm_summary.json'),summary)
        return summary

    def make_fixture(self):
        for contract in CONTRACTS:
            directory=self.path/contract;directory.mkdir()
            fits=[dict(model=m,converged=True,lnL=-1.,parameter_bounds=self.manifest['model_parameter_bounds'][m]) for m in MODELS]
            self.dump(directory/'initial_fitted_parameters.json',fits);self.dump(directory/'fitted_parameters.json',fits)
            self.csv(directory/'models.csv',[dict(model=m,weight=1/6,likelihood_conditioning=contract,precision_mode='strict_f64') for m in MODELS])
            (directory/'ancestral_arrays.npz').write_bytes(b'unit-fixture-presence-only')
            for i,fit in enumerate(fits,1):self.dump(directory/f'initial_fit_{i}.json',fit)
            profiles=[dict(model=m,status='numerically_uncertified' if i==1 else 'complete') for i,m in enumerate(MODELS,1)]
            self.dump(directory/'profile_status.json',profiles)
            for i,fit in enumerate(fits,1):
                if i==1:continue
                self.csv(directory/f'profile_{i}.csv',[dict(param=p) for p in fit['parameter_bounds']])
                for p in fit['parameter_bounds']:self.csv(directory/f'profile_{i}_{p}_curve.csv',[dict(x=0)])
            boot=[]
            def definition(name):
                return dict(model=dict(name=name,params={k:dict(type='free',min=v[0],max=v[1]) for k,v in self.manifest['model_parameter_bounds'][name].items()}),
                            config=dict(likelihood_conditioning=contract,precision='strict_f64'))
            for pair in range(1,4):
                null,alt=MODELS[(pair-1)*2:pair*2]
                header=dict(n_boot=2,n_alt_sim=0,seed=pair,simulation_conditioning='tip_survival',
                    null=definition(null),alt=definition(alt))
                rows=[dict(id=i,rep_id=i,phase='null',status='converged' if i==1 else 'failed') for i in (1,2)]
                self.store(directory/f'bootstrap_{pair}.jsonl',header,rows);self.csv(directory/f'bootstrap_{pair}.csv',rows)
                boot.append(dict(null=null,alt=alt,n_failed=1))
            trees=[];table=[]
            for i in (1,2):
                models=[dict(model=m,model_id=j,tree_id=i,likelihood_conditioning=contract,
                    fit_status='failed' if i==2 and j==1 else 'converged') for j,m in enumerate(MODELS,1)]
                trees.append(dict(id=i,models=models,tree_status='complete' if i==1 else 'incomplete_model_set'))
                table.extend(models)
            self.store(directory/'multitree.jsonl',dict(bsm_maps_per_tree=2,seed=1,trees=['a','b'],max_trees=2,
                match_mode='monophyletic',models=[definition(m) for m in MODELS]),trees)
            self.csv(directory/'per_tree_models.csv',table)
            for name in ('multitree_probability_bounds.csv','tree_parameter_sensitivity.csv','tree_model_winner_frequency.csv'):
                self.csv(directory/name,[dict(status='unit_fixture')])
            self.report['contracts'][contract]=dict(profiles=profiles,bootstrap=boot,single_tree_mapping=self.maps(directory,''),
                multitree=dict(n_planned=2,n_complete=1,n_incomplete=1,joint_mapping=self.maps(directory,'joint_')))

    def check(self):return verify_case_content(self.path,self.manifest,self.report)

    def test_processed_profiles_and_tree_failures_remain_explicit(self):
        result=self.check()
        self.assertEqual(result['status'],'incomplete_profiles')
        for cell in result['contracts'].values():
            self.assertEqual(cell['profiles_uncertified'],1)
            self.assertEqual(cell['trees_incomplete'],1)
            self.assertEqual(cell['bootstrap_failed_counts'],[1,1,1])

    def test_bootstrap_failure_cannot_be_erased_from_summary(self):
        self.report['contracts']['tip_survival']['bootstrap'][0]['n_failed']=0
        with self.assertRaisesRegex(ValueError,'erased'):self.check()

    def test_dropped_bootstrap_replicate_is_not_a_processed_case(self):
        path=self.path/'tip_survival/bootstrap_1.jsonl';lines=path.read_text().splitlines()
        path.write_text('\n'.join(lines[:2])+'\n')
        with self.assertRaisesRegex(ValueError,'every planned'):self.check()

    def test_cladogenesis_row_count_does_not_hide_duplicate_node(self):
        path=self.path/'tip_survival/cladogenetic_events.csv'
        with path.open() as handle:rows=list(csv.DictReader(handle))
        rows[1]['node']=rows[0]['node'];self.csv(path,rows)
        with self.assertRaisesRegex(ValueError,'duplicated internal node'):self.check()

    def test_joint_metadata_cannot_duplicate_map_id(self):
        path=self.path/'tip_survival/joint_map_metadata.csv'
        with path.open() as handle:rows=list(csv.DictReader(handle))
        rows[1]['map']=rows[0]['map'];self.csv(path,rows)
        with self.assertRaisesRegex(ValueError,'map ID'):self.check()

    def test_complete_status_cannot_hide_retained_failures(self):
        self.report['status']='complete'
        with self.assertRaisesRegex(ValueError,'terminal status'):self.check()

    def test_frozen_reference_bounds_must_match_bootstrap_and_posterior_models(self):
        for filename,key in (('bootstrap_1.jsonl','null'),('multitree.jsonl','models')):
            with self.subTest(filename=filename):
                path=self.path/'tip_survival'/filename;original=path.read_bytes()
                lines=[json.loads(s) for s in path.read_text().splitlines()]
                manifest=decode(lines[0]['manifest'])
                model=manifest[key][0] if key=='models' else manifest[key]
                model['model']['params']['d']['min']=1e-12
                self.store(path,manifest,[line['record'] for line in lines[1:]])
                with self.assertRaisesRegex(ValueError,'physical parameter domain'):self.check()
                path.write_bytes(original)


if __name__=='__main__':unittest.main()
