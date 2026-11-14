"""Analytic-reference and genuine molecular-simulation interface checks."""
import unittest
from unittest.mock import patch
from types import SimpleNamespace
from pathlib import Path
import tempfile
import json

import numpy as np
import rustbears as b
from scripts.validate_v12_v3_external_posterior import coalescent_tree,jc69_sequences,newick
from scripts.validate_v12_v4_protocol import exact_moments
from scripts.validate_v12_v4_precision import make_fit,exact_event_counts
from scripts.run_uncertainty_case_study import case_model_template
from rustbears.uncertainty.bootstrap import _nested_model
from scripts.run_v12_revbayes_posterior import diagnostics
from scripts.validate_v12_v3_external_posterior import collect


class ValidationProtocolTests(unittest.TestCase):
    def test_stuck_continuous_chain_cannot_pass_using_other_good_traces(self):
        fake=SimpleNamespace(__version__='test',rhat=lambda *a,**k:1.,ess=lambda *a,**k:3000.)
        traces={'Likelihood':np.random.default_rng(1).normal(size=(2,1000)),'population':np.ones((2,1000))}
        with patch.dict('sys.modules',{'arviz':fake}):result=diagnostics(traces,{})
        self.assertFalse(result['diagnostic_pass'])
        self.assertEqual(result['constant_numeric_traces'],['population'])
        with patch.dict('sys.modules',{'arviz':fake}):
            result=diagnostics({'Likelihood':traces['Likelihood']},{'a,b':np.ones((2,1000))})
        self.assertTrue(result['diagnostic_pass'])
        self.assertEqual(result['constant_traces'][0]['kind'],'clade_indicator')

    def test_posterior_from_another_alignment_is_rejected_before_import(self):
        with tempfile.TemporaryDirectory() as directory:
            folder=Path(directory);(folder/'alignment.nex').write_text('current alignment')
            (folder/'external_run.json').write_text(json.dumps({'status':'completed','sequence_sha256':'another-input'}))
            with patch('scripts.validate_v12_v3_external_posterior._load_dataset',return_value=(folder,{'clock_rate':.01})):
                with self.assertRaisesRegex(ValueError,'this dataset alignment'):
                    collect(SimpleNamespace(dataset=folder))

    def test_empirical_pairs_keep_common_physical_zero_domains(self):
        tree='((a:1,b:1):1,(c:1,d:1):1);'
        geography=dict(a='10',b='01',c='11',d='10')
        for family in ('DEC','DIVALIKE','BAYAREALIKE'):
            null=b.BioGeoBEARS(tree,geography,case_model_template(family),likelihood_conditioning='tip_survival')
            alt=b.BioGeoBEARS(tree,geography,case_model_template(family+'+J'),likelihood_conditioning='tip_survival')
            template,extra=_nested_model(null,alt)
            self.assertEqual(extra,('j',))
            for key in ('d','e'):
                self.assertEqual(null.model.params[key].min,0.)
                self.assertEqual(template.params[key].min,0.)

    def test_first_second_moment_recursion_matches_branch_insertions(self):
        fit=make_fit();moments=exact_moments(fit.model,fit.params)
        means=exact_event_counts(fit.model,fit.params)
        for kind in "de":
            self.assertAlmostEqual(moments[kind]["mean"],means[kind],places=12)
            self.assertAlmostEqual(sum(moments[kind]["branch_means"].values()),means[kind],places=12)
            self.assertGreater(moments[kind]["variance"],0)

    def test_one_area_surviving_tree_has_zero_losses_and_fixed_splits(self):
        template=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=0),e=dict(type="fixed",init=.2))
        model=b.BioGeoBEARS("((a:1,b:1):1,(c:1,d:1):1);",{x:"1" for x in "abcd"},template,
                            likelihood_conditioning="tip_survival")
        moments=exact_moments(model,model.model.resolve())
        self.assertAlmostEqual(moments["e"]["mean"],0,places=13)
        self.assertAlmostEqual(moments["e"]["variance"],0,places=13)
        self.assertAlmostEqual(moments["y"]["mean"],3,places=13)
        self.assertAlmostEqual(moments["y"]["variance"],0,places=12)

    def test_coalescent_is_a_true_ultrametric_tree_and_jc69_has_correct_distance(self):
        tree=coalescent_tree(50,np.random.default_rng(11))
        self.assertEqual(tree.tip_count,50)
        np.testing.assert_allclose(tree.ages[:50],0,atol=1e-12)
        self.assertAlmostEqual(tree.ages[tree.root],10)
        restored=b.read_newick(newick(tree))
        self.assertEqual(set(restored.tip_labels),set(tree.tip_labels))
        pair=b.read_newick("(a:1,b:1);")
        seq=jc69_sequences(pair,100000,.1,np.random.default_rng(12))
        observed=np.mean(np.array(list(seq["a"]))!=np.array(list(seq["b"])))
        truth=.75*(1-np.exp(-4*.2/3))
        self.assertLess(abs(observed-truth),6*np.sqrt(truth*(1-truth)/100000))


if __name__=="__main__":unittest.main()
