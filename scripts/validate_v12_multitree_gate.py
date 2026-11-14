#!/usr/bin/env python3
"""Deterministic small Gate 3 experiments on actual fitted probability models."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from unittest.mock import patch
import warnings

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import rustbears as b
from rustbears.uncertainty.averaging import model_average,node_clades
from rustbears.uncertainty.multitree import fit_over_trees
from rustbears.uncertainty.checkpoint import model_manifest,task_fingerprint,canonical_manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True)
    args=p.parse_args();started=time.monotonic()
    original="((T1:1,T2:1):1,T3:2);";rotated="(T3:2,(T2:1,T1:1):1);"
    alternate="((T1:1,T3:1):1,T2:2);"
    geo={"T1":"10","T2":"10","T3":"01"}
    definitions=[b.model_template(name).with_parameters(d=dict(type="fixed",init=.2),
        e=dict(type="fixed",init=.1),j=dict(type="fixed",init=j)) for name,j in (("DEC",0),("DEC+J",.1))]
    reference=[b.BioGeoBEARS(original,geo,m,likelihood_conditioning="tip_survival").fit() for m in definitions]
    averaged=model_average(reference)
    manifest=dict(experiment="Gate3-small-v1",models=[model_manifest(f) for f in reference],seed=20261115,
                  repeated_trees=50,nonmonophyletic_trees=10,injected_failures=2)
    common=dict(init_strategy="full_cold_start",likelihood_conditioning="tip_survival",num_threads=4,seed=20261115)
    repeated=fit_over_trees([original,rotated]*25,geo,definitions,bsm_maps_per_tree=8,**common)
    difference=float(np.max(np.abs(repeated.complete_case_probs_top-averaged.ancestral_probs_branch_top)))
    assert difference<1e-14
    required={"tree_id","tree_topology_hash","model_id","map_id","stratum_weight","precision_mode","likelihood_conditioning"}
    e4={}
    for name,table in (("per_tree",repeated.per_tree_df),("clado",repeated.joint_bsm.clado_events_df),
                       ("ana",repeated.joint_bsm.ana_events_df),("map_metadata",repeated.joint_bsm.map_metadata)):
        assert required.issubset(table.columns),(name,required-set(table.columns))
        if name!="per_tree":
            assert not table[list(required)].isna().any().any()
            ordered=table.sort_values(["tree_id","model_id","map_id"],kind="stable")
            assert np.array_equal(ordered.index,table.index)
        e4[name]=dict(rows=len(table),all_required_columns_present=True,stable_order=True)
    tree_sample=[]
    for i in range(20):
        tree=b.read_newick(original if i<10 else alternate)
        scale=1+.05*i;tree.lengths*=scale;tree.ages*=scale
        tree_sample.append(tree)
    mixed=fit_over_trees(tree_sample,geo,definitions,summary_tree=original,**common)
    node=node_clades(mixed.summary_tree).index(frozenset(("T1","T2")))
    assert abs(mixed.clade_support_all_trees[node]-.5)<1e-15
    assert mixed.clade_complete_trees_count[node]==10
    real=b.BioGeoBEARS.fit
    failed_ages=[float(tree_sample[i].ages[tree_sample[i].root]) for i in (2,6)]
    def failing(model,*a,**kw):
        if model.model.name=="DEC+J" and any(abs(model.tree.ages[model.tree.root]-x)<1e-13 for x in failed_ages):
            raise RuntimeError("predeclared test-side single-candidate optimizer failure")
        return real(model,*a,**kw)
    with patch.object(b.BioGeoBEARS,"fit",failing),warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        mono=fit_over_trees(tree_sample,geo,definitions,summary_tree=original,on_excess_failure="conservative_bounds",**common)
        mrca=fit_over_trees(tree_sample,geo,definitions,summary_tree=original,match_mode="mrca",on_excess_failure="conservative_bounds",**common)
    assert mono.n_trees_incomplete==2 and mrca.n_trees_incomplete==2
    mono_width=mono.prob_failure_bounds_top[1][node]-mono.prob_failure_bounds_top[0][node]
    mrca_width=mrca.prob_failure_bounds_top[1][node]-mrca.prob_failure_bounds_top[0][node]
    assert np.max(np.abs(mono_width-.2))<1e-15
    assert np.max(np.abs(mrca_width-.1))<1e-15
    bad=mono.per_tree_df[mono.per_tree_df.tree_status=="incomplete_model_set"]
    assert len(bad)==4 and bad.weight.isna().all() and bad.stratum_weight.isna().all()
    result=dict(schema="rustbears.multitree.gate-small.v1",status="passed",seconds=time.monotonic()-started,
        task_fingerprint=task_fingerprint(manifest),manifest=canonical_manifest(manifest),
        rotated_replicate_trees=dict(count=50,max_probability_error=difference,threshold=1e-14),
        nonmonophyletic_clade=dict(total_trees=20,nonmonophyletic=10,observed_support=float(mixed.clade_support_all_trees[node]),complete_matching_trees=int(mixed.clade_complete_trees_count[node])),
        injected_model_failures=dict(failed_tree_ids=[3,7],incomplete_model_sets=2,rows_isolated=len(bad),
            monophyletic_bound_width=mono_width.tolist(),mrca_bound_width=mrca_width.tolist(),
            threshold=1e-15,scope="test-side optimizer exception; every candidate attempted and whole failed tree excluded"),
        e4_provenance=e4,
        limitation="Small deterministic Gate3 checks; does not assert three-platform CI or real posterior-tree optimizer coverage.",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    destination=Path(args.output);destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":result["status"],"seconds":result["seconds"],"max_probability_error":difference,"output":str(destination)}))


if __name__=="__main__":
    main()
