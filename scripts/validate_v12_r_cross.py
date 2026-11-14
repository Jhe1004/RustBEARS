"""Write deterministic Python references for the native R summary tests."""
import argparse
from pathlib import Path
import pandas as pd
import rustbears as b
from rustbears.uncertainty.averaging import model_average, node_clades, topology_hash
from rustbears.uncertainty.diagnostics import bsm_diagnostics
from rustbears.uncertainty.multitree import fit_over_trees

TREE = "((T1:1,T2:1):1,(T3:1,T4:1):1);"
GEOG = {"T1":"10", "T2":"10", "T3":"01", "T4":"01"}


def write_references(output):
    output.mkdir(parents=True, exist_ok=True)
    records=[]
    for family in ("DEC","DEC+J","DIVALIKE","DIVALIKE+J","BAYAREALIKE","BAYAREALIKE+J"):
        model = b.model_template(family).with_parameters(**{k:dict(type="fixed",init=v)
            for k,v in dict(d=.2,e=.1,j=.1 if "+J" in family else 0.).items()})
        fit=b.BioGeoBEARS(TREE,GEOG,model,likelihood_conditioning="tip_survival").fit()
        for node,clade in enumerate(node_clades(fit.model.tree)):
            for state,probability in enumerate(fit.evaluation.ancestral_probs_branch_top[node]):
                records.append(dict(model=family,clade="|".join(sorted(clade)),state=state,
                                    probability=probability,lnL=fit.lnL))
    pd.DataFrame(records).to_csv(output/"six_models.tsv",sep="\t",index=False)
    fits=[]
    for rate in (.2,.3):
        model=b.model_template("DEC").with_parameters(d=dict(type="fixed",init=rate),e=dict(type="fixed",init=.1))
        fits.append(b.BioGeoBEARS(TREE,GEOG,model,likelihood_conditioning="tip_survival").fit())
    average=model_average(fits,criterion="AIC")
    maps=average.run_bsm(num_maps=40,seed=12345)
    maps.clado_events_df.to_csv(output/"average_clado.tsv",sep="\t",index=False)
    maps.ana_events_df.to_csv(output/"average_ana.tsv",sep="\t",index=False)
    maps.map_metadata.to_csv(output/"average_metadata.tsv",sep="\t",index=False)
    bsm_diagnostics(maps,min_mean_count=0).event_summary_df.to_csv(output/"average_diagnostics.tsv",sep="\t",index=False)
    rotated="((T4:1,T3:1):1,(T2:1,T1:1):1);"
    multi=fit_over_trees([TREE,rotated]*5,GEOG,[fits[0].model.model],
        likelihood_conditioning="tip_survival",criterion="AIC",init_strategy="full_cold_start",
        cold_start_audit_fraction=0,bsm_maps_per_tree=8,num_threads=1,seed=12345)
    multi.joint_bsm.clado_events_df.to_csv(output/"multitree_clado.tsv",sep="\t",index=False)
    multi.joint_bsm.ana_events_df.to_csv(output/"multitree_ana.tsv",sep="\t",index=False)
    multi.joint_bsm.map_metadata.to_csv(output/"multitree_metadata.tsv",sep="\t",index=False)
    for target,dependence in (("fixed_sample",None),("posterior_population","independent")):
        bsm_diagnostics(multi.joint_bsm,min_mean_count=0,tree_target=target,tree_dependence=dependence).event_summary_df.to_csv(
            output/f"multitree_{target}_diagnostics.tsv",sep="\t",index=False)
    trees=[TREE,"((甲:1,乙:1):1,Ω:2);"]
    pd.DataFrame([dict(tree=tree,sha256=topology_hash(b.read_newick(tree))) for tree in trees]).to_csv(output/"topology_hashes.tsv",sep="\t",index=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,required=True)
    write_references(parser.parse_args().output)
