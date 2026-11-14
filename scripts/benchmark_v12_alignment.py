#!/usr/bin/env python3
"""Audited same-host scalar alignment and timing, with explicit adapter status.

Only RustBEARS and the preserved original BioGeoBEARS adapter are executable.
Lagrange-ng/RevBayes remain pending until their actual numerical contracts are
implemented and tested. There is no generic-command fallback posing as an
adapter. The default small calibration never reports a publication speedup.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent
sys.path.insert(0,str(ROOT))


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+"\n")


def capabilities():
    return {
        "rustbears":dict(adapter="implemented",version_source="installed Python package + native SHA256",
            contracts=["legacy_unconditioned","tip_survival","unconditional_complete"],
            verified_benchmark_scope="single-layer DEC, Contract-C, legacy_f32, 2 areas, 3 tips"),
        "original_biogeobears":dict(adapter="implemented_readonly_oracle",contracts=["legacy_unconditioned"],
            verified_benchmark_scope="same scalar case after explicit Q/C/root comparison; no Contract-A/B equivalence claim"),
        "lagrange_ng":dict(adapter="not_implemented",executable_in_path=shutil.which("lagrange-ng"),
            documented="DEC; evaluate/optimize; dispersion/extinction; worker controls; Newick+PHYLIP+control file",
            root_measure="unknown_pending_source_and_runtime_check",empty_clado="unknown_pending",
            tip_survival_conditioning="unknown_pending",benchmark_eligible=False,
            primary_source="https://github.com/computations/lagrange-ng"),
        "revbayes":dict(adapter="not_implemented",executable_in_path=shutil.which("rb"),
            documented="dnPhyloCTMCClado accepts tree,Q,cladoProbs,rootFrequencies(Simplex)",
            root_measure="explicit normalized simplex documented; mapping to raw legacy weights pending",
            empty_clado="unknown_pending",tip_survival_conditioning="no explicit survival option established by this documentation; a verified correction is pending",
            benchmark_eligible=False,primary_source="https://revbayes.github.io/documentation/dnPhyloCTMCClado.html")}


def run_measured(command,folder,name):
    stdout,stderr=folder/f"{name}.stdout.log",folder/f"{name}.stderr.log"
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1",
                                   PYTHONDONTWRITEBYTECODE="1",MPLCONFIGDIR="/tmp/rustbears-matplotlib")
    started=time.perf_counter()
    with stdout.open("w") as out,stderr.open("w") as err:
        process=subprocess.Popen(command,stdout=out,stderr=err,cwd=ROOT,env=env)
        if hasattr(os,"wait4"):
            _,status,usage=os.wait4(process.pid,0)
            process.returncode=os.waitstatus_to_exitcode(status)
            rss=int(usage.ru_maxrss*(1 if sys.platform=="darwin" else 1024))
        else:
            process.wait();rss=None
    record=dict(command=command,returncode=process.returncode,wall_seconds=time.perf_counter()-started,
        peak_rss_bytes=rss,peak_rss_scope="whole child process including startup/imports" if rss is not None else "unavailable_on_this_host",
        stdout=str(stdout),stderr=str(stderr),stdout_sha256=sha(stdout),stderr_sha256=sha(stderr),
        executable_sha256=sha(command[0]),thread_environment={k:env[k] for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS")})
    if process.returncode:raise RuntimeError(f"adapter failed ({process.returncode}); see {stderr}")
    return record


def rust_worker(spec_path,output):
    import rustbears as b
    from rustbears._ffi import _resolve_native_library_path
    from scripts.validate_v12_gate1 import ReferenceTree
    spec=json.loads(Path(spec_path).read_text())
    model=b.BioGeoBEARS(spec["tree_file"],spec["geog_file"],"DEC",max_range_size=spec["max_range_size"],precision="legacy_f32")
    params=model.model.resolve(spec["params"])
    from rustbears.bsm import _layers
    layer=list(_layers(model,params));q=layer[0];events=layer[2].reshape((-1,4));n=len(model.states)
    c=np.zeros((n,n,n))
    for a,l,r,w in events:c[int(a),int(l),int(r)]+=w
    for i in range(1,n):c[i]/=c[i].sum()
    model.evaluate(params,calc_ancprobs=False)
    values=[];optimized=None;started=time.perf_counter()
    for _ in range(spec["repeats"]):
        if spec["phase"]=="evaluate":value=model.evaluate(params,calc_ancprobs=False).lnL
        else:
            optimized=model.fit(initial=params,seed=spec["seed"],maxiter=spec["maxiter"])
            value=optimized.lnL
        values.append(value)
    elapsed=time.perf_counter()-started
    result=dict(adapter="rustbears",version=b.__version__,python=sys.version,native_path=str(_resolve_native_library_path()),
        native_sha256=sha(_resolve_native_library_path()),phase=spec["phase"],repeats=spec["repeats"],
        measured_api_seconds=elapsed,lnL=values,state_labels=[",".join(map(str,s)) or "_" for s in model.states],
        Q_colmajor=q.ravel(order="F").tolist(),C_colmajor=c.ravel(order="F").tolist(),n_states=n,
        root_weights=model.root_likelihood_weights.tolist(),contract=model.likelihood_conditioning,empty_clado="zero_row",
        optimized_params=optimized.params if optimized is not None else None,
        optimization_protocol="RustBEARS default multi-start/refinement L-BFGS-B; unlike original optimizer work budget" if optimized else None)
    write_json(output,result)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir")
    p.add_argument("--phase",choices=("evaluate","optimization"),default="evaluate")
    p.add_argument("--repeats",type=int,default=3)
    p.add_argument("--worker",choices=("rustbears",))
    p.add_argument("--spec");p.add_argument("--worker-output")
    p.add_argument("--inventory-only",action="store_true")
    args=p.parse_args()
    if args.worker:
        rust_worker(args.spec,args.worker_output);return
    if not args.output_dir or args.repeats<1:p.error("output-dir and positive repeats are required")
    out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    info=capabilities();write_json(out/"capability_inventory.json",info)
    if args.inventory_only:
        print(json.dumps(dict(status="inventory_written",output=str(out))));return
    m0=next(PROJECT.glob("04_M0_黄金参考集*"))
    rscript=m0/"runtime/r44-final/bin/Rscript";oracle=m0/"runtime/oracle-library"
    if not rscript.is_file() or not (oracle/"BioGeoBEARS").is_dir():
        raise FileNotFoundError("preserved original BioGeoBEARS runtime is unavailable; no adapter substitution is allowed")
    tree=out/"tree.newick";geog=out/"geography.data"
    tree.write_text("((T1:1,T2:1):1,T3:2);\n");geog.write_text("3 2 (A B)\nT1 10\nT2 10\nT3 01\n")
    cases=[]
    for index,(d,e) in enumerate(((.2,.1),(.05,.02),(.4,.2)),1):
        folder=out/f"case_{index}";folder.mkdir(exist_ok=True)
        spec=dict(tree_file=str(tree),geog_file=str(geog),oracle_library=str(oracle),max_range_size=2,
            workdir=str(folder),params=dict(d=d,e=e,j=0),phase=args.phase,repeats=args.repeats,seed=20261115,maxiter=500,
            root_measure="unnormalized unit weights on nonempty states",empty_clado="zero_row",conditioning="legacy_unconditioned",
            tree_sha256=sha(tree),geography_sha256=sha(geog),warmup_calls=1,
            timer_scope="fixed-parameter public likelihood API; startup reported separately" if args.phase=="evaluate" else "full default optimization; algorithms/budgets differ and no optimizer speedup is issued")
        spec_path=folder/"spec.json";write_json(spec_path,spec)
        rust_output=folder/"rustbears.json";oracle_output=folder/"original_biogeobears.json"
        commands={"rustbears":[sys.executable,str(Path(__file__).resolve()),"--worker","rustbears","--spec",str(spec_path),"--worker-output",str(rust_output)],
            "original_biogeobears":[str(rscript),"--vanilla",str(ROOT/"scripts/benchmark_v12_oracle.R"),str(spec_path),str(oracle_output)]}
        timing={name:run_measured(command,folder,name) for name,command in commands.items()}
        rust=json.loads(rust_output.read_text());original=json.loads(oracle_output.read_text())
        native_paths=original["native_paths"]
        for path in (native_paths.values() if isinstance(native_paths,dict) else native_paths):
            if Path(path).is_file():timing["original_biogeobears"].setdefault("native_sha256",{})[path]=sha(path)
        nq=np.max(np.abs(np.asarray(rust["Q_colmajor"])-original["Q_colmajor"]))
        nc=np.max(np.abs(np.asarray(rust["C_colmajor"])-original["C_colmajor"]))
        same=(rust["state_labels"]==original["state_labels"] and rust["root_weights"]==original["root_weights"] and
              rust["contract"]==original["contract"] and rust["empty_clado"]==original["empty_clado"] and nq<=1e-12 and nc<=1e-12)
        error=float(np.max(np.abs(np.atleast_1d(rust["lnL"])-np.atleast_1d(original["lnL"]))))
        cases.append(dict(case=index,spec_sha256=sha(spec_path),mathematical_alignment_verified=bool(same),
            maximum_Q_error=float(nq),maximum_C_error=float(nc),maximum_lnL_error=error,
            scalar_tolerance=1e-8,scalar_alignment_passed=bool(same and error<=1e-8),
            execution=timing,api_seconds={"rustbears":rust["measured_api_seconds"],"original_biogeobears":original["measured_api_seconds"]},
            speedup=None,speedup_status="not_reported_small_alignment_calibration" if args.phase=="evaluate" else "not_reported_unmatched_optimizer_work_protocol"))
        write_json(out/"progress.json",cases)
    result=dict(status="scalar_alignment_complete" if all(c["scalar_alignment_passed"] for c in cases) else "alignment_failed_no_speedup",
        host=dict(platform=platform.platform(),machine=platform.machine(),node=platform.node()),phase=args.phase,
        adapters=info,cases=cases,full_benchmark_completed=False,
        protocol="No speedup until exact tree/geography/state order/Q/C/root measure/empty-state propagation/conditioning, numerical tolerance and timing scope are verified. Evaluation and optimization are separate experiments.",
        source_sha256={name:sha(ROOT/name) for name in ("scripts/benchmark_v12_alignment.py","scripts/benchmark_v12_oracle.R")})
    write_json(out/"results.json",result)
    print(json.dumps(dict(status=result["status"],max_lnL_error=max(c["maximum_lnL_error"] for c in cases),output=str(out))))


if __name__=="__main__":main()
