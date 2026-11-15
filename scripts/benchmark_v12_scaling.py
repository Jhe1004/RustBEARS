#!/usr/bin/env python3
"""Same-host, fresh-process scalar workloads with mandatory numerical alignment.

Four software packages are represented, but their default models are not
assumed equivalent. Standard DEC: original BioGeoBEARS / Rust legacy precision,
and RevBayes / Rust strict precision with a normalized root. Lagrange-ng is a
separate two-area v=2 stratum. Every process constructs one model and evaluates
one fixed likelihood. Thus time includes interpreter startup and input parsing;
it is not a kernel or optimizer speedup. Matrix exports are an untimed preflight.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CASES = {4: (2, 2), 16: (4, 4), 57: (6, 4), 128: (7, 7)}
THREAD_ENV = dict(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def make_inputs(folder, tips, areas):
    from scripts.run_uncertainty_validation import balanced_tree
    (folder / "tree.newick").write_text(balanced_tree(tips) + "\n")
    # Fixed cycling singleton observations, selected without fit or runtime results.
    ranges = {f"T{i+1}": "".join(str(int(j == i % areas)) for j in range(areas)) for i in range(tips)}
    rows = "\n".join(f"{key} {value}" for key, value in ranges.items()) + "\n"
    names = " ".join(chr(65+i) for i in range(areas))
    (folder / "geography.data").write_text(f"{tips} {areas} ({names})\n" + rows)
    (folder / "geography.phy").write_text(f"{tips} {areas}\n" + rows)
    (folder / "geography.nex").write_text(f'#NEXUS\nbegin data;\ndimensions ntax={tips} nchar={areas};\n'
        'format datatype=standard symbols="01";\nmatrix\n' + rows + ';\nend;\n')


def rust_worker(spec_path, output, variant):
    import rustbears as b
    from rustbears._ffi import _resolve_native_library_path
    from rustbears.bsm import _layers
    from rustbears.uncertainty.contracts import RootPriorSpec
    spec = json.loads(Path(spec_path).read_text())
    definition = b.model_template("DEC")
    if variant == "lagrange":
        definition = definition.with_parameters(v=dict(type="fixed", init=2))
    options = dict(max_range_size=spec["max_range_size"],
        likelihood_conditioning="legacy_unconditioned", precision="legacy_f32" if variant == "original" else "strict_f64")
    if variant == "revbayes":
        options["root_prior"] = RootPriorSpec("uniform_nonnull")
    started = time.perf_counter()
    model = b.BioGeoBEARS(spec["tree_file"], spec["geog_file"], definition, **options)
    params = model.model.resolve(spec["params"])
    value = model.evaluate(params, calc_ancprobs=False, num_threads=1).lnL
    elapsed = time.perf_counter() - started
    result = dict(lnL=value, version=b.__version__, construction_evaluation_seconds=elapsed,
        native_sha256=sha(_resolve_native_library_path()), precision=options["precision"],
        root_weights=model.root_likelihood_weights.tolist(),
        state_labels=[",".join(map(str, state)) or "_" for state in model.states])
    if spec["export_components"]:
        layer = list(_layers(model, params)); n = len(model.states)
        c = np.zeros((n, n, n))
        for a, left, right, weight in layer[2].reshape(-1, 4):
            c[int(a), int(left), int(right)] += weight
        c[1:] /= c[1:].sum(axis=(1, 2))[:, None, None]
        np.savez(Path(output).with_suffix(".npz"), Q=layer[0], C=c)
    write_json(output, result)


def rev_script(spec):
    areas = spec["areas"]; n = spec["states"]; maximum = spec["max_range_size"]
    d, e = spec["params"]["d"], spec["params"]["e"]
    dispersal = json.dumps([[0 if i == j else d for j in range(areas)] for i in range(areas)])
    extinction = json.dumps([[e if i == j else 0 for j in range(areas)] for i in range(areas)])
    source = f'''tree <- readTrees("tree.newick")[1]
data01 = readDiscreteCharacterData("geography.nex")
data = formatDiscreteCharacterData(data01,"DEC",numStates={n})
Q := fnDECRateMatrix(dispersalRates={dispersal},extirpationRates={extinction},nullRange="Include",maxRangeSize={maximum})
ep ~ dnDirichlet([1.0,1.0])
ep.clamp(simplex(1,1))
C := fnDECCladoProbs(eventProbs=ep,eventTypes=["s","a"],numCharacters={areas},maxRangeSize={maximum},probType="pattern",widespreadAllopatry=FALSE)
x ~ dnPhyloCTMCClado(tree=tree,Q=Q,cladoProbs=C,rootFrequencies=simplex({','.join(['0']+['1']*(n-1))}),nSites=1,type="NaturalNumbers")
x.clamp(data)
write(x.lnProbability(),file="revbayes.lnL")
'''
    if spec["export_components"]:
        # mnFile's JSON serialization exports full double precision; print(C)
        # hard-codes four decimals and cannot establish matrix equality.
        source += f'''write(data.getStateDescriptions(),file="revbayes.states.txt")
for(i in 1:{n}){{for(j in 1:{n}){{write(i,j,Q[i][j],"\\n",file="revbayes.Q.tsv",append=TRUE)}}}}
monitors=VectorMonitors()
monitors.append(mnFile(C,filename="revbayes.C.jsonl",printgen=1,format="json"))
moves=VectorMoves()
mdl=model(C)
chain=mcmc(mdl,moves,monitors)
chain.run(generations=1)
'''
    return source + "q()\n"


def measured(command, folder, name, timeout):
    """POSIX wait4 measures each process separately, with bounded timeout cleanup."""
    if not hasattr(os, "wait4"):
        raise RuntimeError("This Mac/Linux benchmark requires per-process wait4 RSS")
    cache = folder / "runtime_cache"; cache.mkdir(exist_ok=True)
    env = dict(os.environ, **THREAD_ENV, MPLCONFIGDIR=str(cache), TMPDIR=str(cache),
               PYTHONDONTWRITEBYTECODE="1")
    if name == "lagrange_ng":
        env["DYLD_LIBRARY_PATH"] = str(Path(command[0]).parent / "lib")
    out, err = folder / f"{name}.stdout.log", folder / f"{name}.stderr.log"
    load_before = list(os.getloadavg())
    started = time.perf_counter(); timed_out = False
    with out.open("w") as stdout, err.open("w") as stderr:
        process = subprocess.Popen(command, cwd=folder, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
        while True:
            pid, status, usage = os.wait4(process.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() - started > timeout:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                _, status, usage = os.wait4(process.pid, 0)
                break
            time.sleep(.002)
        process.returncode = os.waitstatus_to_exitcode(status)
    return dict(returncode=process.returncode, timed_out=timed_out, timeout_seconds=timeout,
        wall_seconds=time.perf_counter()-started,
        peak_rss_bytes=int(usage.ru_maxrss)*(1 if sys.platform == "darwin" else 1024),
        cpu_seconds=usage.ru_utime+usage.ru_stime, command=list(map(str, command)),
        executable_sha256=sha(command[0]), stdout_sha256=sha(out), stderr_sha256=sha(err),
        load_before=load_before, load_after=list(os.getloadavg()), thread_environment=THREAD_ENV,
        wall_timer_poll_interval_seconds=.002)


def commands_for(args, folder, spec):
    commands = {}
    for variant in ("original", "revbayes") + (("lagrange",) if spec["states"] == 4 else ()):
        commands[f"rust_{variant}"] = [sys.executable, str(Path(__file__).resolve()), "--worker", variant,
            "--spec", str(folder/"spec.json"), "--worker-output", str(folder/f"rust_{variant}.json")]
    commands["original_biogeobears"] = [str(args.rscript), "--vanilla", str(ROOT/"scripts/benchmark_v12_scaling_oracle.R"),
                                        str(folder/"spec.json"), str(folder/"original_biogeobears.json")]
    (folder/"benchmark.Rev").write_text(rev_script(spec))
    commands["revbayes"] = [str(args.revbayes), "-o", "outputPrecision=17", "-o", "useScaling=true", str(folder/"benchmark.Rev")]
    if spec["states"] == 4:
        d, e = spec["params"]["d"], spec["params"]["e"]
        (folder/"lagrange.conf").write_text('treefile = tree.newick\ndatafile = geography.phy\nareanames = A B\n'
            'mode = evaluate\nworkers = 1\nthreads-per-worker = 1\nexpm-mode = pade\nprefix = lagrange\nstates\n')
        (folder/"lagrange.ckp").write_text(f"iteration,default:d,default:e,final\n0,{d},{e},1\n")
        commands["lagrange_ng"] = [str(args.lagrange), "lagrange.conf"]
    return commands


def read_scores(folder, states, successful):
    scores = {}
    for name in successful:
        if name.startswith("rust_") or name == "original_biogeobears":
            scores[name] = float(json.loads((folder/f"{name}.json").read_text())["lnL"])
        elif name == "revbayes":
            scores[name] = float((folder/"revbayes.lnL").read_text())
        else:
            result = json.loads((folder/"lagrange.results.json").read_text())
            params = result["params"][0]
            if params["dispersion"] != .015 or params["extinction"] != .008:
                raise ValueError("Lagrange-ng did not apply the predeclared rates")
            values = [row["llh"]-math.log(row["ratio"]) for row in result["node-results"][0]["states"] if row["ratio"] > 0]
            if max(values)-min(values) > 1e-8:
                raise ValueError("Lagrange-ng root-state normalizers disagree")
            scores[name] = values[0]
    if not all(math.isfinite(v) for v in scores.values()):
        raise ValueError("Nonfinite scalar likelihood")
    return scores


def align_components(folder, states, successful):
    checks = {}
    if {"rust_original", "original_biogeobears"} <= set(successful):
        rust = json.loads((folder/"rust_original.json").read_text())
        oracle = json.loads((folder/"original_biogeobears.json").read_text())
        native = np.load(folder/"rust_original.npz")
        q = np.asarray(oracle["Q_colmajor"]).reshape((states, states), order="F")
        c = np.asarray(oracle["C_colmajor"]).reshape((states, states, states), order="F")
        checks["original_biogeobears"] = dict(max_Q_error=float(np.max(abs(q-native["Q"]))),
            max_C_error=float(np.max(abs(c-native["C"]))),
            labels_equal=rust["state_labels"] == oracle["state_labels"],
            root_equal=rust["root_weights"] == oracle["root_weights"])
    if {"rust_revbayes", "revbayes"} <= set(successful):
        rust = json.loads((folder/"rust_revbayes.json").read_text())
        bits = [s.strip() for s in (folder/"revbayes.states.txt").read_text().strip("[] \n").split(",")]
        # RevBayes's character converter exposes all 2^areas descriptions even
        # when numStates/maxRangeSize truncate Q/C. The actual process uses the
        # first n cardinality-ordered states; verify the excluded suffix too.
        maximum = CASES[states][1]
        if len(bits) != 2**CASES[states][0] or any(value.count("1") > maximum for value in bits[:states]) or any(value.count("1") <= maximum for value in bits[states:]):
            raise ValueError("RevBayes range-size truncation is not the declared state prefix")
        labels = [",".join(str(i) for i, char in enumerate(value) if char == "1") or "_" for value in bits[:states]]
        if len(labels) != states or set(labels) != set(rust["state_labels"]):
            raise ValueError("RevBayes state descriptions do not match the declared state space")
        permutation = [labels.index(label) for label in rust["state_labels"]]
        q = np.loadtxt(folder/"revbayes.Q.tsv")[:, 2].reshape(states, states)
        c = np.zeros((states, states, states))
        row = next(json.loads(line) for line in (folder/"revbayes.C.jsonl").read_text().splitlines() if "\"Iteration\"" in line and "\"C\"" in line)
        for a, left, right, weight in row["C"]:
            c[a, left, right] += weight
        native = np.load(folder/"rust_revbayes.npz")
        checks["revbayes"] = dict(max_Q_error=float(np.max(abs(q[np.ix_(permutation,permutation)]-native["Q"]))),
            max_C_error=float(np.max(abs(c[np.ix_(permutation,permutation,permutation)]-native["C"]))),
            state_permutation=permutation, labels_equal=True,
            root_equal=bool(np.allclose(rust["root_weights"], [0]+[1/(states-1)]*(states-1), rtol=0, atol=1e-15)))
    if states == 4 and "rust_lagrange" in successful and "lagrange_ng" in successful:
        # Pinned official Workspace split generator, independently verified by
        # the scalar adapter. Duplicate singleton-vicariance terms are retained.
        c = np.zeros((4,4,4)); c[1,1,1] = c[2,2,2] = 1
        for x in (1,2):
            for left, right in ((x,3),(3,x),(x,3^x),(3^x,x)):
                c[3,left,right] += 1/8
        native = np.load(folder/"rust_lagrange.npz")
        d,e = .015,.008
        q = np.array([[0,0,0,0],[e,-e-d,0,d],[e,0,-e-d,d],[0,e,e,-2*e]])
        checks["lagrange_ng"] = dict(max_Q_error=float(np.max(abs(q-native["Q"]))),
            max_C_error=float(np.max(abs(c-native["C"]))), labels_equal=True, root_equal=True,
            explicit_variant="two-area v=2; not standard DEC; source-defined C and runtime scalar audit",
            root_equivalence="empty-root contribution identically zero; nonempty raw unit weights")
    for check in checks.values():
        check["passed"] = bool(check["labels_equal"] and check["root_equal"] and
            check["max_Q_error"] <= 1e-12 and check["max_C_error"] <= 1e-12)
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rscript", type=Path); parser.add_argument("--oracle-library", type=Path)
    parser.add_argument("--revbayes", type=Path); parser.add_argument("--lagrange", type=Path)
    parser.add_argument("--mode", choices=("smoke","formal"), default="smoke")
    parser.add_argument("--gate-status", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--tips", default="19,50,150"); parser.add_argument("--states", default="4,16,57,128")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--worker", choices=("original","revbayes","lagrange"))
    parser.add_argument("--spec", type=Path); parser.add_argument("--worker-output", type=Path)
    args = parser.parse_args()
    if args.worker:
        return rust_worker(args.spec, args.worker_output, args.worker)
    if args.output is None:
        parser.error("--output is required")
    tips = [int(v) for v in args.tips.split(",")]; states_list = [int(v) for v in args.states.split(",")]
    if min(tips) < 3 or not set(states_list) <= CASES.keys() or args.repeats < 1 or args.timeout <= 0:
        parser.error("invalid dimensions, repetitions or timeout")
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    if (output/"manifest.json").exists():
        raise FileExistsError("use a new output directory; benchmark runs are immutable")
    source = {name: sha(ROOT/name) for name in ("scripts/benchmark_v12_scaling.py","scripts/benchmark_v12_scaling_oracle.R")}
    source.update({str(path.relative_to(ROOT)):sha(path) for path in sorted((ROOT/"rustbears").rglob("*.py"))})
    source["scripts/run_uncertainty_validation.py"] = sha(ROOT/"scripts/run_uncertainty_validation.py")
    from rustbears._ffi import _resolve_native_library_path
    native_path = Path(_resolve_native_library_path())
    try:
        commit = subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
        dirty = bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip())
    except (OSError,subprocess.CalledProcessError):
        commit,dirty = None,None
    manifest = dict(schema="rustbears-scaling-v1", mode=args.mode, tips=tips, states=states_list,
        repeats=args.repeats, timeout_per_process_seconds=args.timeout, params=dict(d=.015,e=.008,j=0),
        timing_scope="fresh process: interpreter startup, imports, model construction, one fixed-parameter likelihood; no optimizer or kernel speedup",
        timed_matrix_exports=False, preflight="complete Q/C/root/state-order comparison outside timed workloads",
        order="independent predeclared permutation per cell/repetition; PCG64 seed 20261115; no best-of-time selection", source_sha256=source,
        source_commit=commit,source_dirty=dirty,native_path=str(native_path),native_sha256=sha(native_path),
        host=dict(platform=platform.platform(),machine=platform.machine(),node=platform.node(),logical_cpus=os.cpu_count()),
        thread_environment=THREAD_ENV,
        unsupported="Lagrange-ng standard DEC is not equivalent; only the explicit two-area v=2 stratum is compared")
    if args.mode == "formal" and args.execute:
        if args.gate_status is None:
            parser.error("formal execution requires --gate-status")
        gates = json.loads(args.gate_status.read_text())
        if not all(gates.get(f"gate{i}_passed") is True for i in (2,3)):
            raise ValueError("formal scaling requires completed Gate 2 and Gate 3")
        manifest["gate_status_sha256"] = sha(args.gate_status)
    write_json(output/"manifest.json", manifest)
    if not args.execute:
        return
    for key in ("rscript","revbayes","lagrange","oracle_library"):
        path = getattr(args,key)
        if path is None or not path.exists():
            raise FileNotFoundError(f"required genuine runtime missing: {key}")
        setattr(args,key,path.resolve())
    oracle_hashes = {}
    for name in ("BioGeoBEARS","cladoRcpp","rexpokit"):
        for path in sorted((args.oracle_library/name).rglob("*")):
            if path.is_file():oracle_hashes[str(path.relative_to(args.oracle_library))] = sha(path)
    write_json(output/"original_oracle_files.json",oracle_hashes)
    manifest["original_oracle_files_sha256"] = sha(output/"original_oracle_files.json")
    manifest["runtime_executables"] = {key:dict(path=str(getattr(args,key)),sha256=sha(getattr(args,key))) for key in ("rscript","revbayes","lagrange")}
    write_json(output/"manifest.json",manifest)
    records = []; alignments = {}; order_rng=np.random.default_rng(20261115)
    for states in states_list:
        for n_tips in tips:
            for replicate in range(-1,args.repeats):
                # Component preflight is per state space; score agreement is
                # required for each actual tree and repetition as well.
                if replicate == -1 and states in alignments:
                    continue
                folder = output/f"n{n_tips}_s{states}"/("preflight" if replicate == -1 else f"rep{replicate+1}")
                folder.mkdir(parents=True)
                areas, maximum = CASES[states]
                make_inputs(folder,n_tips,areas)
                spec = dict(tree_file=str(folder/"tree.newick"),geog_file=str(folder/"geography.data"),
                    workdir=str(folder), oracle_library=str(args.oracle_library), max_range_size=maximum,
                    areas=areas, states=states, params=manifest["params"], export_components=replicate == -1)
                write_json(folder/"spec.json",spec)
                commands = commands_for(args,folder,spec); names=list(commands)
                order=[names[int(i)] for i in order_rng.permutation(len(names))]
                execution = {name:measured(commands[name],folder,name,args.timeout) for name in order}
                successful = [name for name,row in execution.items() if row["returncode"] == 0]
                errors = {}
                scores = {}
                for name in successful:
                    try:
                        scores.update(read_scores(folder,states,[name]))
                    except Exception as exc:
                        errors[name] = f"{type(exc).__name__}: {exc}"
                if replicate == -1:
                    alignments[states] = {}
                    for name,variant in (("original_biogeobears","original"),("revbayes","revbayes"),("lagrange_ng","lagrange")):
                        selected = [key for key in (name,f"rust_{variant}") if key in scores]
                        try:
                            alignments[states].update(align_components(folder,states,selected))
                        except Exception as exc:
                            alignments[states][name] = dict(passed=False,error=f"{type(exc).__name__}: {exc}")
                comparisons={}
                for name,variant in (("original_biogeobears","original"),("revbayes","revbayes"),("lagrange_ng","lagrange")):
                    rust = f"rust_{variant}"
                    if name not in commands:
                        continue
                    error = abs(scores[name]-scores[rust]) if name in scores and rust in scores else None
                    aligned = bool(alignments[states].get(name,{}).get("passed",False) and error is not None and error <= 1e-8)
                    comparisons[name] = dict(rust_reference=rust, maximum_lnL_error=error, scalar_tolerance=1e-8,
                        alignment_passed=aligned, fresh_process_time_ratio=execution[name]["wall_seconds"]/execution[rust]["wall_seconds"] if aligned and replicate>=0 else None)
                records.append(dict(tips=n_tips,states=states,replicate=replicate,preflight=replicate==-1,
                    execution_order=order,
                    input_sha256={f:sha(folder/f) for f in ("tree.newick","geography.data")},
                    execution=execution,scores=scores,errors=errors,comparisons=comparisons))
                write_json(output/"progress.json",dict(records=records,alignments=alignments))
                print(json.dumps(dict(tips=n_tips,states=states,replicate=replicate,
                    aligned={k:v["alignment_passed"] for k,v in comparisons.items()},errors=errors)),flush=True)
    source_stable = all(sha(ROOT/name)==value for name,value in source.items()) and sha(native_path)==manifest["native_sha256"]
    oracle_stable = all(sha(args.oracle_library/name)==value for name,value in oracle_hashes.items())
    complete = all(c["alignment_passed"] for row in records for c in row["comparisons"].values()) and source_stable and oracle_stable
    write_json(output/"results.json",dict(status="completed" if complete else "completed_with_unsupported_or_failed_comparisons",
        formal=args.mode=="formal", all_declared_comparisons_aligned=complete, source_stable=source_stable,original_oracle_stable=oracle_stable,
        manifest=manifest,alignments=alignments,records=records))


if __name__ == "__main__":
    main()
