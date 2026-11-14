"""Concrete pinned RevBayes/Lagrange-ng scalar adapters; never a speed benchmark.

Requires the official binaries and source archives isolated under stage 19.
The Lagrange comparison uses an explicit two-area v=2 variant. Standard DEC
is deliberately reported unsupported because its C matrix differs.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT.parent / "19_RustBEARS_v1_2_开发与验证_20261115"
sys.path.insert(0, str(ROOT))
import rustbears as b
from rustbears.bsm import _layers
from rustbears.uncertainty.contracts import RootPriorSpec

TREE = "((T1:1,T2:1):1,T3:2);"
GEOG = {"T1": "10", "T2": "10", "T3": "01"}
BITS = ("00", "10", "01", "11")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def measured(command, work, name, extra_env=None):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    env.update(extra_env or {})
    out, err = work / f"{name}.stdout.log", work / f"{name}.stderr.log"
    start = time.perf_counter()
    with out.open("w") as stdout, err.open("w") as stderr:
        run = subprocess.run(list(map(str, command)), cwd=work, env=env, stdout=stdout, stderr=stderr, timeout=60)
    result = dict(command=list(map(str, command)), returncode=run.returncode,
                  process_seconds=time.perf_counter() - start, timer_scope="entire process; no speedup reported",
                  executable_sha256=sha(command[0]), stdout=str(out), stderr=str(err),
                  stdout_sha256=sha(out), stderr_sha256=sha(err))
    if run.returncode:
        raise RuntimeError(f"{name} failed ({run.returncode}); see {out} and {err}")
    return result


def nexus(path, geography):
    rows = "\n".join(f"{name} {value}" for name, value in geography.items())
    path.write_text(f'#NEXUS\nbegin data;\ndimensions ntax={len(geography)} nchar=2;\n'
                    f'format datatype=standard symbols="01";\nmatrix\n{rows}\n;\nend;\n')


def rust_model(d, e, *, vicariance=1, root="uniform_nonnull", contract="legacy_unconditioned"):
    definition = b.model_template("DEC")
    if vicariance != 1:
        definition = definition.with_parameters(v=dict(type="fixed", init=vicariance))
    # The historical C likelihood permits explicit custom clado weights. New
    # uncertainty APIs intentionally reject this nonstandard v=2 model.
    prior = None if vicariance != 1 and root == "legacy_flat_weights" else RootPriorSpec(root)
    model = b.BioGeoBEARS(TREE, GEOG, definition, precision="strict_f64", likelihood_conditioning=contract,
                         root_prior=prior)
    params = model.model.resolve(dict(d=d, e=e))
    layer = list(_layers(model, params))
    q = layer[0]
    c = np.zeros((4, 4, 4))
    for a, l, r, weight in layer[2].reshape(-1, 4):
        c[int(a), int(l), int(r)] += weight
    c[1:] /= c[1:].sum(axis=(1, 2))[:, None, None]
    return model.evaluate(params, calc_ancprobs=False).lnL, q, c, model.root_likelihood_weights


def independent_score(q, c, root):
    # Fixed three-tip tree: one-hot tips, exact SciPy exponentials, tensor contraction.
    p1, p2 = expm(q), expm(2 * q)
    left = np.einsum("ijk,j,k->i", c, p1[:, 1], p1[:, 1])
    root_like = np.einsum("ijk,j,k->i", c, p1 @ left, p2[:, 2])
    return math.log(float(root @ root_like))


def rev_base(d, e, null_range):
    return f'''tree <- readTrees("tree.newick")[1]
data01 = readDiscreteCharacterData("geography.nex")
data = formatDiscreteCharacterData(data01,"DEC")
Q := fnDECRateMatrix(dispersalRates=[[0.0,{d}],[{d},0.0]],extirpationRates=[[{e},0.0],[0.0,{e}]],nullRange="{null_range}")
ep ~ dnDirichlet([1.0,1.0])
ep.clamp(simplex(1,1))
C := fnDECCladoProbs(eventProbs=ep,eventTypes=["s","a"],numCharacters=2,probType="pattern")
'''


def revbayes_case(binary, work, d, e, audit_c=False):
    records, values = {}, {}
    for null_range in ("Include", "CondSurv"):
        name = f"revbayes_{null_range}"
        source = rev_base(d, e, null_range)
        source += f'''x ~ dnPhyloCTMCClado(tree=tree,Q=Q,cladoProbs=C,rootFrequencies=simplex(0,1,1,1),nSites=1,type="NaturalNumbers")
x.clamp(data)
write(x.lnProbability(),file="{name}.lnL")
P := Q.getTransitionProbabilities(1.0)
for(i in 1:4){{for(j in 1:4){{write(i,j,Q[i][j],P[i][j],"\\n",file="{name}.QP.tsv",append=TRUE)}}}}
write(C,file="{name}.C_display_rounded.txt")
write(data.getStateDescriptions(),file="{name}.states.txt")
'''
        script = work / f"{name}.Rev"
        script.write_text(source)
        qp = work / f"{name}.QP.tsv"
        qp.unlink(missing_ok=True)
        records[name] = measured([binary, "-o", "outputPrecision=17", "-o", "useScaling=false", script.name], work, name)
        array = np.loadtxt(qp)
        values[null_range] = dict(lnL=float((work / f"{name}.lnL").read_text()),
                                 Q=array[:, 2].reshape(4, 4), P=array[:, 3].reshape(4, 4))
    if audit_c:
        # Runtime tomography identifies all 64 C entries, including the null row.
        # A zero-length cherry with a point root distribution has L = C[a,l,r].
        (work / "cherry.newick").write_text("(L:0,R:0);\n")
        for l in range(4):
            for r in range(4):
                nexus(work / f"pair_{l}_{r}.nex", dict(L=BITS[l], R=BITS[r]))
        source = rev_base(d, e, "Include") + '''tree0 <- readTrees("cherry.newick")[1]
rf ~ dnDirichlet([1.0,1.0,1.0,1.0])
rf.clamp(simplex(0,1,0,0))
z ~ dnPhyloCTMCClado(tree=tree0,Q=Q,cladoProbs=C,rootFrequencies=rf,nSites=1,type="NaturalNumbers")
'''
        for a in range(4):
            source += "rf.clamp(simplex(" + ",".join(str(int(i == a)) for i in range(4)) + "))\n"
            for l in range(4):
                for r in range(4):
                    source += f'''data01 = readDiscreteCharacterData("pair_{l}_{r}.nex")
data = formatDiscreteCharacterData(data01,"DEC")
z.clamp(data)
write({a},{l},{r},z.lnProbability(),"\\n",file="C_runtime.tsv",append=TRUE)
'''
        (work / "C_runtime.Rev").write_text(source)
        (work / "C_runtime.tsv").unlink(missing_ok=True)
        records["revbayes_C_runtime"] = measured([binary, "-o", "outputPrecision=17", "-o", "useScaling=false", "C_runtime.Rev"], work, "revbayes_C_runtime")
        values["C"] = np.exp(np.loadtxt(work / "C_runtime.tsv")[:, 3]).reshape(4, 4, 4)
    return records, values


def lagrange_case(binary, work, d, e):
    (work / "geography.phy").write_text("3 2\nT1 10\nT2 10\nT3 01\n")
    (work / "lagrange.conf").write_text(f'''treefile = tree.newick
datafile = geography.phy
areanames = A B
mode = evaluate
dispersion = {d}
extinction = {e}
workers = 1
threads-per-worker = 1
expm-mode = pade
prefix = lagrange
states
''')
    # Official resume interface: this release overwrites direct evaluation rates
    # with Workspace::getDefaultParams(). Verify the actual output rates below.
    (work / "lagrange.ckp").write_text(f"iteration,default:d,default:e,final\n0,{d},{e},1\n")
    record = measured([binary, "lagrange.conf"], work, "lagrange",
                      dict(DYLD_LIBRARY_PATH=str(binary.parent / "lib")))
    result = json.loads((work / "lagrange.results.json").read_text())
    actual = result["params"][0]
    if actual["dispersion"] != d or actual["extinction"] != e:
        raise RuntimeError("Lagrange-ng actual parameters differ from the requested scalar case")
    # JSON gives joint root-state log-likelihood and its normalized posterior.
    # L = joint / posterior recovers the full precision total; stdout is rounded.
    states = result["node-results"][0]["states"]
    totals = [row["llh"] - math.log(row["ratio"]) for row in states if row["ratio"] > 0]
    if max(totals) - min(totals) > 1e-12:
        raise RuntimeError("inconsistent Lagrange root-state normalizers")
    return record, float(totals[0]), result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, default=STAGE / "runtime")
    args = parser.parse_args()
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    rev = args.runtime / "revbayes-v1.4.2/revbayes-v1.4.2/bin/rb"
    lag = args.runtime / "lagrange-ng-v0.7.2-7/lagrange-ng"
    rows = []; native_c = None
    for i, (d, e) in enumerate(((.2, .1), (.05, .02), (.4, .2)), 1):
        work = output / f"case_{i}"; work.mkdir(exist_ok=True)
        (work / "tree.newick").write_text(TREE + "\n"); nexus(work / "geography.nex", GEOG)
        execution, values = revbayes_case(rev, work, d, e, audit_c=i == 1)
        if i == 1:
            native_c = values["C"]
        rust, q, c, root = rust_model(d, e)
        score_a = rust_model(d, e, contract="tip_survival")[0]
        lag_execution, lag_score, lag_result = lagrange_case(lag, work, d, e)
        variant, q_variant, c_variant, root_variant = rust_model(d, e, vicariance=2, root="legacy_flat_weights")
        # Source-defined Lagrange split generator: singleton order causes two
        # occurrences of each AB -> A|B event. This is deliberately NOT deduped.
        lag_c = np.zeros((4, 4, 4)); lag_c[1, 1, 1] = lag_c[2, 2, 2] = 1
        for x in (1, 2):
            for left, right in ((x, 3), (3, x), (x, 3 ^ x), (3 ^ x, x)):
                lag_c[3, left, right] += 1 / 8
        np.savez(work / "model_components.npz", Q=q, C_standard=c, C_revbayes_runtime=native_c,
                 C_lagrange_source=lag_c, C_rust_v2=c_variant, root_revbayes=root, root_lagrange=root_variant,
                 P_include=values["Include"]["P"], P_branch_conditional=values["CondSurv"]["P"])
        row = dict(case=i, d=d, e=e, tree=TREE, states=list(BITS),
            revbayes=dict(configuration='Include; pattern weights s=a; simplex(0,1,1,1); Contract-C',
                lnL=values["Include"]["lnL"], rust_lnL=rust, max_Q_error=float(np.max(abs(q-values["Include"]["Q"]))),
                max_C_error=float(np.max(abs(c-native_c))), max_P_error=float(np.max(abs(expm(q)-values["Include"]["P"]))),
                scalar_error=abs(values["Include"]["lnL"]-rust), independent_error=abs(independent_score(q,c,root)-rust),
                default_CondSurv_lnL=values["CondSurv"]["lnL"], rust_whole_tree_A_lnL=score_a,
                default_CondSurv_is_contract_A=False, default_vs_A_error=abs(values["CondSurv"]["lnL"]-score_a)),
            lagrange_ng=dict(configuration='official checkpoint evaluation; root raw ones; null C row zero',
                parameter_injection='official .ckp resume; actual params verified in result JSON',
                lnL=lag_score, requested_params_verified=True, actual_params=lag_result['params'],
                standard_DEC_status='unsupported_C_mismatch', max_C_error_vs_standard=float(np.max(abs(lag_c-c))),
                standard_DEC_lnL=rust_model(d,e,root='legacy_flat_weights')[0],
                explicit_variant='two areas only: Rust legacy custom v=2, s=y=1, all other DEC values',
                variant_api='historical scalar Contract-C; not an unlocked Tier-1b uncertainty model',
                native_state_to_rust_permutation=[0,2,1,3],
                native_root_weights=[1,1,1,1],effective_root_weights=root_variant.tolist(),
                root_equivalence='root null conditional likelihood is identically zero because C[0,:,:]=0; nonempty raw weights are exactly one in both engines',
                variant_Q_error=float(np.max(abs(q-q_variant))), variant_C_error=float(np.max(abs(lag_c-c_variant))),
                variant_lnL=variant, variant_scalar_error=abs(lag_score-variant),
                independent_lagrange_error=abs(lag_score-independent_score(q,lag_c,root_variant))),
            execution={**execution,'lagrange_ng':lag_execution},speedup=None)
        row['revbayes']['aligned']=all(row['revbayes'][key]<=tol for key,tol in
            [('max_Q_error',1e-12),('max_C_error',1e-12),('max_P_error',1e-12),('scalar_error',1e-8)])
        row['lagrange_ng']['explicit_variant_aligned']=row['lagrange_ng']['variant_C_error']<=1e-12 and row['lagrange_ng']['variant_scalar_error']<=1e-8
        rows.append(row); write_json(output/'progress.json',rows)
    report=dict(status='small_scalar_adapters_verified' if all(r['revbayes']['aligned'] and r['lagrange_ng']['explicit_variant_aligned'] for r in rows) else 'alignment_failed',
        versions=dict(revbayes='v1.4.2',lagrange_ng='v0.7.2-7'),cases=rows,
        standard_DEC_lagrange_comparison='unsupported; no default-DEC speedup may be inferred',
        full_benchmark_completed=False,source_sha256=sha(__file__),
        root_measure='nonempty weights matched; the null root contribution is identically zero under C; no additive score correction',
        empty_clado='zero row in both external kernels; Contract-B completion is not claimed',
        survival='RevBayes CondSurv conditions each branch, not the whole-tree leaf event; Contract-A is unsupported by this adapter')
    write_json(output/'results.json',report)
    print(json.dumps(dict(status=report['status'],max_revbayes_error=max(r['revbayes']['scalar_error'] for r in rows),
        max_lagrange_explicit_variant_error=max(r['lagrange_ng']['variant_scalar_error'] for r in rows),output=str(output))))
    if report['status']!='small_scalar_adapters_verified':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
