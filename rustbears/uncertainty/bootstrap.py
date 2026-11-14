"""Parametric bootstrap of nested-model scores with complete failure accounting."""
from dataclasses import dataclass, field, replace
import time

import numpy as np
import pandas as pd
from scipy.stats import beta, chi2

from .checkpoint import CheckpointStore, model_manifest
from .contracts import ContractMismatchError, resolve_simulation_conditioning
from .guards import clone_model, prepare_comparison_fits, validate_compatible_fits
from .parallel import parallel_scope, positive_int, resolve_parallel
from .rng import TaskIdentity, derive_task_seeds, task_fingerprint
from .simulate import PreparedSimulation


def failure_probability_intervals(exceedances, failed, planned, confidence=.95):
    """Algebraic missing-outcome bounds and a CP envelope for the original tail.

    The Monte Carlo interval is [CP_L(r,B), CP_U(r+failed,B)]. It is not an
    ordinary two-sided interval treating failed replicates as observed hits.
    """
    planned=positive_int(planned,"planned")
    if any(isinstance(x,(bool,np.bool_)) or not isinstance(x,(int,np.integer)) for x in (exceedances,failed)):
        raise TypeError("exceedances and failed must be integers")
    if not 0 <= exceedances <= planned or not 0 <= failed <= planned-exceedances:
        raise ValueError("invalid exceedance or failure counts")
    alpha=1-confidence
    if not 0 < alpha < 1:
        raise ValueError("confidence must lie in (0,1)")
    r=int(exceedances); upper=r+int(failed)
    bounds=((r+1)/(planned+1),(upper+1)/(planned+1))
    cp=(0. if r==0 else float(beta.ppf(alpha/2,r,planned-r+1)),
        1. if upper==planned else float(beta.ppf(1-alpha/2,upper+1,planned-upper)))
    return bounds,cp


def _nested_model(null,alt):
    """Verify parameter nesting and include an exact physical zero if required."""
    from ..model import model_template
    n,a=null.model,alt.model
    default=model_template(a.name)
    for name in ("dist","env","manual","area_size","times","masks"):
        if not np.array_equal(np.asarray(getattr(null,name)),np.asarray(getattr(alt,name))):
            raise ValueError(f"nested models must use identical resolved {name} inputs")
    if n.name.split("+")[0] != a.name.split("+")[0] or n.params.keys()!=a.params.keys():
        raise ValueError("bootstrap_lrt requires nested models with the same base family")
    changes={}
    extra=[]
    for key,p in n.params.items():
        q=a.params[key]
        if p.type=="free" and p.min==p.max: p=replace(p,type="fixed",init=p.min)
        if q.type=="free" and q.min==q.max: q=replace(q,type="fixed",init=q.min)
        if p.type not in ("fixed","free") or q.type not in ("fixed","free"):
            if p.type!=q.type:
                raise ValueError(f"non-nested parameter link: {key}")
            continue
        if q.type=="fixed":
            if p.type!="fixed" or p.init!=q.init:
                raise ValueError(f"alternative fixes null-model parameter {key}")
        elif p.type=="free":
            if p.min<q.min or p.max>q.max:
                raise ValueError(f"alternative bounds do not contain null bounds for {key}")
        else:
            if p.init==0 and q.min>0 and key in ("d","e","a","j"):
                if q.min!=default.params[key].min:
                    raise ValueError(f"custom positive lower bound excludes the null physical zero for {key}")
                changes[key]=dict(min=0.)
            elif not q.min<=p.init<=q.max:
                raise ValueError(f"null value outside alternative bounds for {key}")
            if q.min<q.max:
                extra.append(key)
    if not extra:
        raise ValueError("alternative must add at least one free parameter")
    return a.with_parameters(**changes),tuple(extra)


def _initial(model, values):
    return {key:float(np.clip(values.get(key,p.init),p.min,p.max))
            for key,p in model.model.params.items() if p.type=="free"}


def _fit_options(options,threads,seed,initial,maxiter_multiplier=1):
    out=dict(options)
    out.update(num_threads=threads,seed=np.random.default_rng(seed),initial=initial)
    out["maxiter"]=positive_int(out.get("maxiter",500),"maxiter")*maxiter_multiplier
    return out


class _PairFitError(RuntimeError):
    def __init__(self,message,null_params,alt_params):
        super().__init__(message)
        self.null_params=null_params;self.alt_params=alt_params


def _fit_snapshot(fit):
    """Diagnostic facts only; these never decide whether a replicate is usable."""
    return dict(obj=float(fit.lnL),converged=bool(fit.converged),
                params=dict(fit.params),n_evals=int(fit.n_evals),message=str(fit.message))


def _fit_pair(null_model,alt_model,null_start,alt_start,options,threads,seeds,*,retry=False,null_candidate=None,audit=None):
    """Always fit both models, then inject the fitted null into the alternative."""
    mult=2 if retry else 1
    audit={} if audit is None else audit
    audit.update(raw_obj_null=np.nan,raw_obj_alt=np.nan,raw_d_stat=np.nan,alternative_fits=[])
    nf=null_model.fit(**_fit_options(options,threads,seeds[0],_initial(null_model,null_start),mult))
    audit["null_fit"]=_fit_snapshot(nf)
    if null_candidate is not None and null_candidate.converged and null_candidate.lnL>nf.lnL:
        nf=null_candidate
        audit["retained_null_candidate"]=_fit_snapshot(nf)
    audit["raw_obj_null"]=float(nf.lnL)
    candidates=[]
    errors=[]
    for source,values in (("alternative_start",alt_start),("injected_null_start",{**alt_start,**nf.params})):
        try:
            candidate=alt_model.fit(**_fit_options(options,threads,seeds[1],_initial(alt_model,values),mult))
            candidates.append(candidate)
            audit["alternative_fits"].append(dict(start_source=source,**_fit_snapshot(candidate)))
        except (ValueError,RuntimeError,FloatingPointError) as error:
            errors.append(str(error))
            audit["alternative_fits"].append(dict(start_source=source,error=f"{type(error).__name__}: {error}"))
    eligible=[f for f in candidates if f.converged and np.isfinite(f.lnL)]
    if not eligible:
        best=max(candidates,key=lambda f:f.lnL,default=None)
        if best is not None:
            audit.update(raw_obj_alt=float(best.lnL),raw_d_stat=float(2*(best.lnL-nf.lnL)))
        raise _PairFitError("alternative did not converge: "+"; ".join(errors or [f.message for f in candidates]),
                            nf.params,best.params if best is not None else alt_start)
    af=max(eligible,key=lambda fit:fit.lnL)
    audit.update(raw_obj_alt=float(af.lnL),raw_d_stat=float(2*(af.lnL-nf.lnL)))
    if not nf.converged or not np.isfinite(nf.lnL):
        raise _PairFitError("null model did not converge: "+nf.message,nf.params,af.params)
    if 2*(af.lnL-nf.lnL)<-1e-5:
        raise _PairFitError("negative nested statistic after null-start injection",nf.params,af.params)
    return nf,af


@dataclass
class _ReplicateTask:
    rep_id: int
    phase: str
    null_model: object
    alt_model: object
    geog: object
    null_start: dict
    alt_start: dict
    fit_options: dict
    inner_threads: int
    seeds: tuple
    rng: dict


def _fit_replicate(task):
    # No compiled context or native handle crosses the task boundary.
    nm=clone_model(task.null_model,geog=task.geog,
                   allow_null_tips=task.null_model.config.allow_null_tips or any(not v for v in task.geog.ranges.values()))
    am=clone_model(task.alt_model,geog=task.geog,
                   allow_null_tips=task.alt_model.config.allow_null_tips or any(not v for v in task.geog.ranges.values()))
    errors=[]; fits=None; attempts=[]
    null_start=dict(task.null_start);alt_start=dict(task.alt_start)
    for attempt in range(2):
        audit=dict(attempt=attempt+1,maxiter_multiplier=2 if attempt else 1,
                   initial_null=dict(null_start),initial_alt=dict(alt_start))
        try:
            fits=_fit_pair(nm,am,null_start,alt_start,task.fit_options,
                           task.inner_threads,task.seeds,retry=bool(attempt),audit=audit)
            audit.update(status="converged",error="")
            attempts.append(audit)
            break
        except (ValueError,RuntimeError,FloatingPointError) as error:
            errors.append(f"{type(error).__name__}: {error}")
            audit.update(status="failed",error=errors[-1])
            attempts.append(audit)
            if isinstance(error,_PairFitError):
                null_start=error.null_params;alt_start=error.alt_params
    record=dict(id=task.rep_id,rep_id=task.rep_id,phase=task.phase,
                status="converged" if fits is not None else "failed",retried=bool(errors),
                obj_null=np.nan,obj_alt=np.nan,d_stat=np.nan,error=" | ".join(errors),rng=task.rng,
                params_null=null_start,params_alt=alt_start,attempts=attempts,
                **{key:attempts[-1].get(key,np.nan) for key in ("raw_obj_null","raw_obj_alt","raw_d_stat")})
    arrays={}
    if fits is not None:
        nf,af=fits
        record.update(obj_null=nf.lnL,obj_alt=af.lnL,d_stat=max(0.,2*(af.lnL-nf.lnL)),
                      params_null=nf.params,params_alt=af.params,
                      n_evals_null=nf.n_evals,n_evals_alt=af.n_evals)
        for role,fit in zip(("null","alt"),fits):
            for suffix,values in (("ancestor_top",fit.evaluation.ancestral_probs_branch_top),
                                  ("ancestor_bottom",fit.evaluation.ancestral_probs_branch_bottom)):
                if values is not None:
                    arrays[f"{role}_{suffix}"]=values
    return record,arrays


@dataclass
class BootstrapLRTResult:
    null_model: str
    alt_model: str
    likelihood_conditioning: str
    recomputed_observed_under_new_contract: bool
    root_distribution: str
    precision_mode: str
    d_obs: float
    df_nominal: int
    pval_asymptotic: float
    n_planned: int
    n_converged: int
    n_retried_recovered: int
    n_failed: int
    pval_mc_conservative: float
    pval_failure_bounds: tuple
    pval_mc_ci95: tuple
    replicates_df: pd.DataFrame
    surrogate_power_fixed_threshold: float | None
    elapsed_sec: float
    simulation_conditioning: str = "tip_survival"
    n_exceedances: int = 0
    critical_value: float | None = None
    surrogate_replicates_df: pd.DataFrame | None = None
    surrogate_power_failure_bounds: tuple | None = None
    surrogate_power_mc_ci95: tuple | None = None
    observed_null_fit: object = field(default=None,repr=False)
    observed_alt_fit: object = field(default=None,repr=False)
    replicate_arrays: dict = field(default_factory=dict,repr=False)

    @property
    def effective_parameter_bounds(self):
        return {role:{name:(p.min,p.max) for name,p in fit.model.model.params.items() if p.type=="free"}
                for role,fit in (("null",self.observed_null_fit),("alt",self.observed_alt_fit))}

    def summary(self):
        return pd.DataFrame([dict(null_model=self.null_model,alt_model=self.alt_model,
            likelihood_conditioning=self.likelihood_conditioning,simulation_conditioning=self.simulation_conditioning,
            d_obs=self.d_obs,pval_mc_conservative=self.pval_mc_conservative,
            pval_failure_low=self.pval_failure_bounds[0],pval_failure_high=self.pval_failure_bounds[1],
            pval_mc_ci95_low=self.pval_mc_ci95[0],pval_mc_ci95_high=self.pval_mc_ci95[1],
            n_planned=self.n_planned,n_converged=self.n_converged,n_failed=self.n_failed,
            surrogate_power_fixed_threshold=self.surrogate_power_fixed_threshold)])

    def plot(self,filepath=None):
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots()
        values=self.replicates_df.loc[self.replicates_df.status=="converged","d_stat"]
        ax.hist(values,bins="auto",color="#436c85",alpha=.8)
        ax.axvline(self.d_obs,color="#a3412d",label=f"Observed D = {self.d_obs:.4g}")
        ax.set(xlabel="Bootstrap score difference (2 log units)",ylabel="Replicates",
               title=f"{self.n_failed}/{self.n_planned} failed; conservative p = {self.pval_mc_conservative:.4g}")
        ax.legend();fig.tight_layout()
        if filepath is not None: fig.savefig(filepath,bbox_inches="tight")
        return fig


def bootstrap_lrt(null_fit,alt_fit,n_boot=399,target_conditioning="match_fits",refit_observed=False,
                  root_distribution="model_prior",calc_surrogate_power=False,n_alt_sim=199,alpha=.05,
                  num_threads=8,seed=12345,checkpoint_file=None,*,simulation_conditioning="from_model",
                  resume=False,parallel=None,fit_options=None):
    """Calibrate a nested score difference under a fitted generating model.

    This plug-in test is not finite-sample exact. C scores require an explicit
    generating distribution. Failed replicates stay in the planned denominator;
    the MC interval envelopes both possible failure-tail assignments. Surrogate
    power uses one fixed null threshold, not nested-bootstrap test power.
    """
    started=time.monotonic(); n_boot=positive_int(n_boot,"n_boot")
    if not 0<alpha<1: raise ValueError("alpha must lie in (0,1)")
    if calc_surrogate_power: n_alt_sim=positive_int(n_alt_sim,"n_alt_sim")
    pc=resolve_parallel(parallel,num_threads)
    options=dict(fit_options or {})
    if set(options)-{"optimizer","maxiter"}:
        raise ValueError("fit_options accepts optimizer and maxiter; task seeds/initials/threads are managed")
    with parallel_scope(pc):
        prepared=prepare_comparison_fits(null_fit,alt_fit,target_conditioning=target_conditioning,
                                        refit_observed=refit_observed,
                                        fit_options={**options,"num_threads":pc.inner_threads})
    nf,af=prepared.fits
    nf.model._require_tier1a();af.model._require_tier1a()
    alt_object,extra=_nested_model(nf.model,af.model)
    am=clone_model(af.model,model_object=alt_object)
    # Declare the whole RNG task set before any dispatch; retries reuse keys.
    identities={}
    for phase,count in (("null",n_boot),("alt",n_alt_sim if calc_surrogate_power else 0)):
        for rep in range(1,count+1):
            for role,model_id in (("simulate",nf.model.model.name if phase=="null" else af.model.model.name),
                                  ("fit_null",nf.model.model.name),("fit_alt",af.model.model.name)):
                identities[(phase,rep,role)]=TaskIdentity("bootstrap_"+phase,model_id,"single_tree",rep,role=role)
    observed_ids=[TaskIdentity("bootstrap_observed",m.model.name,"single_tree",role=role)
                  for m,role in ((nf.model,"null"),(af.model,"alt"))]
    seeds=derive_task_seeds(seed,[*identities.values(),*observed_ids])
    # Observed data get the same nested-start safeguard as every replicate.
    with parallel_scope(pc):
        observed_nf,observed_af=_fit_pair(clone_model(nf.model),am,nf.params,af.params,options,pc.inner_threads,
                                         tuple(seeds[t] for t in observed_ids),null_candidate=nf)
    if nf.lnL>observed_nf.lnL and nf.converged: observed_nf=nf
    if af.lnL>observed_af.lnL and af.converged: observed_af=replace(af,model=am)
    if 2*(observed_af.lnL-observed_nf.lnL)<-1e-5:
        raise ValueError("observed alternative is below the null after nested-start injection")
    nf,af=observed_nf,observed_af
    d_obs=max(0.,2*(af.lnL-nf.lnL))
    resolved=resolve_simulation_conditioning(nf.likelihood_conditioning,simulation_conditioning)
    if nf.likelihood_conditioning!="legacy_unconditioned" and resolved!=nf.likelihood_conditioning:
        raise ContractMismatchError("bootstrap simulation must match the fitted probability contract")
    manifest=dict(kind="bootstrap_lrt.v2",null=model_manifest(nf),alt=model_manifest(af),
                  n_boot=n_boot,n_alt_sim=n_alt_sim if calc_surrogate_power else 0,alpha=alpha,
                  seed=int(seed),simulation_conditioning=resolved,root_distribution=root_distribution,
                  fit_options=options,recomputed=prepared.recomputed_observed_under_new_contract,
                  random_tasks=[t.to_dict() for t in identities.values()])
    store=CheckpointStore(checkpoint_file,manifest,resume=resume) if checkpoint_file else None
    if resume and store is None: raise ValueError("resume requires checkpoint_file")
    records={} if store is None else {r["id"]:r for r in store.records}
    arrays={} if store is None else {i:store.load_arrays(i) for i in store.completed_ids}
    for phase,count,generating in (("null",n_boot,nf),("alt",n_alt_sim if calc_surrogate_power else 0,af)):
        if not count: continue
        tasks=[]
        with PreparedSimulation(generating,conditioning=resolved,root_distribution=root_distribution) as sampler:
            for rep in range(1,count+1):
                rid=rep if phase=="null" else n_boot+rep
                if rid in records: continue
                keys=[identities[(phase,rep,role)] for role in ("simulate","fit_null","fit_alt")]
                batch=sampler.draw(1,start_rep=rep,seed=seeds[keys[0]],num_threads=pc.inner_threads)
                rng={role:dict(identity=t.to_dict(),seed=seeds[t],fingerprint=task_fingerprint(seed,t))
                     for role,t in zip(("simulate","fit_null","fit_alt"),keys)}
                tasks.append(_ReplicateTask(rid,phase,nf.model,am,batch.to_tip_ranges(0),nf.params,af.params,
                                            options,pc.inner_threads,tuple(seeds[t] for t in keys[1:]),rng))
        with parallel_scope(pc) as ordered_map:
            for record,values in ordered_map(_fit_replicate,tasks):
                if store is not None: store.append(record,values)
                records[record["id"]]=record;arrays[record["id"]]=values
    frame=pd.DataFrame([records[i] for i in range(1,n_boot+1)]).sort_index(axis=1)
    converged=frame.status.eq("converged")
    failures=int((~converged).sum())
    r=int((frame.loc[converged,"d_stat"]>=d_obs).sum())
    fb,mc=failure_probability_intervals(r,failures,n_boot)
    alt_frame=None;power=None;power_fb=None;power_ci=None;critical=None
    if calc_surrogate_power:
        reference=np.r_[frame.loc[converged,"d_stat"].to_numpy(),np.full(failures,np.inf)]
        critical=float(np.sort(reference)[int(np.ceil((1-alpha)*(n_boot-1)))])
        alt_frame=pd.DataFrame([records[n_boot+i] for i in range(1,n_alt_sim+1)]).sort_index(axis=1)
        alt_ok=alt_frame.status.eq("converged")
        hits=int((alt_frame.loc[alt_ok,"d_stat"]>critical).sum())
        afail=int((~alt_ok).sum())
        power=hits/n_alt_sim
        power_fb=(power,(hits+afail)/n_alt_sim)
        _,power_ci=failure_probability_intervals(hits,afail,n_alt_sim)
    return BootstrapLRTResult(nf.model.model.name,af.model.model.name,nf.likelihood_conditioning,
        prepared.recomputed_observed_under_new_contract,root_distribution,nf.precision_mode,d_obs,len(extra),
        float(chi2.sf(d_obs,len(extra))),n_boot,int(converged.sum()),int((converged & frame.retried).sum()),
        failures,fb[1],fb,mc,frame,power,time.monotonic()-started,resolved,r,critical,alt_frame,power_fb,power_ci,
        nf,af,arrays)
