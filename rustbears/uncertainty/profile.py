"""Profile likelihood/support curves with explicit boundaries and calibration."""
from dataclasses import dataclass, field
import time

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2

from .bootstrap import _fit_options, _initial
from .checkpoint import CheckpointStore, model_manifest, task_fingerprint as input_fingerprint
from .contracts import ContractMismatchError, resolve_simulation_conditioning
from .guards import clone_model, validate_compatible_fits
from .parallel import parallel_scope, positive_int, resolve_parallel
from .rng import TaskIdentity, derive_task_seeds, task_fingerprint
from .simulate import PreparedSimulation

_NONNEGATIVE={"d","e","a","j"}


class ProfileOptimizationError(RuntimeError):
    """A profile ordinate or endpoint could not be certified numerically."""


@dataclass
class ParameterInterval:
    param: str
    interval_kind: str
    mle: float
    lower: float | None
    upper: float | None
    status_lower: str
    status_upper: str
    onesided_upper_from_zero: float | None
    se_hessian_local: float | None
    endpoint_deviance_residual_max: float
    n_evals: int
    profile_curve_df: pd.DataFrame
    critical_value: float = 0.
    disconnected_support_detected: bool = False
    bootstrap_plugin_df: pd.DataFrame | None = None


@dataclass
class ProfileCIResult:
    level: float
    method: str
    interval_kind: str
    likelihood_conditioning: str
    precision_mode: str
    obj_mle: float
    rebased_mle: bool
    updated_fit: object | None
    intervals: dict
    hessian_corr: pd.DataFrame | None
    total_evals: int
    elapsed_sec: float
    interpretation: str = "Connected profile component containing the fitted maximum."

    def to_dataframe(self):
        rows=[]
        for name,p in self.intervals.items():
            rows.append({key:getattr(p,key) for key in (
                "param","interval_kind","mle","lower","upper","status_lower","status_upper",
                "onesided_upper_from_zero","se_hessian_local","endpoint_deviance_residual_max",
                "critical_value","disconnected_support_detected","n_evals")})
        return pd.DataFrame(rows)

    def plot(self,filepath=None,ncols=2):
        import matplotlib.pyplot as plt
        ncols=min(positive_int(ncols,"ncols"),len(self.intervals))
        nrows=(len(self.intervals)+ncols-1)//ncols
        fig,axes=plt.subplots(nrows,ncols,squeeze=False,figsize=(5*ncols,3.5*nrows))
        for ax,(name,p) in zip(axes.flat,self.intervals.items()):
            df=p.profile_curve_df.sort_values("value")
            ax.plot(df.value,df.delta_deviance,color="#436c85")
            ax.axhline(p.critical_value,color="#a3412d",ls="--")
            ax.axvline(p.mle,color=".5",ls=":")
            ax.set(xlabel=name,ylabel="Twice profile score drop",title=p.interval_kind.replace("_"," "))
        for ax in list(axes.flat)[len(self.intervals):]: ax.set_visible(False)
        fig.tight_layout()
        if filepath is not None: fig.savefig(filepath,bbox_inches="tight")
        return fig


@dataclass
class _PointTask:
    model: object
    param: str
    value: float
    initial: dict
    options: dict
    threads: int
    seed: int = 1


def _profile_point(task):
    obj=task.model.model.with_parameters(**{task.param:dict(type="fixed",init=task.value)})
    model=clone_model(task.model,model_object=obj)
    # A strictly interior nuisance point exposes all transitions permitted by
    # this constrained model. Zero likelihood here is a support statement,
    # not an optimizer's failure being silently relabelled as an endpoint.
    support={}
    for key in model.model.free:
        p=model.model.params[key]
        support[key]=(p.min+p.max)/2
    support_ll=model.evaluate(support,calc_ancprobs=False,num_threads=task.threads).lnL
    if np.isneginf(support_ll):
        return dict(param=task.param,value=task.value,profile_obj=-np.inf,converged=True,
                    structural_zero=True,params=model.model.resolve(support),n_evals=1)
    if not np.isfinite(support_ll):
        raise ProfileOptimizationError("nonfinite numerical support probe")
    start=_initial(model,task.initial)
    try:
        fit=model.fit(**_fit_options(task.options,task.threads,task.seed,start))
    except (ValueError,RuntimeError,FloatingPointError):
        fit=model.fit(**_fit_options(task.options,task.threads,task.seed,support,2))
    if not fit.converged or not np.isfinite(fit.lnL):
        fit=model.fit(**_fit_options(task.options,task.threads,task.seed,fit.params,2))
    if not fit.converged or not np.isfinite(fit.lnL):
        raise ProfileOptimizationError(f"profile optimization failed for {task.param}={task.value}")
    return dict(param=task.param,value=task.value,profile_obj=fit.lnL,converged=True,
                structural_zero=False,params=fit.params,n_evals=fit.n_evals+1)


class _Curve:
    def __init__(self,fit,param,options,pc,ordered_map,cache,commit):
        self.fit=fit;self.param=param;self.options=options;self.pc=pc
        self.ordered_map=ordered_map;self.cache=cache;self.commit=commit;self.used={}
        # A rebase changes the nuisance starting basin even when the parameter
        # box is unchanged. Reusing the previous basin's ordinate can splice a
        # stale grid score into newly optimized root-search scores and create
        # an artificial discontinuity. Checkpoint generations must therefore
        # identify both the model and its current reference fit.
        self.generation=input_fingerprint(dict(model=fit.model.model,
                                               reference_params=fit.params,
                                               reference_score=fit.lnL))

    def points(self,values):
        tasks=[];keys=[]
        for value in values:
            value=float(value);key=(self.generation,self.param,value.hex())
            if key not in self.cache and key not in keys:
                tasks.append(_PointTask(self.fit.model,self.param,value,self.fit.params,self.options,self.pc.inner_threads))
                keys.append(key)
        for key,row in zip(keys,self.ordered_map(_profile_point,tasks)):
            row["generation"]=self.generation
            self.cache[key]=row;self.commit(dict(kind="point",**row))
        rows=[]
        for value in values:
            key=(self.generation,self.param,float(value).hex());self.used[key]=self.cache[key];rows.append(self.cache[key])
        return rows

    def point(self,value): return self.points([value])[0]

    def deviance(self,value): return 2*(self.fit.lnL-self.point(value)["profile_obj"])

    def frame(self):
        rows=[]
        for row in sorted(self.used.values(),key=lambda r:r["value"]):
            rows.append({**{k:v for k,v in row.items() if k!="params"},
                         "delta_deviance":max(0.,2*(self.fit.lnL-row["profile_obj"])),
                         **{"nuisance_"+k:v for k,v in row["params"].items() if k!=self.param}})
        frame=pd.DataFrame(rows)
        first=["value","profile_obj","delta_deviance","converged","structural_zero","param","n_evals","generation"]
        return frame.reindex(columns=first+sorted(set(frame.columns)-set(first)))


def _axis_grid(mle,bound,n_grid):
    if bound==mle: return np.array([mle])
    if bound>mle and mle>=0 and bound>100*max(mle,1e-8):
        scale=max(mle,1e-6)
        grid=mle+np.expm1(np.linspace(0,np.log1p((bound-mle)/scale),n_grid))*scale
    elif bound==0 and mle>0:
        grid=np.r_[mle,np.geomspace(mle*.8,max(mle*1e-8,1e-15),max(2,n_grid-2)),0.]
    else:
        grid=np.linspace(mle,bound,n_grid)
    grid[0]=mle;grid[-1]=bound
    return np.unique(grid)[::1 if bound>mle else -1]


def _endpoint(curve,grid,critical,tol,physical_zero=False,side="upper"):
    rows=curve.points(grid)
    drops=np.array([2*(curve.fit.lnL-row["profile_obj"])-critical for row in rows])
    outside=np.flatnonzero(drops>0)
    if not len(outside):
        if physical_zero and grid[-1]==0:
            return 0.,"at_physical_zero",0.,False
        return None,"censored_at_box_"+side,0.,False
    i=int(outside[0])
    if i==0:
        if drops[0]>max(tol,1e-6):
            raise ProfileOptimizationError("profile at fitted maximum is below the supplied score")
        return float(grid[0]),"interior_bracketed",abs(float(drops[0])),False
    left,right=sorted((float(grid[i-1]),float(grid[i])))
    fn=lambda x:curve.deviance(x)-critical
    root=brentq(fn,left,right,xtol=5e-15,rtol=1e-14,maxiter=200)
    residual=abs(fn(root))
    if residual>=tol:
        raise ProfileOptimizationError(f"endpoint deviance residual {residual:g} exceeds {tol:g}")
    disconnected=bool(np.any(drops[i+1:]<=0))
    status="infeasible_zero_approached" if physical_zero and rows[-1]["structural_zero"] else "interior_bracketed"
    return float(root),status,float(residual),disconnected


def _local_hessian(fit,params):
    names=list(fit.model.model.free)
    if not names: return {},None
    theta=np.array([fit.params[k] for k in names])
    h=np.maximum(np.abs(theta)*1e-3,1e-5)
    bounds=np.array([(fit.model.model.params[k].min,fit.model.model.params[k].max) for k in names])
    if np.any(theta-h<=bounds[:,0]) or np.any(theta+h>=bounds[:,1]): return {},None
    with fit.model.prepared_likelihood(fit.params) as prepared:
        def objective(x):
            p=fit.model.model.resolve(dict(zip(names,x)))
            return -fit.model._evaluate(prepared.context,p,False)
        f0=objective(theta);H=np.empty((len(names),len(names)))
        for i in range(len(names)):
            ei=np.zeros(len(names));ei[i]=h[i]
            H[i,i]=(objective(theta+ei)-2*f0+objective(theta-ei))/h[i]**2
            for j in range(i):
                ej=np.zeros(len(names));ej[j]=h[j]
                H[i,j]=H[j,i]=(objective(theta+ei+ej)-objective(theta+ei-ej)-
                               objective(theta-ei+ej)+objective(theta-ei-ej))/(4*h[i]*h[j])
    if not np.isfinite(H).all() or np.linalg.eigvalsh(H).min()<=0: return {},None
    cov=np.linalg.inv(H);se=np.sqrt(np.diag(cov))
    return dict(zip(names,se)),pd.DataFrame(cov/np.outer(se,se),index=names,columns=names)


@dataclass
class _PluginTask:
    model: object
    geog: object
    param: str
    fixed_value: float
    initial: dict
    options: dict
    threads: int
    seeds: tuple
    rep: int
    rng: dict


def _plugin_replicate(task):
    model=clone_model(task.model,geog=task.geog,
                      allow_null_tips=task.model.config.allow_null_tips or any(not v for v in task.geog.ranges.values()))
    full=model.fit(**_fit_options(task.options,task.threads,task.seeds[0],_initial(model,task.initial)))
    point=_profile_point(_PointTask(model,task.param,task.fixed_value,full.params,task.options,task.threads,task.seeds[1]))
    if point["profile_obj"]>full.lnL+1e-6:
        full=model.fit(**_fit_options(task.options,task.threads,task.seeds[0],_initial(model,point["params"]),2))
    if not full.converged or not np.isfinite(full.lnL) or point["profile_obj"]>full.lnL+5e-6:
        raise ProfileOptimizationError("bootstrap-plugin unrestricted fit did not dominate its constraint")
    return dict(rep_id=task.rep,d_stat=max(0.,2*(full.lnL-point["profile_obj"])),
                obj_full=full.lnL,obj_constrained=point["profile_obj"],params_full=full.params,
                params_constrained=point["params"],n_evals=full.n_evals+point["n_evals"],rng=task.rng)


def _plugin_threshold(fit,param,n_boot,level,resolved,seed,options,pc,ordered_map,cache,commit):
    generation=input_fingerprint(dict(params=fit.params,lnL=fit.lnL))
    identities={(rep,role):TaskIdentity("profile_plugin",fit.model.model.name,"single_tree",rep,
                                      role=param+"_"+role)
                for rep in range(1,n_boot+1) for role in ("simulate","full","constrained")}
    seeds=derive_task_seeds(seed,identities.values())
    tasks=[]
    with PreparedSimulation(fit,conditioning=resolved) as sampler:
        for rep in range(1,n_boot+1):
            key=(param,generation,rep)
            if key in cache: continue
            keys=[identities[(rep,role)] for role in ("simulate","full","constrained")]
            batch=sampler.draw(1,start_rep=rep,seed=seeds[keys[0]],num_threads=pc.inner_threads)
            rng={role:dict(identity=t.to_dict(),seed=seeds[t],fingerprint=task_fingerprint(seed,t))
                 for role,t in zip(("simulate","full","constrained"),keys)}
            tasks.append(_PluginTask(fit.model,batch.to_tip_ranges(0),param,fit.params[param],fit.params,
                                     options,pc.inner_threads,tuple(seeds[t] for t in keys[1:]),rep,rng))
    # A failed plug-in calibration is an error; no quantile is computed after
    # silently dropping the hardest simulated datasets.
    for row in ordered_map(_plugin_replicate,tasks):
        cache[(param,generation,row["rep_id"])]=row
        commit(dict(kind="plugin",param=param,generation=generation,**row))
    frame=pd.DataFrame([cache[(param,generation,rep)] for rep in range(1,n_boot+1)]).sort_index(axis=1)
    return float(np.quantile(frame.d_stat,level,method="higher")),frame


def profile_ci(fit,params=None,level=.95,method="asymptotic_chi2",n_boot_plugin=199,
               include_onesided_zero_ref=False,eval_exact_zero=True,multi_basin_guard=True,n_grid=16,
               tol=1e-4,num_threads=1,seed=12345,*,simulation_conditioning="from_model",
               checkpoint_file=None,resume=False,parallel=None,fit_options=None):
    """Return the connected profile component containing the maximum.

    A/B use a regular-asymptotic reference; C always reports support intervals.
    ``bootstrap_plugin`` calibrates at the global point estimate, not by
    pointwise test inversion. Neither mode promises boundary/weak-identification
    coverage. ``onesided_upper_from_zero`` is an explicitly requested reference
    using chi-square(1, 2*level-1), only when the exact zero is feasible.
    """
    started=time.monotonic();validate_compatible_fits((fit,));fit.model._require_tier1a()
    if not 0<level<1: raise ValueError("level must lie in (0,1)")
    if include_onesided_zero_ref and level<=.5: raise ValueError("one-sided reference requires level > .5")
    if method not in ("asymptotic_chi2","bootstrap_plugin"): raise ValueError("unknown profile method")
    n_grid=positive_int(n_grid,"n_grid")
    if n_grid<3: raise ValueError("n_grid must be at least 3")
    if not np.isfinite(tol) or tol<=0: raise ValueError("tol must be positive and finite")
    names=list(fit.model.model.free if params is None else params)
    if not names or len(names)!=len(set(names)) or any(k not in fit.model.model.free for k in names):
        raise ValueError("params must be distinct free parameter names")
    pc=resolve_parallel(parallel,num_threads);options=dict(fit_options or {})
    if set(options)-{"optimizer","maxiter"}: raise ValueError("fit_options accepts optimizer and maxiter")
    resolved=None
    if method=="bootstrap_plugin":
        n_boot_plugin=positive_int(n_boot_plugin,"n_boot_plugin")
        resolved=resolve_simulation_conditioning(fit.likelihood_conditioning,simulation_conditioning)
        if fit.likelihood_conditioning!="legacy_unconditioned" and resolved!=fit.likelihood_conditioning:
            raise ContractMismatchError("profile-plugin simulation must match the probability contract")
    kind=("profile_support_interval" if fit.likelihood_conditioning=="legacy_unconditioned" else
          "bootstrap_plugin_interval" if method=="bootstrap_plugin" else "asymptotic_confidence_interval")
    manifest=dict(kind="profile_ci.v2",fit=model_manifest(fit),params=names,level=level,method=method,
                  n_boot_plugin=n_boot_plugin,include_onesided_zero_ref=include_onesided_zero_ref,
                  eval_exact_zero=eval_exact_zero,multi_basin_guard=multi_basin_guard,n_grid=n_grid,tol=tol,
                  seed=int(seed),simulation_conditioning=resolved,fit_options=options)
    store=CheckpointStore(checkpoint_file,manifest,resume=resume) if checkpoint_file else None
    if resume and store is None: raise ValueError("resume requires checkpoint_file")
    point_cache={};plugin_cache={}
    if store is not None:
        for row in store.records:
            if row["kind"]=="point":
                point_cache[(row["generation"],row["param"],float(row["value"]).hex())]={k:v for k,v in row.items() if k not in ("id","kind")}
            elif row["kind"]=="plugin":
                plugin_cache[(row["param"],row["generation"],row["rep_id"])]= {
                    k:v for k,v in row.items() if k not in ("kind","param","generation","id")}
    def commit(record):
        if store is not None:
            store.append(dict(id=max(store.completed_ids,default=0)+1,**record))
    original=fit;rebased=False;intervals={}
    with parallel_scope(pc) as ordered_map:
        for restart in range(4):
            curves={};candidate=None;intervals={}
            for name in names:
                from ..model import model_template
                curve=_Curve(fit,name,options,pc,ordered_map,point_cache,commit);curves[name]=curve
                p=fit.model.model.params[name];mle=fit.params[name]
                default_min=model_template(fit.model.model.name).params[name].min
                custom_positive=name in _NONNEGATIVE and p.min>0 and p.min!=default_min
                physical=bool(eval_exact_zero and name in _NONNEGATIVE and p.min>=0 and not custom_positive)
                if eval_exact_zero and custom_positive:
                    zero=curve.point(0.)
                    if zero["profile_obj"]>fit.lnL+1e-6:
                        raise ProfileOptimizationError("exact zero improves a fit with a custom positive lower bound; set min=0 or eval_exact_zero=False")
                lower=0. if physical else p.min
                grids=(_axis_grid(mle,lower,n_grid),_axis_grid(mle,p.max,n_grid))
                curve.points(np.r_[grids[0],grids[1]])
                best=max(curve.used.values(),key=lambda row:row["profile_obj"])
                if best["profile_obj"]>fit.lnL+1e-6:
                    if not multi_basin_guard:
                        raise ProfileOptimizationError("profile exceeds fitted maximum; enable multi_basin_guard")
                    candidate=best;break
                calibration=None
                critical=float(chi2.ppf(level,1))
                if method=="bootstrap_plugin":
                    critical,calibration=_plugin_threshold(fit,name,n_boot_plugin,level,resolved,seed,options,
                                                          pc,ordered_map,plugin_cache,commit)
                low=_endpoint(curve,grids[0],critical,tol,physical,side="lower")
                high=_endpoint(curve,grids[1],critical,tol)
                one=None
                if include_onesided_zero_ref and physical and not curve.point(0.)["structural_zero"]:
                    one=_endpoint(curve,grids[1],float(chi2.ppf(2*level-1,1)),tol)[0]
                # Adaptive endpoint searches can discover a basin missed by
                # the initial grid; the same global rebase rule applies.
                best=max(curve.used.values(),key=lambda row:row["profile_obj"])
                if best["profile_obj"]>fit.lnL+1e-6:
                    candidate=best;break
                frame=curve.frame()
                intervals[name]=ParameterInterval(name,kind,mle,low[0],high[0],low[1],high[1],one,None,
                    max(low[2],high[2]),int(frame.n_evals.sum()),frame,critical,low[3] or high[3],calibration)
            if candidate is None: break
            if not multi_basin_guard:
                raise ProfileOptimizationError("profile exceeds fitted maximum; enable multi_basin_guard")
            obj=fit.model.model
            changes={k:dict(min=0.) for k in obj.free if k in _NONNEGATIVE and candidate["params"][k]==0 and obj.params[k].min>0}
            model=clone_model(fit.model,model_object=obj.with_parameters(**changes))
            new=model.fit(**_fit_options(options,pc.inner_threads,seed,_initial(model,candidate["params"]),2))
            if not new.converged or new.lnL<candidate["profile_obj"]-1e-6:
                raise ProfileOptimizationError("global rebase failed to recover the improved profile point")
            fit=new;rebased=True
        else:
            raise ProfileOptimizationError("profile maximum continued changing after three rebases")
        se,corr=_local_hessian(fit,names)
    for name,interval in intervals.items(): interval.se_hessian_local=se.get(name)
    total=sum(int(row["n_evals"]) for row in point_cache.values())+sum(int(row["n_evals"]) for row in plugin_cache.values())
    return ProfileCIResult(level,method,kind,fit.likelihood_conditioning,fit.precision_mode,fit.lnL,rebased,
                           fit if rebased else None,intervals,corr,total,time.monotonic()-started)


def profile_surface(fit,param_x,param_y,grid_x=25,grid_y=25,reoptimize_nuisance=True,num_threads=1,*,parallel=None):
    """Evaluate a two-parameter surface, refitting all remaining free parameters."""
    validate_compatible_fits((fit,));fit.model._require_tier1a()
    if param_x==param_y or any(k not in fit.model.model.free for k in (param_x,param_y)):
        raise ValueError("surface axes must be two distinct free parameters")
    def grid(values,name):
        p=fit.model.model.params[name]
        if isinstance(values,(int,np.integer)) and not isinstance(values,(bool,np.bool_)):
            return np.linspace(p.min,p.max,positive_int(values,"grid size"))
        out=np.asarray(values,dtype=float)
        if out.ndim!=1 or not len(out) or not np.isfinite(out).all() or np.any(out<p.min) or np.any(out>p.max):
            raise ValueError("surface grid must be a finite vector within the parameter bounds")
        return out
    tasks=[];positions=[];pc=resolve_parallel(parallel,num_threads)
    for x in grid(grid_x,param_x):
        for y in grid(grid_y,param_y):
            changes={param_x:dict(type="fixed",init=float(x))}
            if not reoptimize_nuisance:
                changes.update({k:dict(type="fixed",init=fit.params[k]) for k in fit.model.model.free if k not in (param_x,param_y)})
            model=clone_model(fit.model,model_object=fit.model.model.with_parameters(**changes))
            tasks.append(_PointTask(model,param_y,float(y),fit.params,{},pc.inner_threads))
            positions.append((float(x),float(y)))
    with parallel_scope(pc) as ordered_map:
        rows=[{param_x:x,param_y:y,"profile_obj":row["profile_obj"],
               "delta_deviance":2*(fit.lnL-row["profile_obj"]),"converged":row["converged"],
               **{"nuisance_"+k:v for k,v in row["params"].items() if k not in (param_x,param_y)}}
              for (x,y),row in zip(positions,ordered_map(_profile_point,tasks))]
    return pd.DataFrame(rows)
