"""Model configuration, compiled likelihood evaluation and optimization."""
import ast
from dataclasses import dataclass, field, replace
import operator
from collections.abc import Mapping

import numpy as np
from scipy.optimize import differential_evolution, dual_annealing, minimize

from ._ffi import CompiledContext, _EMPTY, kernel
from .io import read_geog, read_layer_matrices, read_newick, read_times
from .operators import _priors, _states, rcpp_areas_list_to_states_list, rcpp_states_list_to_DEmat
from .traits import TraitModel
from .stratified import compile_sections, config_single, effective_times, section_the_tree
from .uncertainty.contracts import CONTRACTS, RootPriorSpec, UnsupportedModelContractError


@dataclass(frozen=True)
class Parameter:
    type: str
    init: float
    min: float
    max: float


_BIN = {ast.Add:operator.add,ast.Sub:operator.sub,ast.Mult:operator.mul,
        ast.Div:operator.truediv,ast.Pow:operator.pow}


def _expression(node, get):
    if isinstance(node, ast.Name):
        return get(node.id)
    if isinstance(node, ast.Constant) and type(node.value) in (float,int):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):
        value=_expression(node.operand,get)
        return -value if isinstance(node.op,ast.USub) else value
    if isinstance(node, ast.BinOp) and type(node.op) in _BIN:
        return _BIN[type(node.op)](_expression(node.left,get),_expression(node.right,get))
    raise ValueError("parameter links allow only names, numbers and arithmetic")


@dataclass(frozen=True)
class BioGeoBEARSModelObject:
    name: str
    params: Mapping[str, Parameter]

    def __post_init__(self):
        linked = {name:ast.parse(p.type.replace("^","**"),mode="eval").body
                  for name,p in self.params.items() if p.type not in ("fixed","free")}
        object.__setattr__(self,"_linked",linked)
        for name,p in self.params.items():
            if not np.isfinite([p.init,p.min,p.max]).all() or p.min>p.max:
                raise ValueError(f"invalid bounds or initial value for {name}")
            if p.type=="free" and not p.min<=p.init<=p.max:
                raise ValueError(f"initial value for {name} outside bounds")
        self.resolve()

    @property
    def free(self):
        return tuple(name for name,p in self.params.items() if p.type=="free")

    def resolve(self, overrides=None):
        supplied={} if overrides is None else dict(overrides)
        if supplied.keys()-self.params.keys():
            raise ValueError(f"unknown parameters: {sorted(supplied.keys()-self.params.keys())}")
        active=set()
        values={}

        def value(name):
            if name not in self.params:
                raise ValueError(f"undefined parameter {name}")
            if name in active:
                raise ValueError(f"circular parameter link at {name}")
            if name in values:
                return values[name]
            active.add(name)
            p=self.params[name]
            x=_expression(self._linked[name],value) if name in self._linked else supplied.get(name,p.init)
            if not np.isfinite(x):
                raise ValueError(f"non-finite parameter {name}")
            values[name]=float(x)
            active.remove(name)
            return x

        for name in self.params:
            value(name)
        return values

    def with_parameters(self, **changes):
        params=dict(self.params)
        for name, change in changes.items():
            params[name] = change if isinstance(change, Parameter) else replace(params[name], **change)
        return BioGeoBEARSModelObject(self.name,params)


def model_template(name):
    display=name.upper().replace(" ","")
    parts=display.split("+")
    base=parts[0]
    if base not in ("DEC","DIVALIKE","BAYAREALIKE") or set(parts[1:])-{"J","X","W","N","U"} or len(set(parts[1:]))!=len(parts[1:]):
        raise ValueError(f"unknown model {name}")
    categories={"DEC":3,"DIVALIKE":2,"BAYAREALIKE":1}[base]
    p = {
        "d":Parameter("free",.01,1e-12,5), "e":Parameter("free",.01,1e-12,5),
        "a":Parameter("fixed",0,1e-12,5), "b":Parameter("fixed",1,1e-12,1),
        "x":Parameter("free" if "X" in parts else "fixed",0,-2.5,2.5),
        "n":Parameter("free" if "N" in parts else "fixed",0,-10,10),
        "w":Parameter("free" if "W" in parts else "fixed",1,-10,10),
        "u":Parameter("free" if "U" in parts else "fixed",0,-10,10),
        "j":Parameter("free" if "J" in parts else "fixed",.01 if "J" in parts else 0,1e-5,categories-1e-5),
        "ysv":Parameter(f"{categories}-j",categories,0,categories),
        "y":Parameter(f"ysv/{categories}",1,0,1),
        "s":Parameter("ysv/3" if base=="DEC" else "fixed",1 if base=="DEC" else 0,0,1),
        "v":Parameter(f"ysv/{categories}" if base!="BAYAREALIKE" else "fixed",1 if base!="BAYAREALIKE" else 0,0,1),
        "mx01":Parameter("fixed",.0001,.0001,.9999),
        "mx01s":Parameter("mx01",.0001,.0001,.9999),
        "mx01v":Parameter("fixed" if base=="DIVALIKE" else "mx01",.5 if base=="DIVALIKE" else .0001,.0001,.9999),
        "mx01j":Parameter("mx01",.0001,.0001,.9999),
        "mx01y":Parameter("fixed" if base=="BAYAREALIKE" else "mx01",.9999 if base=="BAYAREALIKE" else .0001,.0001,.9999),
    }
    normalized = display
    for mod in ("X","W","N","U"):
        normalized = normalized.replace(f"+{mod}", f"+{mod.lower()}")
    return BioGeoBEARSModelObject(normalized,p)


@dataclass
class ModelConfig:
    tree: object
    geog: object
    model: str | BioGeoBEARSModelObject = "DEC"
    max_range_size: int | None = None
    include_null_range: bool = True
    allow_null_tips: bool = False
    use_ambiguities: bool = False
    fixnode: object = None
    fixlikes: object = None
    timeperiods: object = None
    areas_allowed: object = None
    areas_adjacency: object = None
    dispersal_multipliers: object = None
    distances: object = None
    envdistances: object = None
    area_of_areas: object = None
    trait_values: object = None
    trait_model: object = None
    precision: str = "strict_f64"
    min_branchlength: float = 1e-6
    force_sparse: bool = False
    rootedge: bool = False
    likelihood_conditioning: str = "legacy_unconditioned"
    root_prior: RootPriorSpec | None = None


@dataclass
class EvaluationResult:
    lnL: float
    ancestral_probs_branch_top: np.ndarray | None
    ancestral_probs_branch_bottom: np.ndarray | None
    condlikes_table: np.ndarray | None
    downpass_likes_branch_bottom: np.ndarray | None
    states_list: tuple
    ranges_list: tuple
    tree: object = field(repr=False)
    downpass_likes_branch_top: np.ndarray | None = field(default=None,repr=False)
    downpass_table: np.ndarray | None = field(default=None,repr=False)
    likelihood_conditioning: str = "legacy_unconditioned"
    precision_mode: str = "strict_f64"
    root_prior_vec: np.ndarray | None = field(default=None,repr=False)
    log_tip_survival_prob: float | None = None

    def to_dataframe(self):
        import pandas as pd

        if self.ancestral_probs_branch_top is None:
            raise ValueError("this result contains only a likelihood")
        p=self.ancestral_probs_branch_top
        frame=pd.DataFrame(p,columns=self.ranges_list)
        frame.insert(0,"most_likely_state",[self.ranges_list[j] for j in np.argmax(p,axis=1)])
        frame.insert(0,"node_type",self.tree.node_table()["node.type"])
        frame.insert(0,"node",np.arange(len(p))+1)
        return frame


@dataclass
class FitResult:
    lnL: float
    params: dict[str, float]
    n_params: int
    aic: float
    aicc: float
    n_evals: int
    converged: bool
    evaluation: EvaluationResult
    model: object = field(repr=False)
    message: str = ""
    strict_proposal: dict | None = None

    @property
    def likelihood_conditioning(self):
        return self.evaluation.likelihood_conditioning

    @property
    def precision_mode(self):
        return self.evaluation.precision_mode

    @property
    def root_prior_vec(self):
        return self.evaluation.root_prior_vec

    @property
    def log_tip_survival_prob(self):
        return self.evaluation.log_tip_survival_prob

    def summary(self):
        from .stats import compare_models
        return compare_models([self])

    def run_bsm(self, num_maps=100, mode="uniformization", num_threads=1, seed=12345):
        from .bsm import run_bsm
        return run_bsm(self,num_maps,mode,num_threads,seed)

    def profile_ci(self, **options):
        from .uncertainty.profile import profile_ci
        return profile_ci(self, **options)

    def profile_surface(self, param_x, param_y, **options):
        from .uncertainty.profile import profile_surface
        return profile_surface(self, param_x, param_y, **options)

    def bootstrap_lrt(self, alt_fit, **options):
        from .uncertainty.bootstrap import bootstrap_lrt
        return bootstrap_lrt(self, alt_fit, **options)

    def run_bsm_adaptive(self, **options):
        from .uncertainty.diagnostics import run_bsm_adaptive
        return run_bsm_adaptive(self, **options)

    def plot(self, filepath, plot_what="pie", show_corners=True, title=None):
        from .plotting import plot_fit
        return plot_fit(self,filepath,plot_what,show_corners,title)


class BioGeoBEARS:
    def __init__(self, tree, geog=None, model="DEC", **options):
        self.config=tree if isinstance(tree,ModelConfig) else ModelConfig(tree,geog,model,**options)
        c=self.config
        if c.likelihood_conditioning not in CONTRACTS:
            raise ValueError("unknown likelihood_conditioning")
        self.likelihood_conditioning=c.likelihood_conditioning
        self.root_prior_spec=c.root_prior or RootPriorSpec(
            "legacy_flat_weights" if self.likelihood_conditioning=="legacy_unconditioned" else "uniform_nonnull")
        if not isinstance(self.root_prior_spec,RootPriorSpec):
            raise TypeError("root_prior must be a RootPriorSpec")
        if self.likelihood_conditioning=="unconditional_complete" and not c.allow_null_tips:
            raise ValueError("unconditional_complete requires allow_null_tips=True")
        if self.likelihood_conditioning=="tip_survival" and c.allow_null_tips:
            raise ValueError("tip_survival requires allow_null_tips=False")
        self.tree=read_newick(c.tree)
        self.geog=read_geog(c.geog)
        self.model=model_template(c.model) if isinstance(c.model,str) else c.model
        if c.precision not in ("strict_f64","legacy_f32"):
            raise ValueError("precision must be strict_f64 or legacy_f32")
        self.precision=c.precision
        self.areas=self.geog.areas
        self.na=len(self.areas)
        raw=rcpp_areas_list_to_states_list(self.areas,c.max_range_size,c.include_null_range)
        states=tuple(() if s==(-1,) else s for s in raw)
        self.null=int(c.include_null_range)
        self._contract_native=(self.likelihood_conditioning!="legacy_unconditioned" or c.root_prior is not None)
        self.raw_times=read_times(c.timeperiods)
        self.times=effective_times(self.raw_times,float(self.tree.ages[self.tree.root]))
        self.sections=section_the_tree(self.tree,self.times,exact_ages=self._contract_native) if len(self.times)>1 else None
        nlayer=max(1,len(self.times))

        def matrices(source, distance=False):
            return read_layer_matrices(source,self.na,nlayer,missing_diagonal=distance,
                                       areas=self.areas,raw_layers=len(self.raw_times))

        self.dist=matrices(c.distances,True)
        self.env=matrices(c.envdistances,True)
        self.manual=matrices(c.dispersal_multipliers)
        self.allowed_mats=matrices(c.areas_allowed)
        self.adjacency_mats=matrices(c.areas_adjacency)
        if c.area_of_areas is None:
            self.area_size=np.ones((nlayer,self.na))
        else:
            areas=np.asarray(c.area_of_areas,float)
            if len(self.raw_times)>nlayer and areas.shape==(len(self.raw_times),self.na):
                areas=areas[:nlayer]
            self.area_size=np.broadcast_to(areas,(nlayer,self.na)).copy()
            if np.any(self.area_size <= 0) or not np.isfinite(self.area_size).all():
                raise ValueError("area sizes must be positive and finite")
        self.masks=tuple(np.array([not s or
                                   (a[s[0],list(s)].sum()==len(s) and
                                    adj[np.ix_(s,s)].sum()==len(s)**2) for s in states])
                         for a,adj in zip(self.allowed_mats,self.adjacency_mats))
        single_mask = self.masks[0] if self.sections is None else np.ones(len(states),bool)
        self.geographic_states=tuple(s for s,keep in zip(states,single_mask) if keep)
        if self.sections is None:
            if len(self.geographic_states)<=self.null:
                raise ValueError("the spatial constraints exclude all occupied ranges")
            self.masks=(np.ones(len(self.geographic_states),bool),)
        self.states=self.geographic_states
        sep="|" if any(len(a)>1 for a in self.areas) else ""
        self.ranges=tuple(sep.join(self.areas[i] for i in s) or "_" for s in self.states)
        self.tips=self.geog.likelihoods(self.tree.tip_labels,self.states,c.allow_null_tips,c.use_ambiguities)
        self.trait=c.trait_model
        if (c.trait_values is None) != (self.trait is None) or self.trait is not None and not isinstance(self.trait,TraitModel):
            raise ValueError("trait_values requires a TraitModel")
        traits=1
        if self.trait:
            self.trait.arrays(self.model.resolve())
            obs=self.trait.observations(self.tree.tip_labels,c.trait_values)
            traits=len(self.trait.labels)
            self.tips=np.asfortranarray(np.concatenate([self.tips*obs[:,i,None] for i in range(traits)],axis=1))
            self.states=self.geographic_states*traits
            self.ranges=tuple(f"{r}|{t}" for t in self.trait.labels for r in self.ranges)
        def _checked_fixnode_id(raw_id):
            if isinstance(raw_id, (bool, np.bool_)):
                raise ValueError("fixnode IDs must be integers, not booleans")
            if isinstance(raw_id, (int, np.integer)):
                return int(raw_id) - 1
            if isinstance(raw_id, (float, np.floating)):
                if not np.isfinite(raw_id) or not float(raw_id).is_integer():
                    raise ValueError(f"non-integer fixnode ID {raw_id!r}")
                return int(raw_id) - 1
            raise ValueError(f"invalid fixnode ID {raw_id!r}")

        raw_fixes = () if c.fixnode is None else ([c.fixnode] if np.ndim(c.fixnode) == 0 else list(c.fixnode))
        self.fixes=tuple(_checked_fixnode_id(n) for n in raw_fixes)
        self.fixlikes=np.atleast_2d(np.asarray(c.fixlikes,float)) if self.fixes else _EMPTY
        if self.fixes and self.fixlikes.shape==(len(self.fixes),len(states)*traits) and not single_mask.all():
            self.fixlikes=self.fixlikes[:,np.tile(single_mask,traits)]
        if self.fixes and (self.fixlikes.shape!=(len(self.fixes),len(self.states)) or
                           any(n<self.tree.tip_count or n>=len(self.tree.labels) for n in self.fixes)
                           or len(set(self.fixes))!=len(self.fixes)):
            raise ValueError("fixnode must contain unique internal IDs, with one state vector per node")
        if self._contract_native:
            self._require_tier1a()
        root_allowed = np.tile(self.masks[-1], traits)
        self.root_likelihood_weights,self.root_sampling_probs,self.root_weight_log_scale=(
            self.root_prior_spec.resolve(self.states,self.ranges,self.likelihood_conditioning,
                                         allowed_mask=root_allowed))
        self.root_prior_vec=self.root_sampling_probs
        native_root = self.root_likelihood_weights if self._contract_native else _EMPTY
        if self.sections:
            self._geometry,self._prior_width=compile_sections(self.sections,self.geographic_states,self.masks,
                                                               self.null,self.na,self.fixes,self.fixlikes,c.min_branchlength,
                                                               len(self.trait.labels) if self.trait else 0,
                                                               root_prior=native_root,
                                                               exact_ages=self._contract_native)
            self._upper_m=[None]*nlayer
            last=None
            for i in reversed(range(nlayer)):
                self._upper_m[i]=last
                for piece in self.sections.pieces[i]:
                    if not piece.branch:
                        last=len(piece.lengths)
        else:
            self._prior_width=(max(map(len,self.states)),)
            edges=self.tree.edges
            self._geometry=config_single(self.tree.tip_count,len(self.states),edges,self.tree.lengths[edges[:,1]],
                                         self.tips,self.null,self.tree.root_edge if c.rootedge else 0,
                                         () if self.trait else self.states,() if self.trait else self.states[self.null:],
                                         self.na,self.fixes,self.fixlikes,c.min_branchlength,
                                         len(self.trait.labels) if self.trait else 0,
                                         prior=native_root)

    def _require_tier1a(self):
        """Compatibility name for the validated probability-contract model gate."""
        c=self.config
        unsupported=[]
        if self.trait is not None: unsupported.append("trait models")
        if self.fixes: unsupported.append("fixed internal observations")
        if c.rootedge: unsupported.append("root edges")
        if c.use_ambiguities or self.geog.ambiguous: unsupported.append("ambiguous tip observations")
        if (np.any(self.tree.lengths[self.tree.edges[:,1]]<c.min_branchlength) or
            np.any(self.tree.ages[:self.tree.tip_count]>self.tree.fossils_older_than) or
            (len(self.times)>1 and np.any(self.tree.ages[:self.tree.tip_count]>=self.times[0]))):
            unsupported.append("fossil-hook or non-contemporaneous fossil branches")
        if not self.null: unsupported.append("include_null_range=False")
        parts=self.model.name.upper().split("+")
        if parts[0] not in ("DEC","DIVALIKE","BAYAREALIKE") or set(parts[1:])-{"J","X","W","N","U"}:
            unsupported.append("nonstandard model family")
        p=self.model.resolve()
        if any(p[k]!=v or self.model.params[k].type=="free"
               for k,v in {"a":0,"b":1}.items()):
            unsupported.append("nonstandard rate modifiers")
        base_name = parts[0] + ("+J" if "J" in parts else "")
        if parts[0] in ("DEC","DIVALIKE","BAYAREALIKE"):
            base=model_template(base_name)
            if any(self.model.params[k]!=base.params[k] for k in base.params if k not in ("d","e","j","x","w","n","u")):
                unsupported.append("nonstandard cladogenesis parameters")
        for label,mats in (("areas_allowed",self.allowed_mats),("areas_adjacency",self.adjacency_mats)):
            arr=np.asarray(mats)
            if not np.all((arr==0.0) | (arr==1.0)):
                unsupported.append(f"non-binary {label} entries")
        for i in range(len(self.masks)-1):
            if np.any(self.masks[i+1] & ~self.masks[i]):
                unsupported.append("non-monotone state-space disappearance across time strata")
                break
        if unsupported:
            raise UnsupportedModelContractError("probability contracts do not support "+", ".join(unsupported))
        present_mask=self.masks[0]
        if np.any(self.tips[:,~present_mask]>0) or not np.all(self.tips[:,present_mask].sum(axis=1)>0):
            raise UnsupportedModelContractError("tip observations assign mass outside the present stratum's allowed state space")
        for key,matrix in (("x",self.dist),("n",self.env),("w",self.manual),("u",self.area_size)):
            if self.model.params[key].type=="free" and self.model.params[key].min<0 and np.any(np.asarray(matrix)<=0):
                raise ValueError(f"free {key} with negative bounds requires strictly positive modifier entries")

    def _validate_contract_params(self,params,*,force=False):
        if not force and self.likelihood_conditioning=="legacy_unconditioned" and not self._contract_native:
            return
        categories={"DEC":3,"DIVALIKE":2,"BAYAREALIKE":1}[self.model.name.upper().split("+")[0]]
        if params["d"]<0 or params["e"]<0 or not 0<=params["j"]<categories:
            raise ValueError(f"probability contracts require d,e >= 0 and 0 <= j < {categories}")
        modifiable=("d","e","j","x","w","n","u")
        expected=self.model.resolve({k:params[k] for k in modifiable})
        if any(params[k]!=expected[k] for k in params if k not in modifiable):
            raise UnsupportedModelContractError("parameter overrides cannot enable unvalidated rate or cladogenesis modifiers")
        with np.errstate(over="ignore",divide="ignore",invalid="ignore"):
            modifiers=(np.asarray(self.dist)**params["x"]
                       * np.asarray(self.env)**params["n"]
                       * np.asarray(self.manual)**params["w"])
            area_mod=np.asarray(self.area_size)**params["u"]
        if (not np.isfinite(modifiers).all() or np.any(modifiers<0) or
            not np.isfinite(area_mod).all() or np.any(area_mod<0)):
            raise ValueError("spatial/environmental modifier exponents produced invalid rate multipliers")

    def prepared_likelihood(self,params=None,*,num_threads=1):
        """Compile once for repeated tip patterns at fixed model parameters."""
        self._require_tier1a()
        return PreparedLikelihood(self,params,num_threads=num_threads)

    def _layer_arrays(self, i, p, *, upper=False):
        env=self.env[0] if upper and self.sections else self.env[i]
        dist=self.dist[i]**p["x"] * env**p["n"] * self.manual[i]**p["w"]
        prior=tuple(p[f"mx01{kind}"] for kind in "svjy")
        weights=[p[k] for k in "svjy"]
        mult=self._upper_m[i] if upper and not self.config.force_sparse else None
        return (dist*p["d"], self.area_size[i]**p["u"]*p["e"], dist*p["a"],
                weights,dist if mult is None else mult*dist,
                *_priors(self.na,self._prior_width[i],prior,self.precision))

    def _trait_layer(self, i, params, upper=False):
        d,e,a,weights,dispersal,*priors=self._layer_arrays(i,params,upper=upper)
        m,rates,jump=self.trait.arrays(params)
        if self.sections and upper:
            jump=None
        local=tuple(s for s,keep in zip(self.geographic_states,self.masks[i]) if keep)
        n,nt=len(local),len(m)
        q=np.zeros((n*nt,n*nt),order="F")
        parts=[]
        for j,multiplier in enumerate(m):
            slot=slice(j*n,(j+1)*n)
            q[slot,slot]=rcpp_states_list_to_DEmat(self.areas,local,multiplier*d,e,multiplier*a,
                                                  precision=self.precision)
            for k in range(nt):
                if j!=k:
                    target=slice(k*n,(k+1)*n)
                    q[slot,target]=np.eye(n)*rates[j,k]
            part=kernel(4,*_states(local[self.null:]),weights,multiplier*dispersal,*priors,self.na,[],
                        precision=self.precision)
            parts.append(np.column_stack(part))
        np.fill_diagonal(q,0)
        np.fill_diagonal(q,-q.sum(axis=1))
        primary=[]
        extra=[]
        for t,part in enumerate(parts):
            row=part.copy()
            row[:,:3]+=t*n
            if jump is not None:
                split=part[:,:3].astype(int)
                st=local[self.null:]
                j_events=np.array([not set(st[l]+st[r]).issubset(st[p]) for p,l,r in split])
                row[j_events,3]*=jump[t,t]
                for k in range(nt):
                    if t==k or jump[t,k]==0:
                        continue
                    changed=row[j_events].copy()
                    changed[:,1:3]+=(k-t)*n
                    changed[:,3]=part[j_events,3]*jump[t,k]
                    extra.append(changed)
            primary.append(row)
        columns=np.concatenate([*primary,*extra],axis=0)
        denom=kernel(6,columns[:,0],columns[:,3],n*nt-self.null)[0]
        return q,columns,denom,dispersal

    def _layer_buffers(self, i, p, ancestors, full, upper=False):
        sparse=self.config.force_sparse
        if self.trait:
            q,columns,denom,_=self._trait_layer(i,p,upper)
            flags=[0,sparse,ancestors,len(self.trait.labels) if self.sections else full]
            return [flags,q,_EMPTY,_EMPTY,*columns.T,denom,_EMPTY,[0,0,0,0],*([_EMPTY]*11)]
        d,e,a,weights,jump,*priors=self._layer_arrays(i,p,upper=upper)
        return [[0,sparse,ancestors,full],*([_EMPTY]*9),[1,1,sparse,1e-26],d,e,a,
                weights,jump,*priors,self.na]

    def _context(self,num_threads):
        return CompiledContext(self._geometry,stratified=self.sections is not None,num_threads=num_threads)

    def _evaluate(self, ctx, params, full, precision=None, tip_likelihoods=None):
        self._validate_contract_params(params)
        if self.sections:
            args=[[full,full,1,1,1],self.tips if tip_likelihoods is None else tip_likelihoods]
            use_contract=self._contract_native
            for i in range(len(self.times)):
                separate=(not use_contract) and full and (
                    self.trait is not None or not np.array_equal(self.env[i],self.env[0]) or
                    params["j"]!=0 and self._upper_m[i] is not None and not self.config.force_sparse)
                args.extend(([separate],*self._layer_buffers(i,params,False,False)))
                if separate:
                    args.extend(self._layer_buffers(i,params,False,False,True))
        else:
            args=self._layer_buffers(0,params,full,full)
            if params["b"] != 1:
                args.append([params["b"]])
        active_precision=self.precision if precision is None else precision
        if self._contract_native or tip_likelihoods is not None:
            out=ctx.evaluate(args,active_precision,contract=self.likelihood_conditioning,
                             tip_likelihoods=tip_likelihoods)
            survival=float(out[0][-1]) if self.likelihood_conditioning=="tip_survival" else None
        else:
            out=ctx.evaluate(args,active_precision)
            survival=None
        ctx.last_log_tip_survival_prob=survival
        lnL=float(out[0][0])
        if not full:
            return lnL
        n=len(self.states)
        def matrix(x):
            return x.reshape((-1,n),order="F")
        if self.sections:
            cond,down_table,bottom,top,down,td=map(matrix,(out[2],out[3],out[6],out[7],out[8],out[9]))
        else:
            td,cond,down,up_bottom,up_top=map(matrix,(out[2],out[3],out[5],out[6],out[7]))
            with np.errstate(invalid="ignore",divide="ignore"):
                bottom=down*up_bottom
                top=td[:len(up_top)]*up_top
                top[self.tree.root]=td[self.tree.root]
                bottom=bottom/bottom.sum(axis=1,keepdims=True)
                top=top/top.sum(axis=1,keepdims=True)
            down_table=None
        return EvaluationResult(lnL,top,bottom,cond,down,self.states,self.ranges,self.tree,td,down_table,
            self.likelihood_conditioning,active_precision,self.root_sampling_probs.copy(),survival)

    def evaluate(self, params=None, *, num_threads=1, calc_ancprobs=True):
        p=self.model.resolve(params)
        with self._context(num_threads) as ctx:
            value=self._evaluate(ctx,p,calc_ancprobs)
            survival=getattr(ctx,"last_log_tip_survival_prob",None)
        if calc_ancprobs:
            return value
        return EvaluationResult(value,None,None,None,None,self.states,self.ranges,self.tree,
            likelihood_conditioning=self.likelihood_conditioning,precision_mode=self.precision,
            root_prior_vec=self.root_sampling_probs.copy(),log_tip_survival_prob=survival)

    def fit(self, optimizer="L-BFGS-B", num_threads=1, *, initial=None, seed=1, maxiter=500):
        free=self.model.free
        start=self.model.resolve(initial)
        bounds=[(self.model.params[k].min,self.model.params[k].max) for k in free]
        x0=np.asarray([start[k] for k in free])
        if any(not lo<=x<=hi for x,(lo,hi) in zip(x0,bounds)):
            raise ValueError("initial free parameters outside bounds")
        from ._optimization import optimize_likelihood
        with self._context(num_threads) as ctx:
            def score(theta,precision):
                params=self.model.resolve({**start,**dict(zip(free,theta))})
                return self._evaluate(ctx,params,False,precision)
            search=optimize_likelihood(score,list(free),x0,bounds,
                self.tree.ages[self.tree.root],self.precision,optimizer,maxiter,seed)
            if optimizer=="L-BFGS-B":
                from ._refinement import refine_completed_search
                search=refine_completed_search(score,list(free),x0,bounds,self.precision,maxiter,search)
            p=self.model.resolve({**start,**dict(zip(free,search.theta))})
            strict_proposal=None
            if search.proposal is not None:
                strict_proposal=dict(search.proposal)
                theta=strict_proposal.pop("theta")
                strict_proposal["params"]=self.model.resolve({**start,**dict(zip(free,theta))})
            evaluation=self._evaluate(ctx,p,True)
        k=sum(1 for name in free if self.model.params[name].min < self.model.params[name].max)
        n=self.tree.tip_count
        aic=2*k-2*evaluation.lnL
        aicc=aic+2*k*(k+1)/(n-k-1) if n>k+1 else float("inf")
        return FitResult(evaluation.lnL,p,k,aic,aicc,search.evaluations,
                         search.converged,evaluation,self,search.message,strict_proposal)


class PreparedLikelihood:
    """Owned context for fixed-parameter enumeration; tip order follows the tree."""

    def __init__(self,model,params=None,*,num_threads=1):
        self.model=model
        self.params=model.model.resolve(params)
        model._validate_contract_params(self.params)
        self.context=model._context(num_threads)

    def evaluate(self,tip_likelihoods=None,*,calc_ancprobs=False):
        tips=self.model.tips if tip_likelihoods is None else np.asarray(tip_likelihoods,dtype=float)
        if tips.shape!=self.model.tips.shape or not np.isfinite(tips).all() or np.any(tips<0):
            raise ValueError("tip likelihoods must be a finite nonnegative n_tips by n_states matrix")
        return self.model._evaluate(self.context,self.params,calc_ancprobs,tip_likelihoods=tips)

    def evaluate_tip_states(self,tip_states):
        indices=np.asarray(tip_states)
        if indices.ndim!=2 or indices.shape[1]!=self.model.tree.tip_count or indices.dtype.kind not in "iu":
            raise ValueError("tip_states must be an integer n_patterns by n_tips matrix")
        if np.any(indices<0) or np.any(indices>=len(self.model.states)):
            raise ValueError("tip state index outside the model state space")
        out=np.empty(len(indices))
        tips=np.zeros_like(self.model.tips,order="F")
        rows=np.arange(self.model.tree.tip_count)
        for i,pattern in enumerate(indices):
            tips.fill(0)
            tips[rows,pattern]=1
            out[i]=self.evaluate(tips)
        return out

    def close(self):
        self.context.close()

    def __enter__(self):
        self.context.__enter__()
        return self

    def __exit__(self,*exc):
        self.close()
