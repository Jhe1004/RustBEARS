"""Shared guards for comparisons of fitted probability models.

Contract conversion always rebuilds the models from their in-memory observed
tree and geography and refits both models. It never relabels old likelihoods.
"""
from dataclasses import dataclass, replace
import copy
from collections.abc import Mapping

import numpy as np

from .contracts import CONTRACTS, ContractMismatchError


def clone_model(model, *, geog=None, model_object=None, **changes):
    """Create an independent, serializable model with the same in-memory data."""
    from ..model import BioGeoBEARS
    config=copy.deepcopy(replace(model.config,tree=model.tree,geog=model.geog,
                                model=model.model,root_prior=model.root_prior_spec))
    # Inputs may originally have been filenames. A refit uses the matrices
    # that produced this fit, even if a file was edited afterwards.
    for key,values in (("distances",model.dist),("envdistances",model.env),
                       ("dispersal_multipliers",model.manual),("area_of_areas",model.area_size),
                       ("areas_allowed",model.allowed_mats),("areas_adjacency",model.adjacency_mats)):
        if getattr(model.config,key) is not None:
            setattr(config,key,np.array(values,copy=True))
    if model.config.timeperiods is not None:
        config.timeperiods=np.array(model.raw_times,copy=True)
    if geog is not None:
        config.geog=copy.deepcopy(geog)
    if model_object is not None:
        config.model=model_object
    for key,value in changes.items():
        if key not in config.__dataclass_fields__:
            raise ValueError(f"unknown model configuration {key}")
        setattr(config,key,value)
    return BioGeoBEARS(config)


def model_manifest(model):
    """Fingerprint the resolved input, including files already read into memory."""
    return dict(tree=dict(labels=model.tree.labels,children=model.tree.children,
                          lengths=model.tree.lengths,ages=model.tree.ages,root_edge=model.tree.root_edge),
                areas=model.areas,states=model.states,tips=model.tips,
                parameters={key:dict(type=p.type,init=p.init,min=p.min,max=p.max)
                            for key,p in model.model.params.items()},model_name=model.model.name,
                root_likelihood_weights=model.root_likelihood_weights,
                root_sampling_probs=model.root_sampling_probs,root_prior_mode=model.root_prior_spec.mode,
                precision=model.precision,contract=model.likelihood_conditioning,
                dist=model.dist,env=model.env,manual=model.manual,area_size=model.area_size,
                allowed_mats=model.allowed_mats,adjacency_mats=model.adjacency_mats,
                times=model.times,masks=model.masks,fixes=model.fixes,fixlikes=model.fixlikes,
                include_null_range=model.config.include_null_range,
                allow_null_tips=model.config.allow_null_tips,min_branchlength=model.config.min_branchlength,
                force_sparse=model.config.force_sparse,rootedge=model.config.rootedge)


def _tree_signature(tree):
    """Compare rooted topology/branch lengths independently of child rotation."""
    clades={i:(str(label),) for i,label in enumerate(tree.tip_labels)}
    stack=[(tree.root,False)]
    while stack:
        node,visited=stack.pop()
        children=tree.children[node]
        if not children:
            continue
        if visited:
            clades[node]=tuple(sorted(t for child in children for t in clades[child]))
        else:
            stack.append((node,True))
            stack.extend((child,False) for child in children)
    branches=tuple(sorted((clades[node],float(tree.lengths[node]),float(tree.ages[node]))
                          for node in clades if node!=tree.root))
    return (branches,tuple(sorted(tree.tip_labels)),float(tree.ages[tree.root]),float(tree.root_edge))


def _observation_signature(model):
    tip_rows={str(label):tuple(float(v) for v in model.tips[i])
              for i,label in enumerate(model.tree.tip_labels)}
    return (tuple(model.areas),tuple(tuple(int(x) for x in state) for state in model.states),
            tuple(sorted(tip_rows.items())),bool(model.config.allow_null_tips),
            bool(model.config.use_ambiguities))


def validate_compatible_fits(fits,*,require_converged=True):
    """Validate the observed data, root measure, precision and score convention.

    Candidate model parameters may differ. State and area order must agree;
    a rotated but otherwise identical rooted tree is accepted. Callers that
    combine node arrays must still align those arrays by clade.
    """
    from ..model import FitResult
    fits=tuple(fits)
    if not fits:
        raise ValueError("at least one fit is required")
    for fit in fits:
        if not isinstance(fit,FitResult):
            raise TypeError("all inputs must be FitResult objects")
        if require_converged and not fit.converged:
            raise ValueError("comparison requires converged observed fits")
        if not np.isfinite(fit.lnL) or not np.isfinite(list(fit.params.values())).all():
            raise ValueError("comparison requires finite fitted scores and parameters")
        model=fit.model
        expected=model.model.resolve({key:fit.params[key] for key in model.model.free})
        if any(fit.params.get(key)!=value for key,value in expected.items()):
            raise ContractMismatchError("fit fixed/linked parameters do not match its model definition")
        if any(not model.model.params[key].min<=fit.params[key]<=model.model.params[key].max
               for key in model.model.free):
            raise ContractMismatchError("fitted free parameters fall outside their model bounds")
        if (fit.likelihood_conditioning!=model.likelihood_conditioning or
            fit.precision_mode!=model.precision or
            not np.array_equal(fit.root_prior_vec,model.root_sampling_probs)):
            raise ContractMismatchError("fit metadata no longer matches its model")
    first=fits[0]
    base=first.model
    tree=_tree_signature(base.tree)
    observed=_observation_signature(base)
    for fit in fits[1:]:
        model=fit.model
        checks=(
            (fit.likelihood_conditioning==first.likelihood_conditioning,"likelihood contracts"),
            (fit.precision_mode==first.precision_mode,"numerical precision modes"),
            (_tree_signature(model.tree)==tree,"observed trees or branch lengths"),
            (_observation_signature(model)==observed,"observed geography or state spaces"),
            (np.array_equal(model.root_likelihood_weights,base.root_likelihood_weights),"root likelihood weights"),
            (np.array_equal(model.root_sampling_probs,base.root_sampling_probs),"root sampling probabilities"),
        )
        for valid,label in checks:
            if not valid:
                raise ContractMismatchError(f"incompatible {label}")
    return fits


@dataclass(frozen=True)
class PreparedFitComparison:
    null_fit: object
    alt_fit: object
    recomputed_observed_under_new_contract: bool
    original_conditioning: str
    likelihood_conditioning: str

    @property
    def fits(self):
        return (self.null_fit,self.alt_fit)


def prepare_comparison_fits(null_fit,alt_fit,*,target_conditioning="match_fits",
                            refit_observed=False,fit_options=None):
    """Prepare two compatible fits, with explicit opt-in for contract changes.

    The same root specification and model parameter definitions are retained.
    A/B normalize the root weights according to their contracts; switching to
    C preserves an explicitly normalized prior if that was the source prior.
    Original fits and their input models are never mutated.
    """
    from ..model import BioGeoBEARS
    original=validate_compatible_fits((null_fit,alt_fit))
    if not isinstance(refit_observed,(bool,np.bool_)):
        raise TypeError("refit_observed must be a boolean")
    if target_conditioning not in ("match_fits",*CONTRACTS):
        raise ValueError("unknown target_conditioning")
    source=original[0].likelihood_conditioning
    target=source if target_conditioning=="match_fits" else target_conditioning
    if target==source:
        return PreparedFitComparison(*original,False,source,target)
    if not refit_observed:
        raise ContractMismatchError(
            "changing likelihood_conditioning requires refit_observed=True; old scores cannot be relabelled"
        )
    if fit_options is not None and not isinstance(fit_options,Mapping):
        raise TypeError("fit_options must be a mapping")
    # Build both converted models before fitting either, so an impossible
    # observation/contract conversion fails without a partially completed pair.
    models=[]
    for fit in original:
        old=fit.model
        models.append(clone_model(old,likelihood_conditioning=target,
                                  allow_null_tips=(target=="unconditional_complete" or
                                                  target=="legacy_unconditioned" and old.config.allow_null_tips)))
    converted=[]
    for old,model in zip(original,models):
        options=dict(fit_options or {})
        options.setdefault("initial",dict(old.params))
        converted.append(model.fit(**options))
    converted=validate_compatible_fits(converted)
    return PreparedFitComparison(*converted,True,source,target)
