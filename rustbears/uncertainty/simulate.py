"""Native fixed-tree simulation under an explicitly selected distribution."""
from dataclasses import dataclass

import numpy as np

from .._ffi import PreparedSampler
from ..io import TipRanges
from ..operators import _states
from .contracts import resolve_simulation_conditioning, UnsupportedModelContractError


@dataclass
class SimulatedBatch:
    tip_states: np.ndarray
    node_states: np.ndarray
    resolved_conditioning: str
    root_distribution: str
    root_prior_vec: np.ndarray
    precision_mode: str
    areas: tuple
    states: tuple
    tip_labels: tuple
    start_rep: int = 1
    seed: int = 12345

    @property
    def conditioning(self):
        return self.resolved_conditioning

    def to_tip_ranges(self,rep):
        if isinstance(rep,(bool,np.bool_)) or not isinstance(rep,(int,np.integer)) or not 0<=rep<len(self.tip_states):
            raise IndexError("rep must be a zero-based simulated replicate index")
        return TipRanges(self.areas,{label:self.states[int(state)]
            for label,state in zip(self.tip_labels,self.tip_states[rep])})


def sampler_buffers(model,params):
    from ..bsm import _layers,_steps
    steps,layers,hooks,fixes,segments=_steps(model,params)
    return [[model.tree.tip_count,model.tree.root],*_states(model.states),model.tips,
            steps,model.tree.ages,layers,hooks,fixes,segments,*_layers(model,params)]


class PreparedSimulation:
    """A reusable fixed-parameter sampler; batch IDs determine random streams."""

    def __init__(self,fit_or_model,params=None,*,conditioning="from_model",root_distribution="model_prior"):
        from ..model import BioGeoBEARS,FitResult
        if isinstance(fit_or_model,FitResult):
            self.model=fit_or_model.model
            supplied={**fit_or_model.params,**({} if params is None else params)}
        elif isinstance(fit_or_model,BioGeoBEARS):
            self.model=fit_or_model
            supplied=params
        else:
            raise TypeError("fit_or_model must be a BioGeoBEARS model or FitResult")
        m=self.model
        m._require_tier1a()
        self.params=m.model.resolve(supplied)
        m._validate_contract_params(self.params,force=True)
        self.conditioning=resolve_simulation_conditioning(m.likelihood_conditioning,conditioning)
        if root_distribution not in ("model_prior","empirical_posterior"):
            raise ValueError("unknown root_distribution")
        self.root_distribution=root_distribution
        prior=(m.evaluate(self.params).ancestral_probs_branch_top[m.tree.root].copy()
               if root_distribution=="empirical_posterior" else m.root_sampling_probs.copy())
        if self.conditioning=="tip_survival":
            prior[np.array([not s for s in m.states])]=0
        mass=float(prior.sum())
        if not np.isfinite(prior).all() or np.any(prior<0) or mass<=0:
            raise ValueError("simulation root distribution has no admissible finite positive mass")
        self.root_prior_vec=prior/mass
        mode=(4 if root_distribution=="empirical_posterior" else 1) if self.conditioning=="tip_survival" else 2
        self.sampler=PreparedSampler(sampler_buffers(m,self.params),self.root_prior_vec,mode)

    def draw(self,n_reps=100,*,start_rep=1,seed=12345,num_threads=1):
        out=self.sampler.draw(start_rep,n_reps,seed,num_threads)
        m=self.model
        tips=out[0].reshape((-1,m.tree.tip_count)).astype(np.int64)
        nodes=out[1].reshape((-1,len(m.tree.labels))).astype(np.int64)
        return SimulatedBatch(tips,nodes,self.conditioning,self.root_distribution,self.root_prior_vec.copy(),
            m.precision,tuple(m.areas),tuple(m.states),tuple(m.tree.tip_labels),int(start_rep),int(seed))

    def close(self):
        self.sampler.close()

    def __enter__(self):
        self.sampler.__enter__()
        return self

    def __exit__(self,*_):
        self.close()


def simulate_datasets(fit_or_model,params=None,n_reps=100,conditioning="from_model",
                      root_distribution="model_prior",seed=12345,num_threads=1):
    """Return simulated tip and internal states; legacy scores require a distribution choice."""
    with PreparedSimulation(fit_or_model,params,conditioning=conditioning,
                            root_distribution=root_distribution) as sampler:
        return sampler.draw(n_reps,seed=seed,num_threads=num_threads)
