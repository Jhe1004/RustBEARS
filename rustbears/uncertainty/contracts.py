"""Explicit probability contracts and a shared, labelled root specification."""
from dataclasses import dataclass
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np


CONTRACTS = {
    "legacy_unconditioned": 0,
    "tip_survival": 1,
    "unconditional_complete": 2,
}


class ContractMismatchError(ValueError):
    """The requested operation would mix incompatible probability contracts."""


class UnsupportedModelContractError(ValueError):
    """This model and probability-contract combination has not been validated."""


@dataclass(frozen=True)
class RootPriorSpec:
    """Root weights aligned to model states, never inferred from observed tips.

    Mapping keys are integer area tuples or exact range labels (``"_"`` for
    the empty state). Sequence weights follow ``model.states``. Under the
    legacy contract the supplied weight scale is retained by the likelihood;
    simulation always uses a separately normalized probability vector.
    """

    mode: str = "uniform_nonnull"
    custom_weights: object = None

    def __post_init__(self):
        if self.mode not in ("uniform_nonnull", "custom_weights", "legacy_flat_weights"):
            raise ValueError("unknown root prior mode")
        if (self.mode == "custom_weights") != (self.custom_weights is not None):
            raise ValueError("custom_weights must be supplied only for custom_weights mode")
        if self.mode == "custom_weights":
            # Freeze a private snapshot. A caller editing their input dict or
            # array must not silently change the contract's provenance.
            if isinstance(self.custom_weights, Mapping):
                frozen=MappingProxyType({key:float(value) for key,value in self.custom_weights.items()})
            else:
                weights=np.asarray(self.custom_weights,dtype=float)
                if weights.ndim!=1:
                    raise ValueError("custom root weight sequence must be one-dimensional")
                frozen=tuple(float(x) for x in weights)
            object.__setattr__(self,"custom_weights",frozen)

    def resolve(self, states, ranges, contract, allowed_mask=None):
        if contract not in CONTRACTS:
            raise ValueError("unknown likelihood_conditioning")
        occupied = np.array([bool(s) for s in states])
        allowed = None if allowed_mask is None else np.asarray(allowed_mask, dtype=bool)
        if allowed is not None and allowed.shape != (len(states),):
            raise ValueError("root allowed_mask must have one entry per model state")
        if self.mode == "custom_weights":
            if isinstance(self.custom_weights, Mapping):
                index = {tuple(s): i for i, s in enumerate(states)}
                labels = {s: i for i, s in enumerate(ranges)}
                weights = np.zeros(len(states), dtype=float)
                seen = set()
                for key, value in self.custom_weights.items():
                    if isinstance(key, str):
                        if key not in labels:
                            raise ValueError(f"unknown root state label {key!r}")
                        i = labels[key]
                    elif isinstance(key, tuple) and all(
                        isinstance(x, (int, np.integer)) and not isinstance(x, (bool, np.bool_))
                        for x in key
                    ):
                        if len(set(key)) != len(key) or tuple(sorted(key)) not in index:
                            raise ValueError(f"unknown root state {key!r}")
                        i = index[tuple(sorted(key))]
                    else:
                        raise ValueError("root state keys must be range labels or integer area tuples")
                    if i in seen:
                        raise ValueError("the same root state was specified more than once")
                    seen.add(i)
                    weights[i] = value
            else:
                weights = np.asarray(self.custom_weights, dtype=float).copy()
                if weights.shape != (len(states),):
                    raise ValueError("root weight vector must have one entry per model state")
        else:
            weights = occupied.astype(float)
        if not np.isfinite(weights).all() or np.any(weights < 0):
            raise ValueError("root weights must be finite and nonnegative")
        if allowed is not None:
            if self.mode == "custom_weights" and np.any(weights[~allowed] > 0):
                raise ValueError("custom root weights assign positive mass to states disallowed in the root stratum")
            weights[~allowed] = 0.0
        if contract == "tip_survival":
            weights[~occupied] = 0.0
        total = float(weights.sum())
        if not np.isfinite(total) or total <= 0:
            raise ValueError("root weights must have positive finite total mass")
        sampling = weights / total
        if contract == "legacy_unconditioned" and self.mode != "uniform_nonnull":
            likelihood = weights
            log_scale = float(np.log(total))
        else:
            likelihood = sampling.copy()
            log_scale = 0.0
        for array in (likelihood, sampling):
            array.setflags(write=False)
        return likelihood, sampling, log_scale

    def __reduce__(self):
        # MappingProxyType freezes inputs but cannot itself be pickled. Tasks
        # serialize a plain private copy and rebuild the frozen specification.
        weights = dict(self.custom_weights) if isinstance(self.custom_weights, Mapping) else self.custom_weights
        return (type(self), (self.mode, weights))


def resolve_simulation_conditioning(contract, requested):
    if requested == "from_model":
        if contract == "legacy_unconditioned":
            raise ContractMismatchError(
                "legacy_unconditioned has no unique generating distribution; explicitly choose "
                "conditioning='tip_survival' or 'unconditional_complete'"
            )
        return contract
    if requested not in ("tip_survival", "unconditional_complete"):
        raise UnsupportedModelContractError("probability contracts support tip_survival and unconditional_complete simulation")
    return requested
