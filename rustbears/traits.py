"""Discrete trait-by-geography state models."""
from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np


@dataclass(frozen=True)
class TraitModel:
    labels: tuple
    dispersal_multipliers: object
    transition_rates: object
    jump_transition: object = None

    def arrays(self, params):
        n=len(self.labels)
        if n<2 or len(set(self.labels))!=n:
            raise ValueError("a trait model needs at least two unique labels")

        def values(x,shape):
            x=np.asarray(x,dtype=object)
            if x.shape!=shape:
                raise ValueError(f"expected trait model dimensions {shape}")
            result=np.empty(shape)
            for i in np.ndindex(shape):
                v=x[i]
                result[i]=params[v] if isinstance(v,str) and v in params else float(v)
            if not np.isfinite(result).all() or (result<0).any():
                raise ValueError("trait rates and multipliers must be finite and nonnegative")
            return result

        m=values(self.dispersal_multipliers,(n,))
        t=values(self.transition_rates,(n,n))
        np.fill_diagonal(t,0)
        j=None if self.jump_transition is None else values(self.jump_transition,(n,n))
        if j is not None and not np.allclose(j.sum(axis=1),1,atol=1e-12,rtol=0):
            raise ValueError("trait transition probabilities at cladogenetic jumps must sum to one")
        return m,t,j

    def observations(self, labels, source):
        if not isinstance(source, Mapping):
            from .io import read_geog
            r=read_geog(source)
            if tuple(r.areas)!=tuple(map(str,self.labels)):
                raise ValueError("trait file states differ from the trait model")
            source={name:tuple(self.labels[k] for k in ix) for name,ix in r.ranges.items()}
        if set(labels)!=set(source):
            raise ValueError("trait and tree taxon sets differ")
        matrix=np.zeros((len(labels),len(self.labels)))
        index={label:i for i,label in enumerate(self.labels)}
        for i,tip in enumerate(labels):
            record=source[tip]
            if record is None:
                matrix[i,:]=1
            elif isinstance(record,(tuple,list,set)):
                matrix[i,[index[k] for k in record]]=1
            else:
                matrix[i,index[record]]=1
        if not np.all(matrix.sum(axis=1)>0):
            raise ValueError("trait observations cannot have zero support")
        return matrix
