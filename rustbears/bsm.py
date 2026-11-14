"""Whole-tree biogeographic stochastic mapping and structured event tables."""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ._ffi import kernel, PreparedSampler
from .operators import _states, rcpp_states_list_to_DEmat

_EVENT_TYPES = "deajysvt"
_CLADO = ("map","node","left_node","right_node","ancestral_state","left_state","right_state",
          "time_bp","stratum","event_type","source_area","destination_area")
_ANA = ("map","node","stratum","from_state","to_state","time_bp","segment_time",
        "event_type","source_area","destination_area")


@dataclass
class BSMResult:
    clado_events_df: pd.DataFrame
    ana_events_df: pd.DataFrame
    summary_counts: dict
    num_maps: int

    def _tables(self, table):
        groups=dict(tuple(table.groupby("map",sort=False)))
        return [groups.get(i,table.iloc[0:0]).reset_index(drop=True) for i in range(1,self.num_maps+1)]

    @property
    def clado_events_tables(self):
        return self._tables(self.clado_events_df)

    @property
    def ana_events_tables(self):
        return self._tables(self.ana_events_df)

    def plot(self, filepath):
        from .plotting import plot_bsm
        return plot_bsm(self, filepath)


def _steps(model, params):
    t=model.tree
    pair=t.edges.reshape((-1,2,2))
    steps=np.column_stack((pair[:,0,0],pair[:,0,1],pair[:,1,1]))
    lengths=t.lengths if model.sections else t.lengths**params["b"]
    layer_at_node=np.searchsorted(model.times,t.ages,side="right")
    hooks=np.zeros(len(t.labels),dtype=bool)
    hooks[steps[:,0]]=(lengths[steps[:,1:]]<model.config.min_branchlength).any(axis=1)
    segments=[]
    parent=t.parent
    for child,node in enumerate(parent):
        if node<0:
            continue
        old,young=t.ages[node],t.ages[child]
        cuts=[old,*reversed([x for x in model.times if young<x<old]),young]
        for start,end in zip(cuts,cuts[1:]):
            lid=int(np.searchsorted(model.times,(start+end)*.5,side="right"))
            duration=start-end if model.sections else lengths[child]
            if lengths[child]<model.config.min_branchlength:
                duration=0
            segments.append((child,lid,duration,start))
    fixes=np.ones((len(t.labels),len(model.states)),order="F")
    if model.fixes:
        fixes[list(model.fixes),:]=model.fixlikes
    return steps.ravel(), layer_at_node, hooks, fixes, np.asarray(segments).ravel()


def _layers(model, params):
    n=len(model.states)
    for i,mask in enumerate(model.masks):
        if model.trait:
            local,events,_,dispersal=model._trait_layer(i,params)
            mask=np.tile(mask,len(model.trait.labels))
            columns=np.flatnonzero(mask)
            q=np.zeros((n,n),order="F")
            q[np.ix_(columns,columns)]=local
            events[:,:3]=columns[events[:,:3].astype(int)+model.null]
            yield q
            yield mask
            yield events.ravel()
            yield dispersal
            continue
        d,e,a,weights,dispersal,*priors=model._layer_arrays(i,params)
        columns=np.flatnonzero(mask)
        states=[model.states[s] for s in columns]
        q_local=rcpp_states_list_to_DEmat(model.areas,states,d,e,a,precision=model.precision)
        q=np.zeros((n,n),order="F")
        q[np.ix_(columns,columns)]=q_local
        clado_states=states[model.null:]
        idx=columns[model.null:]
        c=kernel(4,*_states(clado_states),weights,dispersal,*priors,model.na,[],
                 precision=model.precision)
        events=np.column_stack([*(idx[x.astype(int)] for x in c[:3]), c[3]])
        yield q
        yield mask
        yield events.ravel()
        yield dispersal


def _frame(raw, cols, model):
    frame=pd.DataFrame(raw.reshape((-1,len(cols))),columns=cols)
    for key in cols:
        if key not in ("time_bp","segment_time"):
            frame[key]=frame[key].astype(np.int64)
    frame["event_type"]=np.asarray(list(_EVENT_TYPES),dtype="U1")[frame["event_type"]]
    for k in ("node","left_node","right_node","stratum"):
        if k in frame:
            frame[k]+=1
    for k in ("source_area","destination_area"):
        # -1 is the absence of a source or destination.
        ids=frame[k].to_numpy()
        frame[k]=np.array((None,*model.areas),dtype=object)[ids+1]
    for k in cols:
        if k.endswith("_state"):
            frame[k.removesuffix("_state")+"_range"]=np.asarray(model.ranges)[frame[k]]
    return frame


def count_ana_clado_events(clado, ana=None, areas=None, num_maps=None):
    if isinstance(clado,BSMResult):
        return clado.summary_counts
    if areas is None or num_maps is None or num_maps<1:
        raise ValueError("provide areas and the total number of maps")
    names=pd.Index(areas)
    if not names.is_unique:
        raise ValueError("area names must be unique")

    def counts(table):
        idx=names.get_indexer(table["source_area"])
        jdx=names.get_indexer(table["destination_area"])
        keep=(idx>=0)&(jdx>=0)
        x=np.zeros((num_maps,len(names),len(names)),dtype=np.int64)
        np.add.at(x,(table["map"].to_numpy()[keep]-1,idx[keep],jdx[keep]),1)
        return x

    def summary(x):
        s=x.std(axis=0,ddof=1) if num_maps>1 else np.full(x.shape[1:],np.nan)
        return pd.DataFrame(x.mean(axis=0),index=names,columns=names),pd.DataFrame(s,index=names,columns=names)

    d=counts(ana[ana["event_type"].isin(["d","a"])])
    j=counts(clado[clado["event_type"]=="j"])
    all_mean,all_sd=summary(d+j)
    d_mean,d_sd=summary(d)
    j_mean,j_sd=summary(j)
    totals=pd.concat([clado[["map","event_type"]],ana[["map","event_type"]]])
    by_map=totals.groupby(["map","event_type"]).size().unstack(fill_value=0).reindex(
        index=range(1,num_maps+1),columns=list(_EVENT_TYPES if "t" in set(totals["event_type"]) else _EVENT_TYPES[:-1]),fill_value=0)
    return dict(all_dispersals_counts_fromto_means=all_mean,all_dispersals_counts_fromto_sds=all_sd,
                ana_dispersals_counts_fromto_means=d_mean,ana_dispersals_counts_fromto_sds=d_sd,
                clado_dispersals_counts_fromto_means=j_mean,clado_dispersals_counts_fromto_sds=j_sd,
                counts_list=by_map, event_means=by_map.mean(),event_sds=by_map.std())


def simulate_source_areas_ana_clado(results):
    """Return tables carrying the source draws performed by native BSM."""
    return dict(clado_events_tables=results.clado_events_tables,ana_events_tables=results.ana_events_tables)


def get_inputs_for_stochastic_mapping(fit):
    return fit


def run_bsm(fit, num_maps=100, mode="uniformization", num_threads=1, seed=12345):
    if mode not in ("uniformization","parallel_substream","rejection"):
        raise ValueError("unknown BSM sampling mode")
    for value,name,limit in ((num_maps,"num_maps",1_000_000),(num_threads,"num_threads",256),
                             (seed,"seed",(1<<53)-1)):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or not (0 if name == "seed" else 1) <= int(value) <= limit
        ):
            raise ValueError(f"invalid {name}")
    num_maps, num_threads, seed = int(num_maps), int(num_threads), int(seed)
    m=fit.model
    if getattr(m,"_contract_native",False):
        from .uncertainty.simulate import sampler_buffers
        m._require_tier1a()
        m._validate_contract_params(fit.params,force=True)
        with PreparedSampler(sampler_buffers(m,fit.params),m.root_sampling_probs,0) as sampler:
            raw=sampler.draw(1,num_maps,seed,num_threads,bsm=True,rejection=mode=="rejection")
    else:
        steps,layers,hooks,fixes,segments=_steps(m,fit.params)
        raw=kernel(26,[m.tree.tip_count,m.tree.root,num_maps,seed,num_threads,mode=="rejection"],
                   *_states(m.states),m.tips,steps,m.tree.ages,layers,hooks,fixes,segments,
                   *_layers(m,fit.params),precision=m.precision)
    clado,ana=_frame(raw[0],_CLADO,m),_frame(raw[1],_ANA,m)
    return BSMResult(clado,ana,count_ana_clado_events(clado,ana,m.areas,num_maps),num_maps)


def runBSM(res, stochastic_mapping_inputs_list=None, nummaps_goal=100, *, bsm_mode="uniformization",
           num_threads=1, seedval=12345):
    return run_bsm(res,nummaps_goal,bsm_mode,num_threads,seedval)
