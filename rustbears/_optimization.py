"""Reference-free, scale-aware bounded likelihood optimization."""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize, differential_evolution, dual_annealing


@dataclass
class SearchResult:
    theta: np.ndarray
    evaluations: int
    converged: bool
    message: str
    proposal: object


def optimize_likelihood(evaluate, names, start, bounds, height, precision, optimizer, maxiter, seed):
    """Uses only the evaluator, supplied bounds/start and tree scale."""
    x0=np.array(start,dtype=float);bounds=np.asarray(bounds,dtype=float).reshape((-1,2))
    if not len(names):
        evaluate(x0,precision)
        return SearchResult(x0,1,True,'fixed parameters',None)
    lo,hi=bounds.T
    if np.all(lo==hi):
        ll=evaluate(lo.copy(),precision)
        if not np.isfinite(ll):
            raise ValueError('no finite likelihood found at equal parameter bounds')
        return SearchResult(lo.copy(),1,True,'equal bounds fixed parameters',None)
    scales=np.array([.05 if k=='j' else max(abs(v),1e-3 if k in ('d','e','a') else .1)
                     for k,v in zip(names,x0)])
    logs=np.array([k in ('d','e','a') and low>0 for k,low in zip(names,lo)])
    loglo=lo.copy();loghi=hi.copy();loglo[logs]=np.log(lo[logs]);loghi[logs]=np.log(hi[logs])
    best={p:[np.inf,x0.copy()] for p in ('strict_f64','legacy_f32')}
    evals=0;traces=[]
    def score(theta,p):
        nonlocal evals
        evals+=1
        try:
            ll=evaluate(np.clip(theta,lo,hi),p)
        except ValueError as error:
            if not any(s in str(error) for s in ('invalid action input','negative likelihood','node_nan:',
                      'invalid likelihood row sum','stratified downpass NaN')):raise
            return np.inf
        value=-float(ll)
        if not np.isfinite(value):return np.inf
        if value<best[p][0]:best[p]=[value,np.clip(theta,lo,hi).copy()]
        return value
    def encoded(theta,kind):
        z=np.array(theta,copy=True)
        if kind=='linear':return z/scales
        z[logs]=np.log(np.maximum(z[logs],lo[logs]));return z
    def decoded(z,kind):
        if kind=='linear':return z*scales
        theta=z.copy();theta[logs]=np.exp(theta[logs]);return theta
    def objective(z,p,kind):
        v=score(decoded(z,kind),p)
        if np.isfinite(v):return v
        anchor=encoded(best[p][1],kind)
        base=best[p][0] if np.isfinite(best[p][0]) else 0.
        return base+1e6+1e5*float(np.sum((z-anchor)**2))
    def local(theta,p,kind):
        start=np.clip(theta,lo,hi).copy()
        if kind=='log':start[logs]=np.maximum(start[logs],np.maximum(lo[logs],1e-8))
        options=dict(maxiter=maxiter,ftol=1e-13,gtol=1e-7,maxls=50,maxfun=max(15000,50*maxiter))
        options['finite_diff_rel_step']=1e-4 if p=='legacy_f32' else 1e-5
        result=minimize(lambda z:objective(z,p,kind),encoded(start,kind),method='L-BFGS-B',jac='3-point',
                        bounds=list(zip(lo/scales,hi/scales)) if kind=='linear' else list(zip(loglo,loghi)),options=options)
        score(decoded(result.x,kind),p)
        traces.append(dict(precision=p,coordinates=kind,success=bool(result.success),message=str(result.message),
                           fun=float(result.fun),iterations=int(getattr(result,'nit',0))))
        return result
    if optimizer=='L-BFGS-B':
        rc=float(np.clip(.15/max(float(height),1.),1e-4,.05))
        probes=[x0.copy()]
        for rates in ({'d':rc,'e':rc,'j':.01},{'d':.2*rc,'e':1e-7,'j':.02},
                      {'d':.1*rc,'e':.1*rc,'j':1e-4}):
            probe=np.clip([rates.get(k,x0[i]) for i,k in enumerate(names)],lo,hi)
            if not any(np.array_equal(probe,x) for x in probes):probes.append(probe)
        ranked=sorted(((score(x,'strict_f64'),i,x) for i,x in enumerate(probes)),key=lambda x:(x[0],x[1]))
        finite=[x for x in ranked if np.isfinite(x[0])]
        if not finite:raise ValueError('no finite likelihood found at the initial point or deterministic probes')
        chosen=[finite[0]]
        if len(finite)>1:
            ei=names.index('e') if 'e' in names else None
            other=[x for x in finite[1:] if ei is not None and (x[2][ei]>1e-4)!=(finite[0][2][ei]>1e-4)]
            chosen.append((other or finite[1:])[0])
        for _,_,x in chosen:local(x,'strict_f64','linear')
        local(best['strict_f64'][1],'strict_f64','log')
        proposal=dict(theta=best['strict_f64'][1].copy(),lnL=-best['strict_f64'][0],
                      converged=any(x['success'] for x in traces),message='StrictF64 scaled-linear multi-basin + log refinement',
                      probe_scale=rc,probes=[dict(theta=x.tolist(),initial=i==0) for i,x in enumerate(probes)],
                      starts=[x[2].tolist() for x in chosen],traces=list(traces))
        if precision=='legacy_f32':
            score(x0,precision);score(proposal['theta'],precision)
            local(best[precision][1],precision,'linear')
        final=best[precision]
        return SearchResult(final[1],evals,any(t['success'] for t in traces if t['precision']==precision),
                            'Deterministic multi-basin search; '+traces[-1]['message'],proposal)
    score(x0,precision)
    active=loglo<loghi
    base_z=encoded(x0,'log')
    def sub_fn(z_sub):
        z_full=base_z.copy()
        z_full[active]=z_sub
        return objective(z_full,precision,'log')
    sub_bounds=list(zip(loglo[active],loghi[active]))
    if optimizer=='differential_evolution':
        result=differential_evolution(sub_fn,sub_bounds,x0=base_z[active],seed=seed,maxiter=maxiter)
    elif optimizer=='dual_annealing':
        result=dual_annealing(sub_fn,sub_bounds,x0=base_z[active],seed=seed,maxiter=maxiter)
    else:raise ValueError('unknown optimizer')
    if not np.isfinite(best[precision][0]):raise ValueError('no finite likelihood found within the parameter bounds')
    return SearchResult(best[precision][1],evals,bool(result.success),str(result.message),None)
