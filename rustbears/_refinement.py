"""Absolute scaled-coordinate differences, including parameters at a bound.

The cached SearchResult is only the optimizer's own completed default-start
multi-basin phase. No external reference parameters enter this continuation.
"""
import numpy as np
from scipy.optimize import minimize
from ._optimization import SearchResult


def refine_completed_search(evaluate,names,start,bounds,precision,maxiter,base):
    if not names or base.proposal is None:return base
    bounds=np.asarray(bounds,dtype=float);lo,hi=bounds.T
    scales=np.array([.05 if k=='j' else max(abs(v),1e-3 if k in ('d','e','a') else .1)
                     for k,v in zip(names,start)])
    low,high=lo/scales,hi/scales
    best={p:[np.inf,np.array(start,dtype=float)] for p in ('strict_f64','legacy_f32')}
    evals=0;traces=[]
    def score(theta,p):
        nonlocal evals
        evals+=1;theta=np.clip(theta,lo,hi)
        try:ll=evaluate(theta,p)
        except ValueError as error:
            if not any(x in str(error) for x in ('invalid action input','negative likelihood','node_nan:',
                       'invalid likelihood row sum','stratified downpass NaN')):raise
            return np.inf
        value=-float(ll)
        if not np.isfinite(value):return np.inf
        if value<best[p][0]:best[p]=[value,theta.copy()]
        return value
    def objective(u,p):
        value=score(u*scales,p)
        if np.isfinite(value):return value
        f,theta=best[p]
        return (f if np.isfinite(f) else 0.)+1e6+1e5*float(np.sum((u-theta/scales)**2))
    def local(theta,p):
        # An absolute step in scaled physical coordinates cannot collapse as
        # e approaches 1e-12. Bound-clipped differences use their actual span.
        step=1e-4 if p=='legacy_f32' else 1e-5
        def jac(u):
            gradient=np.empty(len(u))
            for i in range(len(u)):
                left=u.copy();right=u.copy()
                left[i]=max(low[i],u[i]-step);right[i]=min(high[i],u[i]+step)
                gradient[i]=(objective(right,p)-objective(left,p))/(right[i]-left[i]) if right[i]>left[i] else 0.
            return gradient
        result=minimize(lambda u:objective(u,p),np.clip(theta,lo,hi)/scales,jac=jac,method='L-BFGS-B',
            bounds=list(zip(low,high)),options=dict(maxiter=maxiter,ftol=1e-13,gtol=1e-7,maxls=50,maxfun=max(15000,50*maxiter)))
        score(result.x*scales,p)
        traces.append(dict(precision=p,coordinates='scaled-linear-absolute-difference',absolute_step=step,
            success=bool(result.success),message=str(result.message),fun=float(result.fun),iterations=int(getattr(result,'nit',0))))
    proposal=dict(base.proposal)
    theta=np.asarray(proposal['theta'],dtype=float)
    score(theta,'strict_f64');score(base.theta,'strict_f64')
    local(best['strict_f64'][1],'strict_f64')
    proposal.update(theta=best['strict_f64'][1].copy(),lnL=-best['strict_f64'][0],
        traces=[*proposal.get('traces',[]),*traces],message=proposal['message']+'; absolute-step boundary refinement')
    if precision=='legacy_f32':
        score(base.theta,precision);score(proposal['theta'],precision)
        local(best[precision][1],precision)
    final=best[precision]
    return SearchResult(final[1],base.evaluations+evals,base.converged or any(t['success'] for t in traces if t['precision']==precision),
        base.message+'; absolute-step boundary refinement: '+traces[-1]['message'],proposal)
