"""NumPy interfaces for the shared BioGeoBEARS numerical kernels."""
from functools import lru_cache
from math import comb

import numpy as np
from scipy import sparse

from ._ffi import _EMPTY, kernel


def _checked_area_index(value):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError("state area indices must be integers, not booleans")
    if isinstance(value, (int, np.integer)):
        idx = int(value)
    elif isinstance(value, (float, np.floating)):
        if not np.isfinite(value) or not float(value).is_integer():
            raise ValueError(f"non-integer state area index {value!r}")
        idx = int(value)
    else:
        raise ValueError(f"invalid state area index {value!r}")
    if idx < 0:
        raise ValueError(f"negative state area index {value!r}")
    return idx


def _normalize_state(state):
    if len(state) == 0:
        return (-1,)
    if len(state) == 1:
        first = state[0]
        if not isinstance(first, (bool, np.bool_)) and isinstance(
            first, (int, np.integer, float, np.floating)
        ):
            if np.isnan(first) or (np.isfinite(first) and float(first) == -1.0):
                return (-1,)
    return tuple(_checked_area_index(i) for i in state)


def _states(states):
    normalized = [_normalize_state(x) for x in states]
    return (
        np.fromiter((i for s in normalized for i in s), float),
        np.asarray([len(s) for s in normalized], dtype=float),
    )


def _square(value):
    q = np.asarray(value, dtype=float)
    if q.ndim != 2 or q.shape[0] != q.shape[1] or not q.shape[0]:
        raise ValueError("expected a nonempty square matrix")
    return q


def rcpp_areas_list_to_states_list(areas, maxareas=None, include_null_range=True):
    n = len(areas)
    s = kernel(1, n, n if maxareas is None else maxareas, include_null_range)
    return [tuple(a.astype(int)) for a in s]


def Rcpp_combn_zerostart(n, m):
    return kernel(2, n, m, comb(n, m))[0].astype(int).reshape((m, -1), order="F")


def rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat=None,
                              include_null_range=True, normalize_TF=True, makeCOO_TF=False,
                              min_precision=1e-26, *, precision="strict_f64"):
    n = len(areas_list)
    dmat = _square(dmat)
    if dmat.shape != (n,n):
        raise ValueError("dispersal matrix does not match area count")
    amat = np.zeros_like(dmat) if amat is None else _square(amat)
    out = kernel(3, *_states(states_list), dmat, elist, amat,
                 [n, normalize_TF, makeCOO_TF, min_precision if makeCOO_TF else 1e-15],
                 precision=precision)
    return dict(ia=out[0].astype(int), ja=out[1].astype(int), a=out[2]) if makeCOO_TF else (
        out[0].reshape((len(states_list), len(states_list)), order="F"))


def maxent(constr, states, prior=None, tol=1e-7, lambda_=False, *, precision="strict_f64"):
    states = np.atleast_2d(np.asarray(states, dtype=float))
    if states.shape[1] == 1 and states.shape[0] > 1:
        states = states.T
    targets = np.atleast_2d(np.asarray(constr, dtype=float))
    k,n = states.shape
    if targets.shape[1] != k or np.any(states <= 0) or not np.isfinite(states).all():
        raise ValueError("invalid MaxEnt constraint or state dimensions")
    prior = np.full(n, 1/n) if prior is None else np.asarray(prior, dtype=float)
    prior = np.broadcast_to(prior, (len(targets), n))
    if (prior<0).any() or (prior>1).any() or (prior.sum(axis=1)<=0).any():
        raise ValueError("invalid MaxEnt prior")
    prior = prior / prior.sum(axis=1, keepdims=True)
    fits = [kernel(13, states.T, n, k, target, p, tol, precision=precision)
            for target,p in zip(targets,prior)]
    probabilities = np.asarray([a[0] for a in fits])
    result = dict(prob=probabilities, entropy=np.array([a[1][0] for a in fits]),
                  iter=np.array([int(a[1][1]) for a in fits]), moments=np.array([a[2] for a in fits]),
                  constr=targets, prior=prior)
    if lambda_:
        result["lambda"] = np.array([
            np.linalg.lstsq(np.column_stack([np.ones(n), states.T, np.log(p)]), np.log(v), rcond=None)[0]
            for p,v in zip(prior,probabilities)])
    if len(targets) == 1:
        result = {key: value[0] for key,value in result.items()}
    result["states"] = states
    return result


@lru_cache(maxsize=128)
def _priors(areas, width, constraints, precision):
    out = [np.zeros((areas,areas),order="F") for _ in range(4)]
    for kind,p in enumerate(constraints):
        if not 1e-5 <= p <= 0.99999:
            raise ValueError("size prior constraint must be between 0.00001 and 0.99999")
        for i in range(1,width+1):
            n = i // 2 if kind == 1 else i
            if n:
                out[kind][i-1,:n] = np.round(maxent(p*(n+1), np.arange(1,n+1),
                                                     precision=precision)["prob"],3)
            elif width == 1:
                out[kind][0,0] = 1
    allowed = np.maximum.reduce(out) > 0
    sizes = np.max(allowed * np.arange(1,areas+1), axis=1)
    return (*out, sizes)


def _clado_args(l, s, v, j, y, dmat, priors, maxmin, precision):
    width = max(map(len,l))
    n = width if dmat is None else _square(dmat).shape[0]
    dmat = np.ones((n,n), order="F") if dmat is None else dmat
    default = _priors(n, width, (0.0001,)*4, precision)
    padded = []
    for x,base in zip(priors,default):
        if x is None:
            padded.append(base)
        else:
            p = np.zeros((n,n), order="F")
            x = np.asarray(x)
            nr,nc = min(x.shape[0],n), min(x.shape[1],n)
            p[:nr,:nc] = x[:nr,:nc]
            padded.append(p)
    if maxmin is not None:
        maxmin = np.asarray(maxmin)
        maxmin = np.pad(maxmin, (0,max(0,n-len(maxmin))))
    return (*_states(l), [s,v,j,y], dmat, *padded, default[4] if maxmin is None else maxmin, n, _EMPTY)


def _clado(op, left, right, l, s, v, j, y, dmat, priors, maxmin, denom, kind, precision):
    args = _clado_args(l,s,v,j,y,dmat,priors,maxmin,precision)
    return kernel(op,*args,*((left,right,denom,kind) if op == 5 else ()),precision=precision)


def rcpp_calc_anclikes_sp(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=None,
                          maxent01s=None, maxent01v=None, maxent01j=None, maxent01y=None,
                          max_minsize_as_function_of_ancsize=None, Rsp_rowsums=None, printmat=False,
                          *, precision="strict_f64"):
    denom = np.ones(len(l)) if Rsp_rowsums is None else Rsp_rowsums
    return _clado(5,Rcpp_leftprobs,Rcpp_rightprobs,l,s,v,j,y,dmat,
                  (maxent01s,maxent01v,maxent01j,maxent01y),max_minsize_as_function_of_ancsize,
                  denom,0,precision)[0]


def rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=None,
                                  maxent01s=None, maxent01v=None, maxent01j=None, maxent01y=None,
                                  max_minsize_as_function_of_ancsize=None, printmat=False,
                                  *, precision="strict_f64"):
    return _clado(5,Rcpp_leftprobs,Rcpp_rightprobs,l,s,v,j,y,dmat,
                  (maxent01s,maxent01v,maxent01j,maxent01y),max_minsize_as_function_of_ancsize,
                  _EMPTY,1,precision)[0]


def rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1,
                                           dmat=None, maxent01s=None, maxent01v=None, maxent01j=None,
                                           maxent01y=None, max_minsize_as_function_of_ancsize=None,
                                           printmat=False, *, precision="strict_f64"):
    out = _clado(4,Rcpp_leftprobs,Rcpp_rightprobs,l,s,v,j,y,dmat,
                 (maxent01s,maxent01v,maxent01j,maxent01y),max_minsize_as_function_of_ancsize,
                 _EMPTY,0,precision)
    return [a.astype(int) for a in out[:3]] + out[3:]


def rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1,
                                  dmat=None, maxent01s=None, maxent01v=None, maxent01j=None,
                                  maxent01y=None, max_minsize_as_function_of_ancsize=None,
                                  printmat=False, *, precision="strict_f64"):
    out = _clado(5,Rcpp_leftprobs,Rcpp_rightprobs,l,s,v,j,y,dmat,
                 (maxent01s,maxent01v,maxent01j,maxent01y),max_minsize_as_function_of_ancsize,
                 _EMPTY,2,precision)
    n=len(l)
    return ([a.astype(int) for a in out[:n]], [a.astype(int) for a in out[n:2*n]], out[2*n:])


def rcpp_calc_anclikes_sp_using_COOprobs(Rcpp_leftprobs, Rcpp_rightprobs, RCOO_left_i_list,
                                        RCOO_right_j_list, RCOO_probs_list, Rsp_rowsums,
                                        printmat=False, *, precision="strict_f64"):
    lengths = [len(a) for a in RCOO_left_i_list]
    if lengths != [len(a) for a in RCOO_right_j_list] or lengths != [len(a) for a in RCOO_probs_list]:
        raise ValueError("COO probability lists must have the same shape")
    return kernel(8,Rcpp_leftprobs,Rcpp_rightprobs,np.repeat(np.arange(len(lengths)),lengths),
                  np.concatenate(RCOO_left_i_list),np.concatenate(RCOO_right_j_list),
                  np.concatenate(RCOO_probs_list),Rsp_rowsums,precision=precision)[0]


def rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar, numstates=None, printmat=False):
    anc,left,right,weights = COO_weights_columnar
    if numstates is None:
        numstates = max(map(np.max,(anc,left,right)))+1
    return kernel(6,anc,weights,numstates)[0]


def rcpp_calc_splitlikes_using_cooWeights_columnar(Rcpp_leftprobs, Rcpp_rightprobs, COO_weights_columnar,
                                                  Rsp_rowsums, printmat=False, *, precision="strict_f64"):
    return kernel(7,Rcpp_leftprobs,Rcpp_rightprobs,*COO_weights_columnar,Rsp_rowsums,
                  precision=precision)[0]


rcpp_calc_splitlikes_using_COOweights_columnar = rcpp_calc_splitlikes_using_cooWeights_columnar


def rcpp_convolve(a, b):
    return kernel(9,a,b)[0]


def rcpp_mult2probvect(a, b):
    return kernel(10,a,b)[0]


def expokit_dgpadm_Qmat2(times, Qmat, transpose_needed=True, *, precision="strict_f64"):
    q = _square(Qmat)
    n = len(q)
    if transpose_needed:
        q = q.T
    return kernel(11,q,n,times,precision=precision)[0].reshape((n,n),order="C")


def coo2crs(ia, ja, a, n=None, transpose_needed=False):
    if n is None:
        n = max(np.max(ia),np.max(ja))
    ptr,cols,values,dim = kernel(14,ia,ja,a,n,transpose_needed)
    return dict(ia=ptr.astype(int),ja=cols.astype(int),a=values,dimension=dim.astype(int))


def _action(Qmat, time, vector, transpose, transform, n, precision):
    if sparse.issparse(Qmat):
        csr = Qmat.tocsr()
        return kernel(17,csr.shape[0],csr.indptr+1,csr.indices+1,csr.data,time,vector,transpose,
                      precision=precision)[0]
    if not transform:
        if isinstance(Qmat, dict):
            ia,ja,a = Qmat["ia"],Qmat["ja"],Qmat["a"]
        else:
            ia,ja,a = np.asarray(Qmat).T
        return kernel(18,n,ia,ja,a,time,vector,transpose,precision=precision)[0]
    q = _square(Qmat)
    return kernel(12,q.T if transpose else q,len(q),time,vector,precision=precision)[0]


def _expokit(Qmat, t, vector, transpose, transform, n, anorm, check_for_0_rows,
             precision, symmetric):
    original = _square(Qmat) if transform else None
    if transform:
        q = original.copy()
        if vector is None and precision!="strict_f64" and check_for_0_rows:
            zeros = ~(q!=0).any(axis=1)
            if zeros.any():
                pad=1e-15 if check_for_0_rows is True else float(check_for_0_rows)
                q[zeros,:] = pad
                ix=np.flatnonzero(zeros)
                q[np.ix_(ix,ix)] = -pad*(len(q)-1)
        if anorm is None:
            anorm=np.linalg.norm(q,ord=1)
        q=q.T if transpose else q
        n=len(q)
        rows,cols=np.nonzero(q)
        values=q[rows,cols]
    else:
        if isinstance(Qmat,dict):
            ia,ja,values=Qmat["ia"],Qmat["ja"],Qmat["a"]
        else:
            ia,ja,values=np.asarray(Qmat).T
        rows,cols=np.asarray(ia,dtype=int)-1,np.asarray(ja,dtype=int)-1
        q=np.zeros((n,n),order="F")
        q[rows,cols]=values
        if anorm is None:
            anorm=np.linalg.norm(q if transpose else q.T,ord=1)
    if precision=="strict_f64":
        if symmetric and vector is not None:
            return kernel(12,q,n,t,vector,precision=precision)[0][None,:]
        res=kernel(11,q,n,t,precision=precision)[0]
    else:
        def action(v):
            return kernel(16,q,n,t,v,n-1,.01,anorm,symmetric,rows+1,cols+1,values,
                          precision=precision)[0]
        if symmetric and vector is not None:
            return action(vector)[None,:]
        res=np.concatenate([action(v) for v in np.eye(n)])
    if not symmetric and vector is not None:
        return dict(res=res,w=res[-n:])
    return res.reshape(n,n)


def expokit_dgexpv_Qmat(Qmat, t=2.1, inputprobs_for_fast=None, transpose_needed=True,
                       transform_to_coo_TF=True, coo_n=None, anorm=None, check_for_0_rows=True,
                       *, precision="strict_f64"):
    return _expokit(Qmat,t,inputprobs_for_fast,transpose_needed,transform_to_coo_TF,
                    coo_n,anorm,check_for_0_rows,precision,False)


def expokit_dmexpv_Qmat(Qmat, t=2.1, inputprobs_for_fast=None, transpose_needed=True,
                       transform_to_coo_TF=True, coo_n=None, anorm=None, check_for_0_rows=True,
                       *, precision="strict_f64"):
    return _expokit(Qmat,t,inputprobs_for_fast,transpose_needed,transform_to_coo_TF,
                    coo_n,anorm,check_for_0_rows,precision,True)


def expokit_dgexpv(mat, t, vector, transpose_needed=True, *, precision="strict_f64"):
    """Apply a dense or sparse generator to a vector without densifying sparse input."""
    return {"output_probs": _action(mat,t,vector,transpose_needed,True,None,precision)}
