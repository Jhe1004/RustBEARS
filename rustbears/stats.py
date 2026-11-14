"""Akaike scores and likelihood-ratio calculations."""
import numpy as np
import pandas as pd
from scipy.stats import chi2


def lrttest(LnL_1,LnL_2,numparams1,numparams2,returnwhat="pval"):
    df=numparams1-numparams2
    d=2*(LnL_1-LnL_2)
    p=chi2.sf(d,df) if df>0 else (float(d<=0) if df==0 else np.nan)
    if returnwhat=="all":
        return pd.DataFrame([dict(alt="",null="",LnLalt=LnL_1,LnLnull=LnL_2,DFalt=numparams1,
                                  DFnull=numparams2,DF=df,Dstatistic=d,pval=p,test="chi-squared",tail="one-tailed")])
    return float(p)


def _weights(values):
    a=np.asarray(values,dtype=float)
    ok=np.isfinite(a)
    if not ok.any():
        return np.full_like(a,np.nan)
    w=np.exp(-.5*(a-a[ok].min()))
    return w/np.nansum(w)


def AkaikeWeights_on_summary_table(restable,colname_to_use="AIC",add_to_table=True):
    w=_weights(restable[colname_to_use])
    if not add_to_table:
        return w
    result=restable.copy()
    result[colname_to_use+"_wt"]=w
    return result


def compare_models(results,ref_index=0):
    r=list(results)
    if not r:
        raise ValueError("at least one fitted model is required")
    ref=r[ref_index]
    table=pd.DataFrame([{"model":f.model.model.name,"lnL":f.lnL,"numparams":f.n_params,
                         **f.params,"AIC":f.aic,"AICc":f.aicc,
                         "LRT_D":2*(f.lnL-ref.lnL),"LRT_df":f.n_params-ref.n_params,
                         "LRT_pval":lrttest(f.lnL,ref.lnL,f.n_params,ref.n_params)} for f in r])
    table["AICweight"]=_weights(table["AIC"])
    table["AICcweight"]=_weights(table["AICc"])
    return table
