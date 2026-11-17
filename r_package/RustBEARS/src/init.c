/* GPL >= 2. R allocation and error handling stay on this side of the ABI. */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <inttypes.h>
#ifndef _WIN32
#include <dlfcn.h>
#endif
typedef struct Reply Reply;
extern Reply *bgb_kernel(uint32_t,const double **,const size_t *,size_t,uint32_t);
extern const char *bgb_result_error(const Reply *);
extern size_t bgb_result_count(const Reply *);
extern size_t bgb_result_len(const Reply *,size_t);
extern const double *bgb_result_data(const Reply *,size_t);
extern void bgb_result_free(Reply *);
void bgb_dgemm(int n,const double *a,const double *b,double *c) {
    double one=1.0,zero=0.0;
    F77_CALL(dgemm)("N","N",&n,&n,&n,&one,a,&n,b,&n,&zero,c,&n FCONE FCONE);
}
void bgb_dgemm_rect(int n,int cols,const double *a,const double *b,double *c) {
 double one=1.0,zero=0.0;
 F77_CALL(dgemm)("N","N",&n,&cols,&n,&one,a,&n,b,&n,&zero,c,&n FCONE FCONE);
}
int bgb_eigen(int n,const double *a,double *wr,double *wi,double *vr) {
 double *copy=(double*)malloc((size_t)n*n*sizeof(double));
 if(!copy)return -100;
 memcpy(copy,a,(size_t)n*n*sizeof(double));
 double query=0.0,dummy=0.0;int ldvl=1,lwork=-1,info=0;
 F77_CALL(dgeev)("N","V",&n,copy,&n,wr,wi,&dummy,&ldvl,vr,&n,&query,&lwork,&info FCONE FCONE);
 if(info!=0){free(copy);return info;}
 lwork=(int)query;double *work=(double*)malloc((size_t)lwork*sizeof(double));
 if(!work){free(copy);return -100;}
 memcpy(copy,a,(size_t)n*n*sizeof(double));
 F77_CALL(dgeev)("N","V",&n,copy,&n,wr,wi,&dummy,&ldvl,vr,&n,work,&lwork,&info FCONE FCONE);
 free(copy);free(work);return info;
}
int bgb_dgesv(int n,double *a,double *b) {
 int *piv=(int*)malloc((size_t)n*sizeof(int));if(!piv)return -100;
 int info=0;F77_CALL(dgesv)(&n,&n,a,&n,piv,b,&n,&info);free(piv);return info;
}
static SEXP c_bgb_m0_info(void) {
    SEXP out=PROTECT(Rf_allocVector(VECSXP,2));
    SEXP names=PROTECT(Rf_allocVector(STRSXP,2));
    SET_STRING_ELT(names,0,Rf_mkChar("phase")); SET_STRING_ELT(names,1,Rf_mkChar("rust_abi"));
    SET_VECTOR_ELT(out,0,Rf_mkString("RustBEARS_ABI_v1")); SET_VECTOR_ELT(out,1,Rf_ScalarInteger(1));
    Rf_setAttrib(out,R_NamesSymbol,names); UNPROTECT(2); return out;
}
static SEXP copy_reply(void *ptr) {
    Reply *reply=(Reply*)ptr; size_t count=bgb_result_count(reply);
    SEXP out=PROTECT(Rf_allocVector(VECSXP,(R_xlen_t)count));
    for(size_t i=0;i<count;i++) {
        size_t n=bgb_result_len(reply,i);
        SEXP x=PROTECT(Rf_allocVector(REALSXP,(R_xlen_t)n));
        if(n)memcpy(REAL(x),bgb_result_data(reply,i),n*sizeof(double));
        SET_VECTOR_ELT(out,(R_xlen_t)i,x); UNPROTECT(1);
    }
    UNPROTECT(1); return out;
}
static void free_reply(void *ptr,Rboolean jump) {(void)jump;bgb_result_free((Reply*)ptr);}
static SEXP c_bgb_kernel(SEXP opcode,SEXP args,SEXP precision) {
    if(TYPEOF(args)!=VECSXP)Rf_error("kernel arguments must be a list");
    R_xlen_t n=XLENGTH(args);
    SEXP numeric=PROTECT(Rf_allocVector(VECSXP,n));
    const double **ptrs=(const double**)R_alloc((size_t)(n?n:1),sizeof(double*));
    size_t *lens=(size_t*)R_alloc((size_t)(n?n:1),sizeof(size_t));
    for(R_xlen_t i=0;i<n;i++) {
        SEXP v=PROTECT(Rf_coerceVector(VECTOR_ELT(args,i),REALSXP));
        SET_VECTOR_ELT(numeric,i,v);ptrs[i]=REAL(v);lens[i]=(size_t)XLENGTH(v);UNPROTECT(1);
    }
    Reply *reply=bgb_kernel((uint32_t)Rf_asInteger(opcode),ptrs,lens,(size_t)n,(uint32_t)Rf_asInteger(precision));
    const char *failure_message=bgb_result_error(reply);
    if(failure_message){char message[1024];snprintf(message,sizeof(message),"%s",failure_message);bgb_result_free(reply);UNPROTECT(1);Rf_error("%s",message);}
    SEXP out=R_UnwindProtect(copy_reply,reply,free_reply,reply,NULL);
    UNPROTECT(1);return out;
}
/* Pure arithmetic helpers: safe to call from Rust worker threads. */
double bgb_sum(const double *x,size_t n){long double sum=0;for(size_t i=0;i<n;i++)sum+=(long double)x[i];return (double)sum;}
double bgb_indexed_sum(const double *x,const size_t *indices,size_t n){long double sum=0;for(size_t i=0;i<n;i++)sum+=(long double)x[indices[i]];return (double)sum;}
double bgb_prod4(double a,double b,double c,double d){long double v=a;v*=b;v*=c;v*=d;return (double)v;}
void bgb_dgemv(int n,const double *a,const double *v,double *out,int trans){double one=1.0,zero=0.0;int inc=1;const char *tr=trans?"T":"N";F77_CALL(dgemv)(tr,&n,&n,&one,a,&n,v,&inc,&zero,out,&inc FCONE);}
/* Scoped OpenBLAS limit: no R API, restored before any R allocation/error. */
#ifndef _WIN32
typedef int (*BlasGetThreads)(void);
typedef void (*BlasSetThreads)(int);
static BlasGetThreads blas_get=NULL;
static BlasSetThreads blas_set=NULL;
static int blas_resolved=0;
static void resolve_blas_threads(void){
 if(blas_resolved)return;blas_resolved=1;
 blas_get=(BlasGetThreads)dlsym(RTLD_DEFAULT,"openblas_get_num_threads");
 blas_set=(BlasSetThreads)dlsym(RTLD_DEFAULT,"openblas_set_num_threads");
}
static int blas_enter(void){resolve_blas_threads();if(!blas_get||!blas_set)return 0;int old=blas_get();if(old>1)blas_set(1);return old;}
static void blas_leave(int old){if(old>1&&blas_set)blas_set(old);}
#else
static int blas_enter(void){return 0;}
static void blas_leave(int old){(void)old;}
#endif
typedef struct ContextHandle ContextHandle;
extern ContextHandle *bgb_context_create(const double **,const size_t *,size_t,size_t);
extern ContextHandle *bgb_strat_context_create(const double **,const size_t *,size_t,size_t);
extern const char *bgb_context_error(const ContextHandle *);
extern void bgb_context_free(ContextHandle *);
extern Reply *bgb_context_evaluate(ContextHandle *,const double **,const size_t *,size_t,uint32_t);
extern Reply *bgb_context_evaluate_contract(ContextHandle *,const double **,const size_t *,size_t,uint32_t,uint32_t,const double *,size_t);
static void finalize_context(SEXP ptr){ContextHandle *p=(ContextHandle*)R_ExternalPtrAddr(ptr);R_ClearExternalPtr(ptr);if(p)bgb_context_free(p);}
static void numeric_buffers(SEXP args,const double ***ptrs,size_t **lens,SEXP *keep){
 if(TYPEOF(args)!=VECSXP)Rf_error("tree arguments must be a list");
 R_xlen_t n=XLENGTH(args);*keep=PROTECT(Rf_allocVector(VECSXP,n));
 *ptrs=(const double**)R_alloc((size_t)(n?n:1),sizeof(double*));*lens=(size_t*)R_alloc((size_t)(n?n:1),sizeof(size_t));
 for(R_xlen_t i=0;i<n;i++){SEXP x=PROTECT(Rf_coerceVector(VECTOR_ELT(args,i),REALSXP));SET_VECTOR_ELT(*keep,i,x);(*ptrs)[i]=REAL(x);(*lens)[i]=(size_t)XLENGTH(x);UNPROTECT(1);}
}
static SEXP c_bgb_tree(SEXP ptr,SEXP config,SEXP args,SEXP precision,SEXP threads){
 if(TYPEOF(ptr)!=EXTPTRSXP)Rf_error("invalid private tree handle");
 SEXP tag=R_ExternalPtrTag(ptr);
 int stratified=tag==Rf_install("bgb_M3_context");
 if(!stratified&&tag!=Rf_install("bgb_M2_context"))Rf_error("invalid private tree handle tag");
 ContextHandle *context=(ContextHandle*)R_ExternalPtrAddr(ptr);
 if(!context){
  R_RegisterCFinalizerEx(ptr,finalize_context,TRUE);
  const double **p;size_t *lens;SEXP keep;numeric_buffers(config,&p,&lens,&keep);
  size_t num_threads=(size_t)Rf_asInteger(threads);
  if(stratified)context=bgb_strat_context_create(p,lens,(size_t)XLENGTH(config),num_threads);
  else context=bgb_context_create(p,lens,(size_t)XLENGTH(config),num_threads);
  const char *failure=bgb_context_error(context);
  if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_context_free(context);UNPROTECT(1);Rf_error("%s",msg);}
  R_SetExternalPtrAddr(ptr,context);UNPROTECT(1);
 }
 const double **p;size_t *lens;SEXP keep;numeric_buffers(args,&p,&lens,&keep);
 uint32_t mode=(uint32_t)Rf_asInteger(precision);
 int blas_threads=blas_enter();
 Reply *reply=bgb_context_evaluate(context,p,lens,(size_t)XLENGTH(args),mode);
 blas_leave(blas_threads);
 const char *failure=bgb_result_error(reply);
 if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_result_free(reply);UNPROTECT(1);Rf_error("%s",msg);}
 SEXP out=R_UnwindProtect(copy_reply,reply,free_reply,reply,NULL);UNPROTECT(1);return out;
}
static SEXP c_bgb_tree_release(SEXP ptr){
 if(TYPEOF(ptr)!=EXTPTRSXP)Rf_error("invalid private tree handle");
 SEXP tag=R_ExternalPtrTag(ptr);
 if(tag!=Rf_install("bgb_M2_context")&&tag!=Rf_install("bgb_M3_context"))Rf_error("invalid private tree handle");
 finalize_context(ptr);return R_NilValue;
}
static SEXP c_bgb_tree_handle(void){
 SEXP ptr=PROTECT(R_MakeExternalPtr(NULL,Rf_install("bgb_M2_context"),R_NilValue));
 UNPROTECT(1);return ptr;
}
static SEXP c_bgb_strat_handle(void){
 SEXP ptr=PROTECT(R_MakeExternalPtr(NULL,Rf_install("bgb_M3_context"),R_NilValue));
 UNPROTECT(1);return ptr;
}
/* v1.2 additive ABI. Legacy entry points and index order remain unchanged. */
static SEXP c_bgb_tree_contract(SEXP ptr,SEXP config,SEXP args,SEXP precision,SEXP threads,SEXP contract,SEXP tips){
 if(TYPEOF(ptr)!=EXTPTRSXP)Rf_error("invalid private tree handle");
 SEXP tag=R_ExternalPtrTag(ptr);
 int stratified=tag==Rf_install("bgb_M3_context");
 if(!stratified&&tag!=Rf_install("bgb_M2_context"))Rf_error("invalid private tree handle tag");
 int code=Rf_asInteger(contract), nth=Rf_asInteger(threads);
 if(code<0 || code>2 || nth<1 || nth>256)Rf_error("invalid contract or thread count");
 ContextHandle *context=(ContextHandle*)R_ExternalPtrAddr(ptr);
 if(!context){
  R_RegisterCFinalizerEx(ptr,finalize_context,TRUE);
  const double **p;size_t *lens;SEXP keep;numeric_buffers(config,&p,&lens,&keep);
  if(stratified)context=bgb_strat_context_create(p,lens,(size_t)XLENGTH(config),(size_t)nth);
  else context=bgb_context_create(p,lens,(size_t)XLENGTH(config),(size_t)nth);
  const char *failure=bgb_context_error(context);
  if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_context_free(context);UNPROTECT(1);Rf_error("%s",msg);}
  R_SetExternalPtrAddr(ptr,context);UNPROTECT(1);
 }
 const double **p;size_t *lens;SEXP keep;numeric_buffers(args,&p,&lens,&keep);
 SEXP tip_values=PROTECT(Rf_coerceVector(tips,REALSXP));
 uint32_t mode=(uint32_t)Rf_asInteger(precision);
 int blas_threads=blas_enter();
 Reply *reply=bgb_context_evaluate_contract(context,p,lens,(size_t)XLENGTH(args),
     mode,(uint32_t)code,REAL(tip_values),(size_t)XLENGTH(tip_values));
 blas_leave(blas_threads);
 const char *failure=bgb_result_error(reply);
 if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_result_free(reply);UNPROTECT(2);Rf_error("%s",msg);}
 SEXP out=R_UnwindProtect(copy_reply,reply,free_reply,reply,NULL);UNPROTECT(2);return out;
}
typedef struct SamplerHandle SamplerHandle;
extern SamplerHandle *bgb_sampler_prepare(const double **,const size_t *,size_t,const double *,size_t,uint32_t);
extern const char *bgb_sampler_error(const SamplerHandle *);
extern void bgb_sampler_release(SamplerHandle *);
extern Reply *bgb_sampler_draw_bsm(SamplerHandle *,uint64_t,size_t,uint64_t,size_t,uint32_t);
extern Reply *bgb_sampler_simulate_tips(SamplerHandle *,uint64_t,size_t,uint64_t,size_t);
static void finalize_sampler(SEXP ptr){
 SamplerHandle *p=(SamplerHandle*)R_ExternalPtrAddr(ptr);R_ClearExternalPtr(ptr);
 if(p)bgb_sampler_release(p);
}
static SamplerHandle *checked_sampler(SEXP ptr){
 if(TYPEOF(ptr)!=EXTPTRSXP || R_ExternalPtrTag(ptr)!=Rf_install("bgb_v12_sampler"))
  Rf_error("invalid PreparedTreeSampler handle");
 SamplerHandle *p=(SamplerHandle*)R_ExternalPtrAddr(ptr);
 if(!p)Rf_error("PreparedTreeSampler is closed");
 return p;
}
static uint64_t exact_uint(SEXP value,const char *name,double low,double high){
 if(XLENGTH(value)!=1 || TYPEOF(value)==LGLSXP)Rf_error("invalid %s",name);
 double x=Rf_asReal(value);
 if(!R_FINITE(x)||x!=floor(x)||x<low||x>high)Rf_error("invalid %s",name);
 return (uint64_t)x;
}
static uint64_t exact_seed(SEXP value){
 if(TYPEOF(value)==STRSXP && XLENGTH(value)==1 && STRING_ELT(value,0)!=NA_STRING){
  const char *s=CHAR(STRING_ELT(value,0));
  if(!*s||*s=='-'||*s=='+'||*s==' '||*s=='\t')Rf_error("invalid seed");
  int base=(s[0]=='0'&&(s[1]=='x'||s[1]=='X'))?16:10;
  char *end;errno=0;unsigned long long x=strtoull(s,&end,base);
  if(errno==ERANGE||end==s||*end)Rf_error("seed must be an exact uint64 decimal or hexadecimal string");
  return (uint64_t)x;
 }
 return exact_uint(value,"seed",0,9007199254740991.0);
}
static SEXP c_bgb_seed_decimal(SEXP value){
 char text[32];snprintf(text,sizeof(text),"%" PRIu64,exact_seed(value));return Rf_mkString(text);
}
static SEXP c_bgb_sampler_prepare(SEXP args,SEXP roots,SEXP condition){
 int code=Rf_asInteger(condition);
 if(code<0||(code>2&&code!=4))Rf_error("invalid simulation conditioning");
 const double **p;size_t *lens;SEXP keep;numeric_buffers(args,&p,&lens,&keep);
 SEXP weights=PROTECT(Rf_coerceVector(roots,REALSXP));
 SEXP ptr=PROTECT(R_MakeExternalPtr(NULL,Rf_install("bgb_v12_sampler"),R_NilValue));
 R_RegisterCFinalizerEx(ptr,finalize_sampler,TRUE);
 int blas_threads=blas_enter();
 SamplerHandle *sampler=bgb_sampler_prepare(p,lens,(size_t)XLENGTH(args),REAL(weights),(size_t)XLENGTH(weights),(uint32_t)code);
 blas_leave(blas_threads);
 const char *failure=bgb_sampler_error(sampler);
 if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_sampler_release(sampler);UNPROTECT(3);Rf_error("%s",msg);}
 R_SetExternalPtrAddr(ptr,sampler);UNPROTECT(3);return ptr;
}
static SEXP checked_sampler_reply(Reply *reply){
 const char *failure=bgb_result_error(reply);
 if(failure){char msg[1024];snprintf(msg,sizeof(msg),"%s",failure);bgb_result_free(reply);Rf_error("%s",msg);}
 return R_UnwindProtect(copy_reply,reply,free_reply,reply,NULL);
}
static SEXP c_bgb_sampler_draw(SEXP ptr,SEXP start,SEXP count,SEXP seed,SEXP threads,SEXP rejection){
 SamplerHandle *sampler=checked_sampler(ptr);
 uint64_t first=exact_uint(start,"start",1,9007199254740991.0), s=exact_seed(seed);
 size_t n=(size_t)exact_uint(count,"num_maps",1,1000000), nth=(size_t)exact_uint(threads,"num_threads",1,256);
 int r=Rf_asInteger(rejection);if(r<0||r>1)Rf_error("invalid rejection mode");
 int blas_threads=blas_enter();
 Reply *reply=bgb_sampler_draw_bsm(sampler,first,n,s,nth,(uint32_t)r);
 blas_leave(blas_threads);return checked_sampler_reply(reply);
}
static SEXP c_bgb_sampler_simulate(SEXP ptr,SEXP start,SEXP count,SEXP seed,SEXP threads){
 SamplerHandle *sampler=checked_sampler(ptr);
 uint64_t first=exact_uint(start,"start",1,9007199254740991.0), s=exact_seed(seed);
 size_t n=(size_t)exact_uint(count,"n_reps",1,1000000), nth=(size_t)exact_uint(threads,"num_threads",1,256);
 int blas_threads=blas_enter();
 Reply *reply=bgb_sampler_simulate_tips(sampler,first,n,s,nth);
 blas_leave(blas_threads);return checked_sampler_reply(reply);
}
static SEXP c_bgb_sampler_release(SEXP ptr){
 if(TYPEOF(ptr)!=EXTPTRSXP || R_ExternalPtrTag(ptr)!=Rf_install("bgb_v12_sampler"))Rf_error("invalid PreparedTreeSampler handle");
 finalize_sampler(ptr);return R_NilValue;
}
/* Only the main R thread enters this RNG scope. The callbacks use no R
 * allocation, error, or interrupt APIs while the Rust stack is live. */
extern Reply *bgb_bsm_reject(const double *,size_t,size_t,size_t,double,size_t,
                             double (*)(double),size_t (*)(const double *,double *,int *,size_t));
static double bgb_bsm_wait(double rate) {
    /* stats::rexp passes scale = 1 / rate to Rmath. */
    return exp_rand() * (1.0 / rate);
}
static size_t bgb_bsm_pick(const double *weights,double *p,int *order,size_t n) {
    double input_sum=bgb_sum(weights,n),sum=0.0;
    for(size_t i=0;i<n;i++) {
        p[i]=weights[i]/input_sum;
        sum+=p[i];
        order[i]=(int)i;
    }
    for(size_t i=0;i<n;i++)p[i]/=sum;
    revsort(p,order,(int)n);
    double u=unif_rand(),mass=0.0;
    size_t j=0;
    for(;j+1<n;j++){
        mass+=p[j];
        if(u<=mass)break;
    }
    return (size_t)order[j];
}
static SEXP copy_bsm_reply(void *ptr) {
    PutRNGstate();
    if(bgb_result_error((Reply*)ptr)) return R_NilValue;
    return copy_reply(ptr);
}
static SEXP c_bgb_bsm_legacy(SEXP q,SEXP n_arg,SEXP from_arg,SEXP to_arg,SEXP duration,SEXP attempts_arg) {
    int n=Rf_asInteger(n_arg),from=Rf_asInteger(from_arg),to=Rf_asInteger(to_arg);
    int attempts=Rf_asInteger(attempts_arg);
    double t=Rf_asReal(duration);
    if(TYPEOF(q)!=REALSXP || n<=0 || from<1 || from>n || to<1 || to>n ||
       attempts<=0 || (size_t)XLENGTH(q)!=(size_t)n*(size_t)n) return R_NilValue;
    const double *values=REAL(q);
    GetRNGstate();
    Reply *reply=bgb_bsm_reject(values,(size_t)n,(size_t)(from-1),(size_t)(to-1),t,
                                (size_t)attempts,bgb_bsm_wait,bgb_bsm_pick);
    /* Protect the Rust reply even if PutRNGstate or R output allocation jumps. */
    return R_UnwindProtect(copy_bsm_reply,reply,free_reply,reply,NULL);
}

/* Match R's row-wise prod() reduction without an intermediate R matrix. */
static SEXP c_bgb_bsm_products(SEXP a,SEXP l,SEXP r,SEXP w) {
    if(TYPEOF(a)!=REALSXP || TYPEOF(l)!=REALSXP || TYPEOF(r)!=REALSXP || TYPEOF(w)!=REALSXP)
        return R_NilValue;
    R_xlen_t n=XLENGTH(a);
    if(XLENGTH(l)!=n || XLENGTH(r)!=n || XLENGTH(w)!=n) return R_NilValue;
    const double *pa=REAL(a),*pl=REAL(l),*pr=REAL(r),*pw=REAL(w);
    SEXP out=PROTECT(Rf_allocVector(REALSXP,n));
    double *values=REAL(out);
    for(R_xlen_t i=0;i<n;i++) {
        if(!R_FINITE(pa[i]) || !R_FINITE(pl[i]) || !R_FINITE(pr[i]) || !R_FINITE(pw[i])) {
            UNPROTECT(1);
            return R_NilValue;
        }
        values[i]=bgb_prod4(pa[i],pl[i],pr[i],pw[i]);
        if(!R_FINITE(values[i])) {
            UNPROTECT(1);
            return R_NilValue;
        }
    }
    UNPROTECT(1);
    return out;
}

static const R_CallMethodDef entries[]={
 {"c_bgb_seed_decimal",(DL_FUNC)&c_bgb_seed_decimal,1},
 {"c_bgb_tree_contract",(DL_FUNC)&c_bgb_tree_contract,7},
 {"c_bgb_sampler_prepare",(DL_FUNC)&c_bgb_sampler_prepare,3},
 {"c_bgb_sampler_draw",(DL_FUNC)&c_bgb_sampler_draw,6},
 {"c_bgb_sampler_simulate",(DL_FUNC)&c_bgb_sampler_simulate,5},
 {"c_bgb_sampler_release",(DL_FUNC)&c_bgb_sampler_release,1},
 {"c_bgb_bsm_products",(DL_FUNC)&c_bgb_bsm_products,4},
 {"c_bgb_bsm_legacy",(DL_FUNC)&c_bgb_bsm_legacy,6},
 {"c_bgb_m0_info",(DL_FUNC)&c_bgb_m0_info,0},
 {"c_bgb_kernel",(DL_FUNC)&c_bgb_kernel,3},
 {"c_bgb_tree_handle",(DL_FUNC)&c_bgb_tree_handle,0},
 {"c_bgb_tree",(DL_FUNC)&c_bgb_tree,5},
 {"c_bgb_strat_handle",(DL_FUNC)&c_bgb_strat_handle,0},
 {"c_bgb_tree_release",(DL_FUNC)&c_bgb_tree_release,1},
 {NULL,NULL,0}
};
void attribute_visible R_init_RustBEARS(DllInfo *dll){R_registerRoutines(dll,NULL,entries,NULL,NULL);R_useDynamicSymbols(dll,FALSE);R_forceSymbols(dll,TRUE);}
