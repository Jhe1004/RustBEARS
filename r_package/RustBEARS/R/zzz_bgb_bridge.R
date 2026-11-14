# GPL >= 2. Private bridge; public signatures remain in the original R shell.
.bgb_precision <- function() {
    value<-tolower(gsub('_','',getOption('RustBEARS.precision',getOption('BioGeoBEARS.precision','StrictF64')),fixed=TRUE))
    if(length(value)!=1L||!value %in% c('legacyf32','strictf64'))stop('Invalid RustBEARS.precision / BioGeoBEARS.precision option')
    if(value=='strictf64')return(1L)
    fused<-getOption('RustBEARS.legacy_fma',getOption('BioGeoBEARS.legacy_fma',grepl('aarch64|arm64',R.version$arch)))
    if(!is.logical(fused)||length(fused)!=1L||is.na(fused))stop('Invalid RustBEARS.legacy_fma / BioGeoBEARS.legacy_fma option')
    if(fused)2L else 0L
}
.bgb_bsm_engine <- function() getOption('RustBEARS.BSM.engine', getOption('BioGeoBEARS.BSM.engine'))
.bgb_detect_cores <- function(requested=NULL) {
    detected<-parallel::detectCores()
    if(length(detected)==1L&&is.finite(detected)&&detected>=1)return(as.integer(detected))
    # Restricted containers may not expose the hardware count. Honor a valid
    # explicit budget instead of making an NA comparison or creating a pool.
    if(length(requested)==1L&&is.numeric(requested)&&is.finite(requested)&&requested>=1&&requested<=256)
        return(as.integer(requested))
    1L
}
.bgb_stratified_engine <- function() getOption('RustBEARS.stratified.engine', getOption('BioGeoBEARS.stratified.engine'))
trim <- function(s) base::trimws(s)
.bgb_composite_optim <- function(par, fn, ..., method, lower, upper, control) {
    .bgb_composite_stats(par, fn, ..., method=method, lower=lower, upper=upper, control=control)
}
.bgb_kernel <- function(op,args) .Call(.bgb_c_bgb_kernel,as.integer(op),lapply(args,as.double),.bgb_precision())
.bgb_state_buffers <- function(x)list(unlist(x,use.names=FALSE),lengths(x))
.bgb_clado_call <- function(.native_symbol,...) {
    symbol<-.native_symbol
    a<-list(...)
    if(symbol=='cpp_areas_list_to_states_list')return(lapply(.bgb_kernel(1,list(length(a[[1]]),a[[2]],a[[3]])),as.integer))
    if(symbol=='cpp_combn_zerostart')return(matrix(as.integer(.bgb_kernel(2,a)[[1]]),nrow=a[[2]]))
    if(symbol %in% c('cpp_states_list_to_DEmat','cpp_states_list_to_DEmat_COO')) {
        coo<-symbol=='cpp_states_list_to_DEmat_COO'
        args<-c(.bgb_state_buffers(a[[2]]),a[3:5],list(c(length(a[[1]]),a[[6]],coo,if(coo)a[[7]] else 1e-15)))
        r<-.bgb_kernel(3,args)
        if(!coo)return(matrix(r[[1]],nrow=length(a[[2]])))
        return(list(ia=as.integer(r[[1]]),ja=as.integer(r[[2]]),a=r[[3]]))
    }
    if(symbol %in% c('convolve3cpp','mult2probvect'))return(.bgb_kernel(if(symbol=='convolve3cpp')9 else 10,a)[[1]])
    if(symbol=='cpp_calc_rowsums_for_COOweights_columnar')return(.bgb_kernel(6,a)[[1]])
    if(symbol=='cpp_calc_splitlikes_using_COOweights_columnar')return(.bgb_kernel(7,a)[[1]])
    if(symbol=='cpp_calc_anclikes_sp_using_COOprobs') {
        anc<-rep(seq_along(a[[4]])-1L,lengths(a[[4]]))
        return(.bgb_kernel(8,c(a[2:3],list(anc,unlist(a[[4]],use.names=FALSE),unlist(a[[5]],use.names=FALSE),unlist(a[[6]],use.names=FALSE),a[[7]])))[[1]])
    }
    na<-nrow(a[[9]])
    priors<-lapply(a[10:13],function(x){if(identical(dim(x),c(na,na)))return(x);out<-matrix(0,na,na);nr<-min(na,nrow(x));nc<-min(na,ncol(x));out[seq_len(nr),seq_len(nc)]<-x[seq_len(nr),seq_len(nc),drop=FALSE];out})
    maxmin<-a[[14]];if(length(maxmin)<na)maxmin<-c(maxmin,rep(0,na-length(maxmin)))
    common<-c(.bgb_state_buffers(a[[4]]),list(as.numeric(unlist(a[5:8])),a[[9]]),priors,list(maxmin,na,numeric(0)))
    if(symbol=='cpp_calc_anclikes_sp_COOweights_faster') {
        r<-.bgb_kernel(4,common);r[1:3]<-lapply(r[1:3],as.integer);return(r)
    }
    kind<-switch(symbol,cpp_calc_anclikes_sp=0,cpp_calc_anclikes_sp_rowsums=1,cpp_calc_anclikes_sp_COOprobs=2,stop('Unknown clado native symbol'))
    r<-.bgb_kernel(5,c(common,list(a[[2]],a[[3]],if(kind==0)a[[15]] else numeric(0),kind)))
    if(kind!=2)return(r[[1]])
    n<-length(a[[4]]);list(lapply(r[seq_len(n)],as.integer),lapply(r[n+seq_len(n)],as.integer),r[2*n+seq_len(n)])
}
.bgb_itscale <- function(states,n,k,target,prior,tol) {
    r<-.bgb_kernel(13,list(as.numeric(t(states)),n,k,target,prior,tol))
    list(prob=r[[1]],entropy=r[[2]][1],niter=as.integer(r[[2]][2]),moments=r[[3]])
}
.bgb_coomat <- function(ia,ja,a,n) {
    out<-matrix(0,n,n)
    if(length(a))for(i in seq_along(a))out[ia[i],ja[i]]<-a[i]
    out
}
.bgb_rexp_call <- function(.native_symbol,...) {
    symbol<-.native_symbol
    a<-list(...)
    if(symbol=='R_rexpokit_as_coo') {
        q<-a[[1]];ij<-which(t(q)!=0,arr.ind=TRUE)
        return(unname(cbind(ij[,2],ij[,1],q[cbind(ij[,2],ij[,1])])) )
    }
    if(symbol=='R_dgpadm') {
        r<-.bgb_kernel(11,list(a[[4]],a[[2]],a[[3]]))[[1]]
        return(list(wsp=r,ind=1L))
    }
    n<-a[[1]];q<-.bgb_coomat(a[[11]],a[[12]],a[[13]],n)
    action<-function(v){
        if(.bgb_precision()!=1L).bgb_kernel(16,list(q,n,a[[3]],v,a[[2]],a[[5]],a[[6]],symbol%in%c('R_dmexpv','R_mydmexpv'),a[[11]],a[[12]],a[[13]]))[[1]]
        else .bgb_kernel(12,list(q,n,a[[3]],v))[[1]]
    }
    if(symbol %in% c('R_mydmexpv','R_mydgexpv'))return(action(a[[4]]))
    if(symbol %in% c('R_dmexpv','R_dgexpv')) {
        p<-if(.bgb_precision()!=1L)as.numeric(vapply(seq_len(n),function(j){v<-numeric(n);v[j]<-1;action(v)},numeric(n)))else .bgb_kernel(11,list(q,n,a[[3]]))[[1]]
        return(list(res=p,w=matrix(p,n,n)[,n]))
    }
    stop('Unknown rexpokit native symbol: ',symbol)
}
# The kexpmv subset actually used by BioGeoBEARS; CRS is 1 based.
expokit_dgexpv <- function(mat,t,vector,transpose_needed=TRUE,transform_to_crs=TRUE,crs_n=NULL,anorm=NULL,mxstep=10000,tol=1e-10) {
    if(transform_to_crs){
        q<-as.matrix(mat);if(transpose_needed)q<-t(q)
        out<-.bgb_kernel(12,list(q,nrow(q),t,vector))[[1]]
    }else{
        n<-if(is.null(crs_n))length(mat$ia)-1L else crs_n
        out<-.bgb_kernel(17,list(n,mat$ia,mat$ja,mat$a,t,vector,transpose_needed))[[1]]
    }
    list(output_probs=out)
}
expokit_dmexpv <- function(mat,t,vector,transpose_needed=TRUE,transform_to_crs=TRUE,crs_n=NULL,anorm=NULL,mxstep=10000,tol=1e-10)expokit_dgexpv(mat,t,vector,transpose_needed,transform_to_crs,crs_n,anorm,mxstep,tol)
