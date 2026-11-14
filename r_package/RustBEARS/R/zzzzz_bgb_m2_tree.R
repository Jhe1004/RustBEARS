# GPL >= 2. Private orchestration; original public formals and objects survive.
.bgb_m1_calc_loglike_sp_prebyte <- calc_loglike_sp_prebyte
.bgb_m2_runtime <- new.env(parent=emptyenv())
.bgb_m2_runtime$contexts <- list()
.bgb_m2_runtime$threads <- 1L
.bgb_m2_runtime$native_calls <- 0L
.bgb_m2_runtime$compatibility_fallbacks <- 0L
.bgb_m2_runtime$context_builds <- 0L
.bgb_m2_runtime$defer_model <- TRUE
.bgb_m2_runtime$prior_cache <- list()
.bgb_m2_clear_contexts <- function(){for(entry in .bgb_m2_runtime$contexts).Call(.bgb_c_bgb_tree_release,entry$handle);.bgb_m2_runtime$contexts<-list();invisible(NULL)}
.bgb_m2_deferred_Q <- function(areas_list,states_list,dmat,elist,amat=NULL,include_null_range=TRUE,normalize_TF=TRUE,makeCOO_TF=FALSE,min_precision=1e-26){
    a<-as.list(environment())
    if(!isTRUE(.bgb_m2_runtime$defer_model))return(do.call(rcpp_states_list_to_DEmat,a))
    structure(a,class='bgb_private_Q_parameters')
}
.bgb_m2_priors <- function(na,inp){
    key<-list(na,inp$maxent01s_param,inp$maxent01v_param,inp$maxent01j_param,inp$maxent01y_param,.bgb_precision())
    hits<-which(vapply(.bgb_m2_runtime$prior_cache,function(x)identical(x$key,key),logical(1)))
    if(length(hits))return(.bgb_m2_runtime$prior_cache[[hits[1]]]$value)
    saw_warning<-FALSE
    value<-withCallingHandlers({
    p<-list(relative_probabilities_of_subsets(na,inp$maxent01s_param,NA_val=0),relative_probabilities_of_vicariants(na,inp$maxent01v_param,NA_val=0),relative_probabilities_of_subsets(na,inp$maxent01j_param,NA_val=0),relative_probabilities_of_subsets(na,inp$maxent01y_param,NA_val=0))
    maxprob<-matrix(mapply(max,p[[1]],p[[2]],p[[3]],p[[4]],MoreArgs=list(na.rm=TRUE)),nrow=nrow(p[[1]]),ncol=ncol(p[[1]]));maxprob[maxprob>0]<-1;maxprob[maxprob<=0]<-0
    list(priors=p,maxmin=apply(maxprob,1,maxsize))
    },warning=function(w){saw_warning<<-TRUE})
    if(!saw_warning&&all(is.finite(unlist(value)))) .bgb_m2_runtime$prior_cache<-c(list(list(key=key,value=value)),head(.bgb_m2_runtime$prior_cache,7))
    value
}

.bgb_m2_clado <- function(x,n) {
    empty <- list(anc=numeric(),left=numeric(),right=numeric(),weight=numeric(),denom=numeric(),matrix=numeric(),kind=0L)
    if(!x$use_cpp || is.null(x$spPmat_inputs)) {
        if(!is.null(x$spPmat)) {empty$kind<-3L;empty$matrix<-as.numeric(x$spPmat)}
        return(empty)
    }
    inp<-x$spPmat_inputs
    inp$l[inp$l==c('_')]<-NULL;inp$l[inp$l==c('-')]<-NULL;inp$l[inp$l==c('-1')]<-NULL
    na<-max(lengths(inp$l));ones<-rep(1,n-as.integer(x$include_null_range))
    cached<-.bgb_m2_priors(na,inp);p<-cached$priors;maxmin<-cached$maxmin
    if(x$cppSpMethod==3&&is.null(x$m)){
        dmat_mat<-as.matrix(inp$dmat);storage.mode(dmat_mat)<-'double'
        actual_na<-nrow(dmat_mat)
        padded<-lapply(p,function(a){out<-matrix(0,actual_na,actual_na);nr<-min(nrow(a),actual_na);nc<-min(ncol(a),actual_na);out[seq_len(nr),seq_len(nc)]<-a[seq_len(nr),seq_len(nc),drop=FALSE];out})
        maxmin<-c(maxmin,rep(0,max(0,actual_na-length(maxmin))))
        empty$kind<-1L
        empty$model<-list(states=.bgb_state_buffers(inp$l),na=actual_na,parameters=c(inp$s,inp$v,inp$j,inp$y),dmat=dmat_mat,priors=padded,maxmin=maxmin)
        return(empty)
    }
    args<-list(Rcpp_leftprobs=ones,Rcpp_rightprobs=ones,l=inp$l,s=inp$s,v=inp$v,j=inp$j,y=inp$y,dmat=inp$dmat,maxent01s=p[[1]],maxent01v=p[[2]],maxent01j=p[[3]],maxent01y=p[[4]],max_minsize_as_function_of_ancsize=maxmin,printmat=FALSE)
    if(x$cppSpMethod==2) {
        if(!is.null(x$m))stop('method 2 with trait model requires original validation')
        coo<-do.call(rcpp_calc_anclikes_sp_COOprobs,args)
        ans<-list(anc=rep(seq_along(coo[[1]])-1L,lengths(coo[[1]])),left=unlist(coo[[1]],use.names=FALSE),right=unlist(coo[[2]],use.names=FALSE),weight=unlist(coo[[3]],use.names=FALSE),denom=vapply(coo[[3]],sum,numeric(1)),matrix=numeric(),kind=1L)
    } else {
        coo<-do.call(rcpp_calc_anclikes_sp_COOweights_faster,c(args,list(m=x$m,m_null_range=x$include_null_range,jts_matrix=x$jts_matrix)))
        ans<-list(anc=coo[[1]],left=coo[[2]],right=coo[[3]],weight=coo[[4]],denom=rcpp_calc_rowsums_for_COOweights_columnar(coo,printmat=FALSE),matrix=numeric(),kind=1L)
    }
    ans
}
.bgb_m2_calc <- function(frame) {
    x<-mget(names(formals(.bgb_m1_calc_loglike_sp_prebyte)),envir=frame,inherits=FALSE)
    contract<-.bgb_v12_runtime$contract
    root_spec<-.bgb_v12_runtime$root_prior
    original<-function(){if(contract!=0L||!is.null(root_spec))stop('UnsupportedModelContractError: this input is outside the validated native probability contract');.bgb_m2_runtime$compatibility_fallbacks<-.bgb_m2_runtime$compatibility_fallbacks+1L;if(inherits(x$Qmat,'bgb_private_Q_parameters'))x$Qmat<-do.call(rcpp_states_list_to_DEmat,unclass(x$Qmat));do.call(.bgb_m1_calc_loglike_sp_prebyte,x)}
    # Stratified treepieces retain their original orchestration until M3.
    if(isTRUE(x$stratified))return(original())
    boolean_args<-x[c('use_cpp','input_is_COO','sparse','include_null_range','calc_ancprobs','rootedge','stratified')]
    if(!all(vapply(boolean_args,function(a)is.logical(a)&&length(a)==1L&&!is.na(a),logical(1))))return(original())
    if(length(x$cppSpMethod)!=1L||!x$cppSpMethod %in% 1:3||length(x$printlevel)!=1L||is.na(x$printlevel)||isTRUE(x$input_is_COO)&&!isTRUE(x$sparse))return(original())
    # Keep historical validation/error diagnostics on exceptional inputs.
    valid<-is.matrix(x$tip_condlikes_of_data_on_each_state)&&inherits(x$phy,'phylo')&&length(x$return_what)==1L&&x$return_what %in% c('loglike','rootprobs','nodelikes','all')
    if(!valid||x$cppSpMethod==1||x$printlevel>=3||x$calc_ancprobs&&x$cppSpMethod!=3)return(original())
    if(!is.null(x$fixnode)&&length(x$fixnode)>1L&&(!identical(order(x$fixnode),seq_along(x$fixnode))||!is.matrix(x$fixlikes)||nrow(x$fixlikes)!=length(x$fixnode)))return(original())
    n<-ncol(x$tip_condlikes_of_data_on_each_state);nt<-length(x$phy$tip.label);nn<-nt+x$phy$Nnode
    if(!is.null(x$probs_of_states_at_root)&&length(x$probs_of_states_at_root)!=n)return(original())
    if(!is.null(x$states_allowed_TF)&&length(x$states_allowed_TF)!=n)return(original())
    prep_warnings<-character()
    prep<-tryCatch(withCallingHandlers(.bgb_m2_clado(x,n),warning=function(w){prep_warnings<<-c(prep_warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
    if(length(prep_warnings))return(original())
    if(inherits(prep,'error')){.bgb_m2_runtime$last_compatibility_reason<-conditionMessage(prep);return(original())}
    if(x$calc_ancprobs&&prep$kind!=1L)return(original())
    if(contract!=0L) {
        if(length(x$m)||length(x$fixnode)||isTRUE(x$rootedge)||!is.null(x$states_allowed_TF)||
           any(x$phy$edge.length<x$min_branchlength)||is.null(prep$model))
            stop('UnsupportedModelContractError: contract A/B currently requires a single-layer unmasked binary tree without hooks or fixed nodes')
        weights<-prep$model$parameters
        if(!.bgb_v12_standard_clado(weights))
            stop('UnsupportedModelContractError: only the six standard cladogenesis models are validated')
    }
    if(!is.null(root_spec)||contract!=0L) {
        states<-x$spPmat_inputs$l
        if(length(states)==n-1L&&isTRUE(x$include_null_range))states<-c(list(NA_integer_),states)
        resolved<-.bgb_v12_resolve_root(root_spec,states,contract)
        x$probs_of_states_at_root<-resolved$likelihood_weights
        .bgb_v12_runtime$last_root<-resolved
    }
    phy<-reorder(x$phy,'pruningwise')
    threads<-.bgb_m2_runtime$threads
    precision<-.bgb_precision()
    rootlen<-if(is.null(phy$root.edge))0 else phy$root.edge
    if(x$rootedge&&(length(rootlen)!=1L||!is.finite(rootlen)))return(original())
    config<-list(c(nt,n,nn,x$min_branchlength,x$include_null_range,x$rootedge,rootlen,prep$kind,x$cppSpMethod,length(x$m)),phy$edge-1L,phy$edge.length,x$tip_condlikes_of_data_on_each_state,as.numeric(x$fixnode-1L),as.numeric(x$fixlikes),as.numeric(x$states_allowed_TF),as.numeric(x$probs_of_states_at_root))
    deferred<-inherits(x$Qmat,'bgb_private_Q_parameters')
    qs<-if(deferred).bgb_state_buffers(x$Qmat$states_list)else list(numeric(),numeric())
    cs<-if(is.null(prep$model))list(numeric(),numeric())else prep$model$states
    config<-c(config,qs,cs,list(c(if(deferred)length(x$Qmat$areas_list)else 0,if(is.null(prep$model))0 else prep$model$na)))
    signature<-list(config=config,states=if(is.null(x$spPmat_inputs))NULL else x$spPmat_inputs$l,constraints=if(is.null(x$spPmat_inputs))NULL else x$spPmat_inputs[c('maxent01s_param','maxent01v_param','maxent01j_param','maxent01y_param')],trait=x$m,jts=x$jts_matrix,precision=precision,sparse=x$sparse,input_is_COO=x$input_is_COO,method=x$cppSpMethod,threads=threads,use_cpp=x$use_cpp)
    entries<-.bgb_m2_runtime$contexts
    hit<-which(vapply(entries,function(e)identical(e$signature,signature),logical(1)))
    if(length(hit)){entry<-entries[[hit[1]]];.bgb_m2_runtime$contexts<-c(list(entry),entries[-hit[1]])} else {
        entry<-list(signature=signature,config=config,handle=.Call(.bgb_c_bgb_tree_handle))
        if(length(entries)>=4L).Call(.bgb_c_bgb_tree_release,entries[[length(entries)]]$handle)
        entries<-c(list(entry),head(entries,3));.bgb_m2_runtime$contexts<-entries
        .bgb_m2_runtime$context_builds<-.bgb_m2_runtime$context_builds+1L
    }
    if(deferred){q<-rows<-cols<-numeric()}else if(x$input_is_COO) {
        if(!is.data.frame(x$Qmat)||ncol(x$Qmat)!=3L)return(original())
        q<-x$Qmat[[3]];rows<-x$Qmat[[1]];cols<-x$Qmat[[2]]
    } else {q<-as.numeric(x$Qmat);rows<-cols<-numeric()}
    full<-x$return_what %in% c('all','nodelikes')
    a<-list(c(x$input_is_COO,x$sparse,x$calc_ancprobs,full||prep$kind==3L),q,rows,cols,prep$anc,prep$left,prep$right,prep$weight,prep$denom,prep$matrix)
    qp<-if(deferred){
        qdmat<-as.matrix(x$Qmat$dmat);storage.mode(qdmat)<-'double'
        qamat<-if(is.null(x$Qmat$amat))matrix(0,nrow(qdmat),ncol(qdmat))else{m<-as.matrix(x$Qmat$amat);storage.mode(m)<-'double';m}
        list(c(1,x$Qmat$normalize_TF,x$Qmat$makeCOO_TF,x$Qmat$min_precision),qdmat,as.double(x$Qmat$elist),qamat)
    }else list(c(0,0,0,0),numeric(),numeric(),numeric())
    cp<-if(!is.null(prep$model))c(list(prep$model$parameters,prep$model$dmat),prep$model$priors,list(prep$model$maxmin,prep$model$na))else rep(list(numeric()),8)
    a<-c(a,qp,cp)
    result<-tryCatch(if(contract==0L&&is.null(root_spec)).Call(.bgb_c_bgb_tree,entry$handle,entry$config,a,precision,as.integer(threads)) else
        .Call(.bgb_c_bgb_tree_contract,entry$handle,entry$config,a,precision,as.integer(threads),contract,numeric()),error=function(e)e)
    if(inherits(result,'error')){if(contract!=0L||!is.null(root_spec))stop(result);.bgb_m2_runtime$last_compatibility_reason<-conditionMessage(result);return(original())}
    .bgb_v12_runtime$last_log_survival<-if(contract==1L)tail(result[[1]],1L) else NULL
    .bgb_m2_runtime$native_calls<-.bgb_m2_runtime$native_calls+1L
    num_short<-sum(phy$edge.length<x$min_branchlength)
    if(x$printlevel>=1L&&num_short>0){
        cat('\nRunning calc_loglike_sp():\nThis run of calc_loglike_sp() has a min_branchlength of: ',x$min_branchlength,'\n',sep='')
        cat('Branches shorter than this will be assumed to be connected to the tree with\nsympatric events (i.e., members of fossil lineages on ~0 length branches.)\n')
        cat('This tree has ',num_short,' branches < ',x$min_branchlength,'.\n\n',sep='')
    }
    # Original explicit-matrix/direct-product messages are unconditional.
    if(prep$kind==3L||prep$kind==0L&&x$use_cpp){
        cond<-if(prep$kind==3L)matrix(result[[4]],as.integer(result[[1]][2]),n)else NULL
        for(i in seq(1L,nrow(phy$edge),by=2L))if(!any(phy$edge.length[c(i,i+1L)]<x$min_branchlength)){
            if(prep$kind==3L){print('node_likelihood:');print(cond[phy$edge[i,1],])}
            else{print('WARNING #2: (use_cpp==TRUE && is.null(spPmat_inputs)==TRUE)');print('use_cpp=TRUE, no spPmat_inputs, direct multiplication')}
        }
    }
    if(x$printlevel>=1L&&x$sparse)cat(paste(rep('LR',x$phy$Nnode),collapse=''),'...segments done\n',sep='')
    if(x$calc_ancprobs){
        if(x$printlevel>=1L)cat('\nUppass starting for marginal ancestral states estimation!\n')
        if(x$sparse)cat(paste(rep('U',nrow(phy$edge)+as.integer(x$rootedge&&rootlen>0)),collapse=''))
        if(x$printlevel>=1L)cat('\nUppass completed for marginal ancestral states estimation!\n')
    }
    if(!full)return(result[[1]][1])
    nr<-as.integer(result[[1]][2]);down<-matrix(result[[3]],nr,n)
    if(nr>nn)rownames(down)<-c(rep('',nn),'relative_probs_of_each_state_at_bottom_of_root_branch')
    # These names deliberately keep the historical return-value semantics.
    if(x$return_what=='nodelikes')return(down)
    cond<-matrix(result[[4]],nr,n)
    if(nr>nn)rownames(cond)<-c(rep('',nn),'conditional_likelihoods_at_bottom_of_root_branch')
    out<-list(computed_likelihoods_at_each_node=result[[2]],relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS=down,condlikes_of_each_state=cond)
    if(x$calc_ancprobs){
        bottom<-matrix(result[[6]],nn,n);upbot<-matrix(result[[7]],nn,n);up<-matrix(result[[8]],nn,n)
        out$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS<-bottom
        out$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS<-upbot
        out$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS<-up
        marginal<-bottom*upbot;marginal<-marginal/rowSums(marginal)
        top<-down*up;root<-phy$edge[nrow(phy$edge),1];top[root,]<-down[root,];top<-top/rowSums(top)
        out$ML_marginal_prob_each_state_at_branch_bottom_below_node<-marginal
        out$ML_marginal_prob_each_state_at_branch_top_AT_node<-top
    }
    out$relative_probs_of_each_state_at_bottom_of_root_branch<-result[[5]]
    out$total_loglikelihood<-result[[1]][1]
    class(out)<-'calc_loglike_sp_results'
    out
}
calc_loglike_sp_prebyte <- .bgb_m1_calc_loglike_sp_prebyte
body(calc_loglike_sp_prebyte) <- quote(.bgb_m2_calc(environment()))
calc_loglike_sp <- compiler::cmpfun(calc_loglike_sp_prebyte)

# A run-local thread setting, restored even if the original optimizer errors.
.bgb_m1_bears_optim_run <- bears_optim_run
bears_optim_run <- .bgb_m1_bears_optim_run
body(bears_optim_run) <- quote({
    oldthreads<-.bgb_m2_runtime$threads
    on.exit({.bgb_m2_runtime$threads<-oldthreads},add=TRUE)
    hint<-BioGeoBEARS_run_object$num_cores_to_use
    .bgb_m2_runtime$threads<-if(is.numeric(hint)&&length(hint)==1L&&is.finite(hint)&&hint>=1)as.integer(hint)else 1L
    .bgb_m1_bears_optim_run(BioGeoBEARS_run_object,skip_optim,skip_optim_option)
})
