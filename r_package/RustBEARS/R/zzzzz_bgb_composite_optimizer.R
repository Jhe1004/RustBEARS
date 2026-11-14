# GPL >= 2. Reference-free deterministic bounded search; original result shells.
.bgb_composite_search <- function(par, fn, dots, lower, upper, maxit=500L) {
    if(is.null(names(par))) {
        table<-dots$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
        names(par)<-rownames(table)[table$type=="free"]
    }
    par <- stats::setNames(as.numeric(par), names(par)); keys <- names(par)
    lower <- rep_len(lower,length(par)); upper <- rep_len(upper,length(par))
    old <- options(RustBEARS.precision="StrictF64", BioGeoBEARS.precision="StrictF64")
    on.exit(options(old),add=TRUE)
    dims <- length(par); count <- 0L
    scale <- vapply(seq_along(par),function(i) if(keys[i]=="j") .05 else
        max(abs(par[i]),if(keys[i] %in% c("d","e","a")) .001 else .1),numeric(1))
    logs <- keys %in% c("d","e","a") & lower>0
    best <- new.env(parent=emptyenv()); best$value <- Inf; best$par <- par
    raw <- function(theta) {
        count <<- count+1L
        val <- tryCatch(do.call(fn,c(list(stats::setNames(theta,keys)),dots)),error=function(e) {
            if(grepl("invalid action input|negative likelihood|node_nan:|invalid likelihood row sum|stratified downpass NaN",conditionMessage(e))) return(NA_real_)
            stop(e)
        })
        if(length(val)!=1L || !is.finite(val) || val < -1e10) return(Inf)
        value <- -as.numeric(val)
        if(value<best$value) {best$value<-value;best$par<-stats::setNames(theta,keys)}
        value
    }
    encode <- function(theta,kind) {
        if(kind=="linear") return(theta/scale)
        theta[logs]<-log(pmax(theta[logs],lower[logs]));theta
    }
    decode <- function(z,kind) {
        if(kind=="linear") return(z*scale)
        z[logs]<-exp(z[logs]);z
    }
    objective <- function(z,kind) {
        val <- raw(pmin(upper,pmax(lower,decode(z,kind))))
        if(is.finite(val)) return(val)
        base <- if(is.finite(best$value)) best$value else 0
        base+1e6+1e5*sum((z-encode(best$par,kind))^2)
    }
    local <- function(theta,kind) {
        if(kind=="log")theta[logs]<-pmax(theta[logs],pmax(lower[logs],1e-8))
        low <- if(kind=="linear") lower/scale else encode(lower,"log")
        high <- if(kind=="linear") upper/scale else encode(upper,"log")
        result <- stats::optim(encode(theta,kind),function(z)objective(z,kind),method="L-BFGS-B",lower=low,upper=high,
            control=list(maxit=maxit,factr=450,pgtol=1e-7,ndeps=rep(1e-5,dims)))
        raw(pmin(upper,pmax(lower,decode(result$par,kind))))
        result
    }
    phy <- dots$phy
    height <- if(!is.null(phy)) max(ape::node.depth.edgelength(phy),1) else 1
    rc <- min(.05,max(1e-4,.15/height)); probes <- list(par)
    rates <- list(c(d=rc,e=rc,j=.01),c(d=.2*rc,e=1e-7,j=.02),c(d=.1*rc,e=.1*rc,j=1e-4))
    for(rate in rates) {
        x<-par;ks<-intersect(keys,names(rate));x[ks]<-rate[ks]
        x<-pmin(upper,pmax(lower,x));names(x)<-keys
        if(!any(vapply(probes,function(y)identical(as.numeric(x),as.numeric(y)),logical(1))))probes[[length(probes)+1L]]<-x
    }
    scores <- vapply(probes,raw,numeric(1)); order <- order(scores,seq_along(scores))
    order <- order[is.finite(scores[order])]
    if(!length(order))stop("no finite likelihood at the initial point or deterministic probes")
    chosen <- order[1]
    if(length(order)>1L) {
        others <- order[-1]; ei<-match("e",keys)
        opposite <- if(!is.na(ei)) others[vapply(others,function(i)(probes[[i]][ei]>1e-4)!=(probes[[order[1]]][ei]>1e-4),logical(1))] else integer()
        chosen<-c(chosen,if(length(opposite))opposite[1] else others[1])
    }
    searches <- lapply(chosen,function(i)local(probes[[i]],"linear"))
    searches[[length(searches)+1L]]<-local(best$par,"log")
    list(par=best$par,value=-best$value,count=count,convergence=if(any(vapply(searches,function(x)x$convergence==0L,logical(1))))0L else 1L,
         scale=rc,probes=probes,starts=chosen,searches=searches)
}

.bgb_composite_stats <- function(par,fn,...,method,lower,upper,control) {
    original_names<-names(par)
    if(method!="L-BFGS-B" || is.null(control$fnscale) || control$fnscale>=0 || !length(par))
        return(stats::optim(par,fn,...,method=method,lower=lower,upper=upper,control=control))
    dots<-list(...)
    if(is.null(names(par))) {
        table<-dots$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
        names(par)<-rownames(table)[table$type=="free"]
    }
    proposal<-.bgb_composite_search(par,fn,dots,lower,upper,if(is.null(control$maxit))500L else control$maxit)
    choices<-list(par,proposal$par)
    values<-vapply(choices,function(x)do.call(fn,c(list(x),dots)),numeric(1))
    settings<-control;settings$ndeps<-rep(1e-4,length(par));settings$factr<-450;settings$pgtol<-1e-7
    refined<-stats::optim(choices[[which.max(values)]],fn,...,method=method,lower=lower,upper=upper,control=settings)
    if(max(values)>refined$value) {refined$par<-choices[[which.max(values)]];refined$value<-max(values)}
    refined$message<-paste("StrictF64 proposal; deterministic multi-basin; target-precision refinement",refined$message)
    names(refined$par)<-original_names
    refined
}

.bgb_optimx <- function(par,fn,...,lower,upper,method,itnmax=NULL,hessian=NULL,control=list()) {
    original_names<-names(par)
    dots<-list(...);fn_dots<-dots;fn_dots$gr<-NULL;fn_dots$hess<-NULL
    if(is.null(names(par))) {
        table<-dots$BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
        names(par)<-rownames(table)[table$type=="free"]
    }
    call <- function(start) {
        names(start)<-original_names
        do.call(optimx::optimx,c(list(par=start,fn=fn,lower=lower,upper=upper,method=method,
                itnmax=itnmax,hessian=if(is.null(hessian))FALSE else hessian,control=control),dots))
    }
    if(!isTRUE(control$maximize) || !length(par))return(call(par))
    proposal<-.bgb_composite_search(par,fn,fn_dots,lower,upper,if(is.null(itnmax))500L else max(itnmax))
    choices<-list(par,proposal$par)
    values<-vapply(choices,function(x)do.call(fn,c(list(x),fn_dots)),numeric(1))
    result<-call(choices[[which.max(values)]])
    values_result<-as.numeric(result$value)
    best<-which.max(values_result)
    if(!length(best)||!is.finite(values_result[best]))best<-1L
    if(!is.finite(values_result[best]) || max(values)>values_result[best]) {
        selected<-choices[[which.max(values)]]
        result[best,seq_along(par)]<-as.list(unname(selected))
        result$value[best]<-max(values)
        result$convcode[best]<-proposal$convergence
    }
    attr(result,"composite_proposal")<-list(par=proposal$par,value=proposal$value,probes=proposal$probes,starts=proposal$starts,evaluations=proposal$count)
    result
}
