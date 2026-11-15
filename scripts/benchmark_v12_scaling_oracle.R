# Original BioGeoBEARS 1.1.4 scalar workload. Matrix exports are preflight only.
args <- commandArgs(TRUE)
spec <- jsonlite::read_json(args[[1]], simplifyVector=TRUE)
.libPaths(c(spec$oracle_library, .libPaths()))
suppressPackageStartupMessages({library(ape);library(optimx);library(GenSA)
    library(rexpokit);library(cladoRcpp);library(snow);library(parallel);library(BioGeoBEARS)})
stopifnot(normalizePath(find.package('BioGeoBEARS')) == normalizePath(file.path(spec$oracle_library,'BioGeoBEARS')))
stopifnot(environmentName(environment(get('rcpp_states_list_to_DEmat',asNamespace('BioGeoBEARS')))) == 'cladoRcpp')
stopifnot(!'RustBEARS' %in% loadedNamespaces())
stopifnot(as.character(packageVersion('BioGeoBEARS')) == '1.1.4')
setwd(spec$workdir)
elapsed <- system.time({
    x <- define_BioGeoBEARS_run(trfn=spec$tree_file,geogfn=spec$geog_file,
        max_range_size=spec$max_range_size,print_optim=FALSE,printlevel=-1,
        tmpwd=spec$workdir,calc_ancprobs=FALSE,num_cores_to_use=1)
    x$wd <- spec$workdir
    for (name in c('d','e','j')) x$BioGeoBEARS_model_object@params_table[name,c('init','est')] <- spec$params[[name]]
    x <- readfiles_BioGeoBEARS_run(x)
    value <- bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_loglike')
})[['elapsed']]
result <- list(lnL=value,version=as.character(packageVersion('BioGeoBEARS')),R_version=R.version.string,
    library_path=find.package('BioGeoBEARS'),construction_evaluation_seconds=elapsed)
if (isTRUE(spec$export_components)) {
    mats <- get_Qmat_COOmat_from_BioGeoBEARS_run_object(x,BioGeoBEARS_model_object=x$BioGeoBEARS_model_object,
        max_range_size=spec$max_range_size,include_null_range=TRUE)
    n <- length(mats$states_list);coo <- mats$COO_weights_columnar;cc <- array(0,c(n,n,n))
    for(i in seq_along(coo[[4]]))cc[coo[[1]][i]+2L,coo[[2]][i]+2L,coo[[3]][i]+2L] <- coo[[4]][i]
    for(i in 2:n)cc[i,,] <- cc[i,,]/sum(cc[i,,])
    result$Q_colmajor <- as.numeric(mats$Qmat);result$C_colmajor <- as.numeric(cc)
    result$state_labels <- vapply(mats$states_list,function(s)if(anyNA(s)||!length(s)||any(s<0))'_' else paste(s,collapse=','),character(1))
    # Effective weights: the original empty-root contribution is exactly zero
    # with the shared C[0,:,:]=0 and an internal root. No score shift is applied.
    result$root_weights <- c(0,rep(1,n-1))
}
jsonlite::write_json(result,args[[2]],pretty=TRUE,auto_unbox=TRUE,digits=17,null='null')
