# Concrete original BioGeoBEARS 1.1.4 adapter. Historical libraries are read-only.
args<-commandArgs(TRUE)
spec<-jsonlite::read_json(args[[1]],simplifyVector=TRUE)
output<-args[[2]]
.libPaths(c(spec$oracle_library,.libPaths()))
suppressPackageStartupMessages({library(ape);library(optimx);library(GenSA);library(rexpokit)
    library(cladoRcpp);library(snow);library(parallel);library(BioGeoBEARS)})
stopifnot(normalizePath(find.package('BioGeoBEARS'))==normalizePath(file.path(spec$oracle_library,'BioGeoBEARS')))
stopifnot(environmentName(environment(get('rcpp_states_list_to_DEmat',asNamespace('BioGeoBEARS'))))=='cladoRcpp')
stopifnot(!'RustBEARS' %in% loadedNamespaces())
dir.create(spec$workdir,recursive=TRUE,showWarnings=FALSE);setwd(spec$workdir)
x<-define_BioGeoBEARS_run(trfn=spec$tree_file,geogfn=spec$geog_file,max_range_size=spec$max_range_size,
    print_optim=FALSE,printlevel=-1,tmpwd=spec$workdir,calc_ancprobs=FALSE,num_cores_to_use=1)
x$wd<-spec$workdir;x$BioGeoBEARS_model_object@params_table['d',c('init','est')]<-spec$params$d
x$BioGeoBEARS_model_object@params_table['e',c('init','est')]<-spec$params$e
x$BioGeoBEARS_model_object@params_table['j',c('init','est')]<-0
x<-readfiles_BioGeoBEARS_run(x)
mats<-get_Qmat_COOmat_from_BioGeoBEARS_run_object(x,BioGeoBEARS_model_object=x$BioGeoBEARS_model_object,
    max_range_size=spec$max_range_size,include_null_range=TRUE)
n<-length(mats$states_list);coo<-mats$COO_weights_columnar;cc<-array(0,c(n,n,n))
for(i in seq_along(coo[[4]]))cc[coo[[1]][i]+2L,coo[[2]][i]+2L,coo[[3]][i]+2L]<-coo[[4]][i]
for(i in 2:n)cc[i,,]<-cc[i,,]/sum(cc[i,,])
labels<-vapply(mats$states_list,function(s)if(anyNA(s)||!length(s)||any(s<0))'_' else paste(s,collapse=','),character(1))
values<-numeric();optimized<-NULL
invisible(bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_loglike'))
elapsed<-system.time({
  for(i in seq_len(spec$repeats)) {
    if(spec$phase=='evaluate')value<-bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_loglike')
    else {
      optimized<-bears_optim_run(x,skip_optim=FALSE)
      value<-optimized$total_loglikelihood
    }
    values<-c(values,value)
  }
})[['elapsed']]
dlls<-getLoadedDLLs();paths<-vapply(dlls,function(d)d[['path']],character(1))
paths<-paths[grepl('cladoRcpp|rexpokit|libR',paths)]
result<-list(adapter='original_biogeobears',version=as.character(packageVersion('BioGeoBEARS')),
    R_version=R.version.string,library_path=find.package('BioGeoBEARS'),native_paths=paths,
    phase=spec$phase,repeats=spec$repeats,measured_api_seconds=elapsed,lnL=values,
    state_labels=labels,Q_colmajor=as.numeric(mats$Qmat),C_colmajor=as.numeric(cc),n_states=n,
    root_weights=c(0,rep(1,n-1)),contract='legacy_unconditioned',empty_clado='zero_row',
    optimized_params=if(is.null(optimized))NULL else as.list(optimized$outputs@params_table[c('d','e'),'est']),
    optimization_protocol=if(spec$phase=='evaluate')NULL else 'original bears_optim_run default; endpoint and work budget must be checked separately before any optimization speed comparison')
jsonlite::write_json(result,output,pretty=TRUE,auto_unbox=TRUE,digits=17,null='null')
