# Real-data regression: invalidate every point after a new nuisance-parameter basin.
suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64',RustBEARS.threads=1L)
directory<-tempfile('r-profile-rebase-');dir.create(directory)
x<-define_BioGeoBEARS_run(
    trfn=system.file('extdata/Psychotria_5.2_ultrametric.newick',package='RustBEARS'),
    geogfn=system.file('extdata/Psychotria_geog.data',package='RustBEARS'),
    max_range_size=4,print_optim=FALSE,printlevel=-1,tmpwd=directory,calc_ancprobs=TRUE)
x$likelihood_conditioning<-'tip_survival';x$num_cores_to_use<-1L;x$use_optimx<-FALSE
x$BioGeoBEARS_model_object<-RustBEARS:::.bgb_u_model_definition(x$BioGeoBEARS_model_object,'BAYAREALIKE+J')
for(k in c('d','e','j'))x$BioGeoBEARS_model_object@params_table[k,c('init','est')]<-
    c(d=1e-12,e=1e-12,j=.10811626692419736)[k]
x<-readfiles_BioGeoBEARS_run(x);x$wd<-directory
fit<-NULL;invisible(capture.output(fit<-bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_all')))
stopifnot(abs(fit$total_loglikelihood-(-24.260680130280434))<1e-10)
result<-profile_ci(fit,n_grid=8,tol=1e-4,num_threads=1)
table<-as.data.frame(result)
stopifnot(result$rebased_mle,abs(result$obj_mle-(-23.35568996585))<1e-6,
    all(vapply(result$intervals,function(i)i$endpoint_deviance_residual_max<1e-4,logical(1))))
print(table)
cat(sprintf('R_PROFILE_REBASE_PASS initial=%.17g updated=%.17g\n',fit$total_loglikelihood,result$obj_mle))
