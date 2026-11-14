suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64')
# Source in an isolated child environment so development tests never unlock or
# mutate a loaded package namespace. Installed export checks run separately.
scope<-new.env(parent=asNamespace('RustBEARS'))
sys.source('r_package/RustBEARS/R/zzzzzzzz_bgb_uncertainty.R',envir=scope)
sys.source('r_package/RustBEARS/R/zzzzzzzzz_bgb_profile_bootstrap.R',envir=scope)
evalq({
folder<-tempfile('rustbears-statistics-');dir.create(folder);on.exit<-NULL
make_fit<-function(tree,geography,contract='unconditional_complete',free='e',bounds=c(0,2),root=NULL) {
    prefix<-tempfile('case-',tmpdir=folder);treefn<-paste0(prefix,'.tree');geogfn<-paste0(prefix,'.data')
    writeLines(tree,treefn);writeLines(geography,geogfn)
    x<-define_BioGeoBEARS_run(trfn=treefn,geogfn=geogfn,max_range_size=if(grepl(' 2 ',geography[1],fixed=TRUE))2 else 1,
        print_optim=FALSE,printlevel=-1,tmpwd=folder,calc_ancprobs=TRUE)
    x$likelihood_conditioning<-contract;x$allow_null_tips<-contract=='unconditional_complete';x$root_prior<-root
    p<-x$BioGeoBEARS_model_object@params_table
    for(k in c('d','e','j')){p[k,'type']<-'fixed';p[k,c('init','est')]<-0}
    if(grepl(' 2 ',geography[1],fixed=TRUE))p['e',c('init','est')]<-.2
    for(k in free){p[k,'type']<-'free';p[k,c('min','max')]<-bounds;p[k,c('init','est')]<-.1}
    x$BioGeoBEARS_model_object@params_table<-p;x<-readfiles_BioGeoBEARS_run(x);x$wd<-folder
    out<-NULL;invisible(capture.output(out<-bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_all')))
    .bgb_s_fit(out)
}
one<-c('3 1 (X)','A 1','B 1','C 1')
fit<-make_fit('((A:1,B:1):3,C:4);',one)
profile<-profile_ci(fit,params='e',n_grid=5,include_onesided_zero_ref=TRUE)
i<-profile$intervals$e
stopifnot(i$lower==0,i$status_lower=='at_physical_zero',abs(i$upper-stats::qchisq(.95,1)/18)<1e-5,
    abs(i$onesided_upper_from_zero-stats::qchisq(.9,1)/18)<1e-5,i$endpoint_deviance_residual_max<1e-4)
legacy<-make_fit('((A:1,B:1):3,C:4);',one,'legacy_unconditioned')
stopifnot(profile_ci(legacy,'e',n_grid=4)$interval_kind=='profile_support_interval')
complete_empty<-make_fit('(A:1,B:1);',c('2 1 (X)','A 0','B 1'))
legacy_empty<-complete_empty;legacy_empty$inputs$likelihood_conditioning<-'legacy_unconditioned'
legacy_empty<-.bgb_u_refit(legacy_empty,optimize=FALSE)
stopifnot(inherits(try(prepare_tree_sampler(legacy_empty),silent=TRUE),'try-error'))
sampler_empty<-prepare_tree_sampler(complete_empty)
stopifnot(nrow(draw_prepared_bsm(sampler_empty,2L)$clado_events_df)==2L)
close_prepared_sampler(sampler_empty)
small<-fit;small$outputs@params_table['e','max']<-.01
stopifnot(is.na(profile_ci(small,'e',n_grid=4)$intervals$e$upper),
          profile_ci(small,'e',n_grid=4)$intervals$e$status_upper=='censored_at_box_upper')

two<-c('2 2 (X Y)','A 10','B 01')
structural<-make_fit('(A:1,B:1);',two,free='d',root=define_root_prior('custom_weights',c(0,1,0,0)))
zero<-.bgb_s_point(structural,list(d=0))
stopifnot(is.infinite(zero$profile_obj),zero$profile_obj<0,zero$structural_zero)
sp<-profile_ci(structural,'d',n_grid=5)
stopifnot(sp$intervals$d$lower>0,sp$intervals$d$status_lower=='infeasible_zero_approached',
          sp$intervals$d$endpoint_deviance_residual_max<1e-4)
founder<-make_fit('(A:1,B:1);',two,free='j',root=define_root_prior('custom_weights',c(0,1,0,0)))
founder_point<-.bgb_s_point(founder,list(j=.2))
stopifnot(abs(founder_point$params[['ysv']]-2.8)<1e-14,abs(founder_point$params[['y']]-2.8/3)<1e-14)

# A deliberately suboptimal finite result must be rebased before intervals.
stale<-.bgb_u_refit(fit,start=c(e=.1),optimize=FALSE)
rebased<-profile_ci(stale,'e',n_grid=4)
stopifnot(rebased$rebased_mle,rebased$obj_mle>stale$total_loglikelihood+1e-6)

# Canonical observation identity accepts rotations and data-row permutations.
rotated<-fit;rotated$input_snapshot$tree_text<-'(C:4,(B:1,A:1):3);'
rotated$input_snapshot$geography_text<-c(one[1],one[4:2])
stopifnot(identical(.bgb_u_data_signature(fit),.bgb_u_data_signature(rotated)))
.bgb_u_compatible(list(fit,rotated))

bounds<-.bgb_s_failure_intervals(2,3,10)
stopifnot(identical(bounds$failure_bounds,c(3/11,6/11)),
          abs(bounds$mc_ci95[1]-stats::qbeta(.025,2,9))<1e-14,abs(bounds$mc_ci95[2]-stats::qbeta(.975,6,5))<1e-14)
all_failed<-.bgb_s_failure_intervals(0,10,10)
stopifnot(identical(all_failed$mc_ci95,c(0,1)),all_failed$failure_bounds[2]==1)

# Shared Python/R SHA-256 vectors, including Unicode canonical escaping.
vector<-.bgb_u_task_seed(12345,'bootstrap_null','DEC','1',outer_rep=1,role='simulate')
stopifnot(vector$fingerprint=='f9aa047750e7cb9f46518889bfba1ef2c54f399782a2049ea1e4f076404d8806')
unicode<-.bgb_u_task_seed(12345,'审计','模型一','树二')
stopifnot(unicode$fingerprint=='dd61e9dda555282ef34e596bfd1ab21963295e562b12f80f5169bfae18be786c')
task<-list(purpose='test',model_id='DEC',tree_id='1')
stopifnot(inherits(try(.bgb_s_tasks(12345,list(task,task)),silent=TRUE),'try-error'))

null<-.bgb_u_refit(fit,fixed=list(e=0),optimize=FALSE)
path<-file.path(folder,'bootstrap.rds')
boot<-bootstrap_lrt(null,fit,n_boot=3,n_alt_sim=2,calc_surrogate_power=TRUE,checkpoint_file=path)
stopifnot(boot$n_planned==3,boot$n_converged+boot$n_failed==3,boot$pval_mc_conservative==1,
          length(boot$replicate_artifacts)==5,file.exists(paste0(path,'.manifest.json')))
file.copy(paste0(path,'.manifest.json'),'/tmp/r_v12_checkpoint_manifest.json',overwrite=TRUE)
continued<-bootstrap_lrt(null,fit,n_boot=3,n_alt_sim=2,calc_surrogate_power=TRUE,checkpoint_file=path,resume=TRUE)
stopifnot(identical(boot$replicates_df,continued$replicates_df),identical(boot$surrogate_replicates_df,continued$surrogate_replicates_df))
stopifnot(inherits(try(bootstrap_lrt(null,fit,n_boot=4,checkpoint_file=path,resume=TRUE),silent=TRUE),'try-error'))
stopifnot(all(vapply(boot$replicate_artifacts,function(r)r$status!='converged'||
    !is.null(r$fits$alt$ML_marginal_prob_each_state_at_branch_top_AT_node),logical(1))))

plugin<-profile_ci(fit,'e',method='bootstrap_plugin',n_boot_plugin=3,n_grid=4)
stopifnot(nrow(plugin$intervals$e$bootstrap_plugin_df)==3,plugin$method=='bootstrap_plugin')

surface_fit<-make_fit('((A:1,B:1):3,C:4);',one,free=c('d','e'))
surface<-profile_surface(surface_fit,'d','e',grid_x=c(0,.1),grid_y=c(0,.2),reoptimize_nuisance=FALSE)
stopifnot(nrow(surface)==4,max(abs(surface$profile_obj+9*surface$e))<1e-12)

floor_alt<-fit;floor_alt$outputs@params_table['e','min']<-1e-12
floor_boot<-bootstrap_lrt(null,floor_alt,n_boot=1)
stopifnot(floor_boot$observed_alt_fit$outputs@params_table['e','min']==0)
custom<-fit;custom$outputs@params_table['e','min']<-.001
stopifnot(inherits(try(bootstrap_lrt(null,custom,n_boot=1),silent=TRUE),'try-error'))
legacy_null<-.bgb_u_refit(legacy,fixed=list(e=0),optimize=FALSE)
stopifnot(inherits(try(bootstrap_lrt(legacy_null,legacy,n_boot=1,target_conditioning='unconditional_complete'),silent=TRUE),'try-error'))
changed<-bootstrap_lrt(legacy_null,legacy,n_boot=1,target_conditioning='unconditional_complete',refit_observed=TRUE)
stopifnot(changed$recomputed_observed_under_new_contract,changed$likelihood_conditioning=='unconditional_complete')

# One deliberately unrecoverable replicate remains in the planned denominator.
pair_original<-.bgb_s_pair;injections<-0L
.bgb_s_pair<-function(null,alt,num_threads,maxit=500L,geogfn=NULL,audit=NULL) {
    if(!is.null(geogfn)){injections<<-injections+1L;if(injections<=2L)stop('injected replicate failure')}
    pair_original(null,alt,num_threads,maxit,geogfn,audit)
}
failed<-bootstrap_lrt(null,fit,n_boot=3)
.bgb_s_pair<-pair_original
stopifnot(failed$n_failed==1,failed$n_converged==2,nrow(failed$replicates_df)==3,
          identical(failed$pval_failure_bounds,c(3/4,1)),grepl('injected',failed$replicates_df$error[1]))
# Signed failed statistics remain diagnostic facts, never inferential observations.
fit_original<-.bgb_s_fit
.bgb_s_fit<-function(res,fixed=list(),start=.bgb_u_params(res),geogfn=NULL,num_threads=1L,maxit=500L,ancestors=TRUE) {
    value<-fit_original(res,fixed,start,geogfn,num_threads,maxit,ancestors)
    if(!is.null(geogfn))value$total_loglikelihood<-if(res$outputs@params_table['e','type']=='free')-2 else -1
    value
}
raw_path<-file.path(folder,'failed_raw.rds')
raw_failed<-bootstrap_lrt(null,fit,n_boot=1,checkpoint_file=raw_path)
.bgb_s_fit<-fit_original
record<-raw_failed$replicate_artifacts[['null:1']]
stopifnot(raw_failed$n_failed==1L,is.na(raw_failed$replicates_df$d_stat),record$raw_obj_null== -1,
    record$raw_obj_alt== -2,record$raw_d_stat== -2,length(record$attempts)==2L,
    all(vapply(record$attempts,function(a)a$status=='failed'&&a$raw_d_stat== -2&&length(a$alternative_fits)==2L,logical(1))),
    identical(raw_failed$pval_mc_ci95,c(0,1)),identical(raw_failed$pval_failure_bounds,c(.5,1)),
    readRDS(raw_path)$records[['null:1']]$raw_d_stat== -2)
# The actual optimizer-failure path attaches raw values before throwing.
search_original<-.bgb_composite_search
.bgb_composite_search<-function(theta,...)list(par=theta,convergence=1L,value=-3,count=7L)
search_error<-tryCatch(.bgb_s_fit(fit),error=function(e)e)
.bgb_composite_search<-search_original
stopifnot(inherits(search_error,'rustbears_fit_failure'),search_error$fit_snapshot$obj== -3,
    !search_error$fit_snapshot$converged,search_error$fit_snapshot$n_evals==7L,
    is.finite(search_error$fit_snapshot$params[['e']]))
cat('R_V12_STATISTICS_PASS\n')
},scope)
