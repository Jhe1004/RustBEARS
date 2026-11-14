suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64')
folder<-tempfile('rustbears-v12-r-');dir.create(folder)
writeLines('((T1:1,T2:1):1,T3:2);',file.path(folder,'tree.newick'))
writeLines(c('3 2 (A B)','T1 10','T2 10','T3 01'),file.path(folder,'geog.data'))
make_fit<-function(contract,root=NULL) {
    x<-define_BioGeoBEARS_run(trfn=file.path(folder,'tree.newick'),geogfn=file.path(folder,'geog.data'),
        max_range_size=2,print_optim=FALSE,printlevel=-1,tmpwd=folder,calc_ancprobs=TRUE)
    x$likelihood_conditioning<-contract;x$root_prior<-root
    x$allow_null_tips<-contract=='unconditional_complete'
    for(k in c('d','e','j')) {
        value<-c(d=.2,e=.1,j=.1)[k]
        x$BioGeoBEARS_model_object@params_table[k,c('init','est')]<-value
    }
    x$BioGeoBEARS_model_object@params_table['j','type']<-'free'
    x<-readfiles_BioGeoBEARS_run(x);x$wd<-folder
    bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_all')
}
fits<-lapply(c('legacy_unconditioned','tip_survival','unconditional_complete'),make_fit)
names(fits)<-c('C','A','B')
stopifnot(all(vapply(fits,function(x)is.finite(x$total_loglikelihood),logical(1))))
stopifnot(abs(fits$C$total_loglikelihood-fits$B$total_loglikelihood-log(3))<1e-12)
stopifnot(abs(fits$A$total_loglikelihood-fits$B$total_loglikelihood+fits$A$log_tip_survival_prob)<1e-12)
prior<-define_root_prior('custom_weights',c(0,4,1,2))
custom<-make_fit('tip_survival',prior)
stopifnot(max(abs(custom$resolved_root_prior$sampling_probs-c(0,4,1,2)/7))<1e-15)
for(contract in c('A','B')) {
    sampler<-prepare_tree_sampler(fits[[contract]],'from_model')
    full<-simulate_datasets(sampler,n_reps=20L,seed=12345,return_node_states=TRUE)
    left<-simulate_datasets(sampler,n_reps=8L,seed=12345,start=1)
    right<-simulate_datasets(sampler,n_reps=12L,seed=12345,start=9,num_threads=2L)
    stopifnot(identical(full$tip_states,rbind(left$tip_states,right$tip_states)))
    if(contract=='A')stopifnot(all(full$tip_states>0))
    close_prepared_sampler(sampler);close_prepared_sampler(sampler)
    stopifnot(inherits(try(simulate_datasets(sampler,n_reps=1L),silent=TRUE),'try-error'))
}
stopifnot(inherits(try(simulate_datasets(fits$C,n_reps=1L),silent=TRUE),'try-error'))
sampler<-prepare_tree_sampler(fits$A)
maps<-draw_prepared_bsm(sampler,num_maps=20L,seed=12345)
stopifnot(nrow(maps$clado_events_df)==40L)
close_prepared_sampler(sampler)
empty_root_fit<-make_fit('unconditional_complete',define_root_prior('custom_weights',c(.3,.2,.2,.3)))
from_b<-simulate_datasets(empty_root_fit,n_reps=20L,conditioning='tip_survival')
stopifnot(all(from_b$tip_states>0))
writeLines(c('3 2 (A B)','T1 01','T2 10','T3 01'),file.path(folder,'geog.data'))
stopifnot(inherits(try(prepare_tree_sampler(fits$A),silent=TRUE),'try-error'))
writeLines(c('3 2 (A B)','T1 10','T2 10','T3 01'),file.path(folder,'geog.data'))
cat('R_V12_CONTRACTS_PASS\n')
