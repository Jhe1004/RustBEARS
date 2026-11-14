# Executable small example. Three bootstrap replicates are a smoke check only.
suppressPackageStartupMessages(library(RustBEARS))
options(RustBEARS.precision='StrictF64',RustBEARS.threads=1L)
work<-tempfile('rustbears-uncertainty-example-');dir.create(work)
writeLines('((a:1,b:1):1,(c:1,d:1):1);',file.path(work,'tree.nwk'))
writeLines(c('4 2 (A B)','a 10','b 10','c 01','d 01'),file.path(work,'ranges.data'))
run<-define_BioGeoBEARS_run(trfn=file.path(work,'tree.nwk'),geogfn=file.path(work,'ranges.data'),
    max_range_size=2,calc_ancprobs=TRUE,print_optim=FALSE,printlevel=-1,tmpwd=work)
run$likelihood_conditioning<-'tip_survival'
run$root_prior<-define_root_prior('uniform_nonnull')
run$num_cores_to_use<-1L;run$use_optimx<-FALSE
run<-readfiles_BioGeoBEARS_run(run);run$wd<-work
dec<-bears_optim_run(run)
runj<-run
runj$BioGeoBEARS_model_object@params_table['j','type']<-'free'
runj$BioGeoBEARS_model_object@params_table['j',c('init','est')]<-.01
decj<-bears_optim_run(runj)
intervals<-profile_ci(dec,params='e',n_grid=6L)
test<-bootstrap_lrt(dec,decj,n_boot=3L,seedval=42L,
    checkpoint_file=file.path(work,'bootstrap.rds'))
# AICc is undefined for a three-parameter model on four tips; explicitly use AIC.
average<-model_average(list(DEC=dec,DECJ=decj),criterion='AIC',num_maps=100L,seedval=42L)
precision<-bsm_diagnostics(average$joint_bsm,tree_target='fixed_sample')
simulated<-simulate_biogeog_datasets(dec,n_reps=3L,seedval=42L)
# Identical trees exercise the API; they do not represent tree uncertainty.
trees<-rep('((a:1,b:1):1,(c:1,d:1):1);',2L)
across_trees<-fit_over_trees(run,trees,models=c('DEC','DEC+J'),criterion='AIC',
    likelihood_conditioning='tip_survival',cold_start_audit_fraction=0,
    bsm_maps_per_tree=20L,checkpoint_file=file.path(work,'trees.rds'))
stopifnot(test$n_planned==3L,length(test$pval_failure_bounds)==2L,
    length(test$pval_mc_ci95)==2L,across_trees$n_trees_completed==2L,
    all(simulated$tip_states>0),abs(sum(average$weights)-1)<1e-12)
cat('R_UNCERTAINTY_EXAMPLE_PASS\n')
