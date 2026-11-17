# Compare the installed native R frontend with deterministic Python references.
suppressPackageStartupMessages(library(RustBEARS))
args<-commandArgs(TRUE)
if(length(args)!=1L)stop('usage: Rscript scripts/validate_v12_r_cross.R REFERENCE_DIRECTORY')
reference<-args[1]
env<-new.env(parent=globalenv())
invisible(capture.output(sys.source('tests/test_r_v12_summaries.R',envir=env)))
read_reference<-function(name)read.delim(file.path(reference,name),check.names=FALSE,stringsAsFactors=FALSE,
    fileEncoding='UTF-8',na.strings=c('','NA'))
equal_frames<-function(actual,expected,tolerance=2e-12) {
    stopifnot(nrow(actual)==nrow(expected),all(names(expected) %in% names(actual)))
    for(key in names(expected)) {
        a<-actual[[key]];b<-expected[[key]]
        if(is.numeric(b)) {
            if(!identical(is.na(a),is.na(b))||any(abs(a-b)>tolerance,na.rm=TRUE))stop(paste('numeric mismatch:',key))
        } else if(is.logical(a)) {
            if(!identical(tolower(as.character(a)),tolower(as.character(b))))stop(paste('logical mismatch:',key))
        } else if(!identical(as.character(a),as.character(b)))stop(paste('text mismatch:',key))
    }
}
hashes<-read_reference('topology_hashes.tsv')
for(i in seq_len(nrow(hashes)))stopifnot(identical(RustBEARS:::.bgb_u_topology_hash(ape::read.tree(text=hashes$tree[i])),hashes$sha256[i]))
expected<-read_reference('six_models.tsv');maximum<-0
for(family in unique(expected$model)) {
    definition<-RustBEARS:::.bgb_u_model_definition(env$x$BioGeoBEARS_model_object,family)
    values<-c(d=.2,e=.1,j=if(grepl('+J',family,fixed=TRUE)).1 else 0)
    for(k in names(values)) {
        definition@params_table[k,c('init','est')]<-values[k]
        definition@params_table[k,'type']<-'fixed'
    }
    directory<-tempfile('r-cross-family-');dir.create(directory)
    fit<-RustBEARS:::.bgb_u_fit_tree(env$x,ape::read.tree(text=env$tree),definition,directory,optimize=FALSE)
    members<-RustBEARS:::.bgb_u_tree_members(ape::read.tree(text=env$tree))
    clades<-vapply(members,function(x)paste(sort(x),collapse='|'),character(1))
    subset<-expected[expected$model==family,]
    actual<-fit$ML_marginal_prob_each_state_at_branch_top_AT_node[cbind(match(subset$clade,clades),subset$state+1L)]
    error<-max(abs(actual-subset$probability),abs(fit$total_loglikelihood-subset$lnL))
    stopifnot(error<2e-12);maximum<-max(maximum,error)
}
equal_frames(env$avg$joint_bsm$clado_events_df,read_reference('average_clado.tsv'))
equal_frames(env$avg$joint_bsm$ana_events_df,read_reference('average_ana.tsv'))
equal_frames(env$avg$joint_bsm$map_metadata_df,read_reference('average_metadata.tsv'))
diagnostic<-env$diag$event_summary_df
diagnostic$event_type<-diagnostic$statistic
diagnostic$assessed<-diagnostic$assessment=='assessed'
expected<-read_reference('average_diagnostics.tsv')
equal_frames(diagnostic,expected[,intersect(names(expected),names(diagnostic)),drop=FALSE])
equal_frames(env$multi$joint_bsm$clado_events_df,read_reference('multitree_clado.tsv'))
equal_frames(env$multi$joint_bsm$ana_events_df,read_reference('multitree_ana.tsv'))
equal_frames(env$multi$joint_bsm$map_metadata_df,read_reference('multitree_metadata.tsv'))
for(target in c('fixed_sample','posterior_population')) {
    diagnostic<-bsm_diagnostics(env$multi$joint_bsm,min_mean_count=0,tree_target=target,
        tree_dependence=if(target=='posterior_population')'independent' else NULL)$event_summary_df
    diagnostic$event_type<-diagnostic$statistic;diagnostic$assessed<-diagnostic$assessment=='assessed'
    expected<-read_reference(paste0('multitree_',target,'_diagnostics.tsv'))
    equal_frames(diagnostic,expected[,intersect(names(expected),names(diagnostic)),drop=FALSE])
}
chain<-fit_over_trees(env$x,c(rep(c(env$tree,env$rotated),30L),env$tree),models=list(fixed=env$x$BioGeoBEARS_model_object),
    likelihood_conditioning='tip_survival',criterion='AIC',init_strategy='full_cold_start',cold_start_audit_fraction=0,
    bsm_maps_per_tree=4L,seedval=54321L)
equal_frames(chain$joint_bsm$clado_events_df,read_reference('chain_clado.tsv'))
equal_frames(chain$joint_bsm$ana_events_df,read_reference('chain_ana.tsv'))
equal_frames(chain$joint_bsm$map_metadata_df,read_reference('chain_metadata.tsv'))
chain_diag<-bsm_diagnostics(chain$joint_bsm,min_mean_count=0,tree_target='posterior_population',tree_dependence='mcmc_chain')$event_summary_df
chain_diag$event_type<-chain_diag$statistic;chain_diag$assessed<-chain_diag$assessment=='assessed'
expected_chain<-read_reference('chain_mcmc_diagnostics.tsv')
equal_frames(chain_diag,expected_chain[,intersect(names(expected_chain),names(chain_diag)),drop=FALSE])
cat(sprintf('R_PYTHON_CROSS_PASS six_models_max_abs_error=%.17g; exact task hashes and all 40 model-stratified plus 80 multitree event histories\n',maximum))
