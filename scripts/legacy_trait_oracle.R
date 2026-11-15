#!/usr/bin/env Rscript
# Independent upstream evaluator. Only plain inputs cross process boundaries;
# upstream S4 classes and kernels are constructed in this clean R session.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('usage: legacy_trait_oracle.R input.rds output.rds')
if ('RustBEARS' %in% loadedNamespaces()) stop('oracle process already loaded RustBEARS')
plain <- function(x) {
    if (isS4(x) || is.environment(x) || is.function(x) || typeof(x) %in% c('externalptr','weakref','language')) return(FALSE)
    if (is.list(x)) return(all(vapply(x, plain, logical(1))))
    is.null(x) || is.atomic(x)
}
input <- readRDS(args[[1L]])
if (!is.list(input) || !plain(input)) stop('oracle input must contain only plain data, never cached S4/closures/native pointers')
if (length(input$oracle_library)) .libPaths(c(normalizePath(input$oracle_library, mustWork=TRUE), .libPaths()))
suppressPackageStartupMessages({
    library(ape); library(optimx); library(GenSA); library(rexpokit)
    library(cladoRcpp); library(snow); library(parallel); library(BioGeoBEARS)
})
if ('RustBEARS' %in% loadedNamespaces()) stop('upstream oracle unexpectedly loaded RustBEARS')
if (as.character(packageVersion('BioGeoBEARS')) != '1.1.4') stop('oracle requires the pinned BioGeoBEARS 1.1.4 source')
if (environmentName(environment(get('rcpp_states_list_to_DEmat',asNamespace('BioGeoBEARS')))) != 'cladoRcpp') stop('oracle does not resolve the original cladoRcpp kernel')
if (length(input$oracle_library) && normalizePath(find.package('BioGeoBEARS')) != normalizePath(file.path(input$oracle_library,'BioGeoBEARS'))) stop('wrong upstream library was loaded')
if (!requireNamespace('digest', quietly=TRUE)) stop('digest is required for oracle provenance')
hash <- function(path) digest::digest(file=path, algo='sha256')
configuration <- input$run_config
if (!is.list(configuration) || is.null(names(configuration)) || anyDuplicated(names(configuration))) stop('run_config must be a named list')
forbidden <- c('BioGeoBEARS_model_object','trait','likelihood_conditioning','root_prior','root_prior_spec','allow_null_tips')
if (any(names(configuration) %in% forbidden)) stop('run_config contains package-bound objects or nonlegacy conventions')
for (field in c('trfn','geogfn')) {
    path <- configuration[[field]]
    if (!is.character(path) || length(path)!=1L || !file.exists(path)) stop('missing input file: ',field)
    configuration[[field]] <- normalizePath(path,mustWork=TRUE)
}
parameters <- input$params_table
if (!is.data.frame(parameters) || !all(c('type','init','est','min','max') %in% names(parameters)) || anyDuplicated(rownames(parameters))) stop('invalid plain parameter table')
if (!is.list(input$candidates) || !length(input$candidates) || is.null(names(input$candidates)) || anyDuplicated(names(input$candidates))) stop('candidates must be a nonempty named list')
run <- BioGeoBEARS::define_BioGeoBEARS_run()
for (name in names(configuration)) run[name] <- configuration[name]
# Rebuild trait observations in the ORIGINAL namespace, not by deserializing a
# caller's class whose package field may route its methods back to RustBEARS.
if (length(input$trait_fn)) {
    run <- BioGeoBEARS::add_trait_to_BioGeoBEARS_run_object(run, trait_fn=normalizePath(input$trait_fn,mustWork=TRUE))
} else if (!is.null(input$trait_df)) {
    if (!is.data.frame(input$trait_df)) stop('trait_df must be a plain data.frame')
    run$trait <- BioGeoBEARS::define_tipranges_object(tmpdf=input$trait_df)
    if (is.null(run$trait_Pmat_txt)) stop('plain trait data require the explicit trait_Pmat_txt model')
}
definition <- BioGeoBEARS::define_BioGeoBEARS_model_object()
definition@params_table <- parameters
run$BioGeoBEARS_model_object <- definition
run$num_cores_to_use <- 1L
run$cluster_already_open <- FALSE
run$calc_ancprobs <- TRUE
run$print_optim <- FALSE
run$printlevel <- -1L
# This only selects the return_all container for an evaluation; no optimizer is
# run and no external package can silently improve the candidate parameters.
run$use_optimx <- FALSE
run <- BioGeoBEARS::readfiles_BioGeoBEARS_run(run)
free <- rownames(parameters)[parameters$type == 'free']
top_name <- 'relative_probs_of_each_state_at_branch_top_AT_node_UPPASS'
bottom_name <- 'relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS'
scores <- lapply(input$candidates, function(candidate) {
    if (!is.numeric(candidate) || any(!is.finite(candidate)) || anyDuplicated(names(candidate)) ||
        !setequal(names(candidate),free) || length(candidate)!=length(free)) stop('each candidate must supply every free parameter exactly once')
    current <- run
    p <- parameters
    p[free,'init'] <- candidate[free]
    p[free,'est'] <- candidate[free]
    current$BioGeoBEARS_model_object@params_table <- p
    current$BioGeoBEARS_model_object <- BioGeoBEARS::calc_linked_params_BioGeoBEARS_model_object(current$BioGeoBEARS_model_object,update_init=TRUE)
    started <- proc.time()[['elapsed']]
    result <- BioGeoBEARS::bears_optim_run(current,skip_optim=TRUE,skip_optim_option='return_all')
    elapsed <- proc.time()[['elapsed']] - started
    if (!is.finite(result$total_loglikelihood)) stop('nonfinite upstream score')
    if (!isTRUE(all.equal(unname(result$outputs@params_table[free,'est']),unname(candidate[free]),tolerance=0))) stop('upstream evaluation changed candidate parameters')
    if (!is.matrix(result[[top_name]]) || !is.matrix(result[[bottom_name]])) stop('upstream did not return ancestral matrices')
    list(lnL=unname(result$total_loglikelihood),ancestral_top=result[[top_name]],
         ancestral_bottom=result[[bottom_name]],params_table=result$outputs@params_table,
         elapsed_seconds=elapsed,candidate_within_declared_bounds=all(candidate[free]>=p[free,'min'] & candidate[free]<=p[free,'max']))
})
packages <- c('BioGeoBEARS','cladoRcpp','rexpokit','ape','optimx','GenSA')
metadata <- list(schema='rustbears.independent_legacy_oracle.v1',
    input_sha256=hash(args[[1L]]),script_sha256=hash(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[[1L]])),
    R=R.version.string,loaded_rustbears='RustBEARS' %in% loadedNamespaces(),
    upstream_remote_sha=packageDescription('BioGeoBEARS')$RemoteSha,
    version=setNames(lapply(packages,function(p)as.character(packageVersion(p))),packages),
    package_path=setNames(lapply(packages,find.package),packages),
    native_files=lapply(c('cladoRcpp','rexpokit'),function(p) {
        files <- list.files(file.path(find.package(p),'libs'),recursive=TRUE,full.names=TRUE)
        files <- files[!dir.exists(files)]
        setNames(lapply(files,hash),files)
    }),evaluation='upstream native kernels, fixed supplied parameters, no optimization')
if (metadata$loaded_rustbears) stop('oracle loaded RustBEARS during evaluation')
saveRDS(list(scores=scores,metadata=metadata),args[[2L]],version=2)
