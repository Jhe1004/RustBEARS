args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L)
branch <- normalizePath(args[1]); original <- normalizePath(args[2])
out <- normalizePath(args[3]); rscript <- normalizePath(args[4])
library(RustBEARS)
stopifnot(!interactive())
ns <- asNamespace('RustBEARS')
# The retained APE-5 coordinate helper uses these unqualified APE internals.
# Supply its declared legacy calling context in this short-lived R process;
# do not replace plotting functions, native symbols, or package namespaces.
ape_bindings <- c('node_height','node_depth','node_depth_edgelength',
                  'node_height_clado','.PlotPhyloEnv')
for (name in ape_bindings) assign(name,getFromNamespace(name,'ape'),envir=.GlobalEnv)
plotenv <- new.env(parent=ns)
targets <- c('plot_simDEC_DECJ_inferences','plot_inf_accuracy_vs_SSEsims_v2',
             'plot_inf_accuracy_vs_SSEsims_v3','plot_SSEsims_treestats_pt1','plot_BSM_fns')
source_files <- file.path(branch,'r_package/RustBEARS/R',
                         c('BioGeoBEARS_SSEsim_makePlots_v1.R','BioGeoBEARS_stochastic_mapping_v1.R'))
original_defs <- list()
for (p in source_files) for (expr in parse(p)) {
  if (is.call(expr) && as.character(expr[[1]]) %in% c('<-','=') &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% targets) {
    name <- as.character(expr[[2]])
    eval(expr,envir=plotenv)
    original_defs[[name]] <- list(body=body(plotenv[[name]]),formals=formals(plotenv[[name]]))
  }
}
stopifnot(setequal(names(original_defs),targets))

# Reuse saved three-tip fits and one saved stratified map; no likelihood or
# simulation routine is called. Copy input files to this new validation folder.
inputdir <- file.path(out,'inputs');dir.create(inputdir)
for (f in c('tree.newick','geog.data')) {
  stopifnot(file.copy(file.path(original,f),file.path(inputdir,f),overwrite=FALSE))
}
read_fit <- function(f) {
  e<-new.env();load(file.path(original,f),envir=e);value<-e$res
  value$inputs$trfn<-file.path(inputdir,'tree.newick')
  value$inputs$geogfn<-file.path(inputdir,'geog.data')
  value
}
dec<-read_fit('3taxa_DEC_M3_areasallowed_v1.Rdata')
decj<-read_fit('3taxa_DEC+J_M3_areasallowed_v1.Rdata')
# This older SSE display helper expects the pre-GenSA optim container sign.
# Normalize only the engineering display copies to their already stored LnL;
# retain the original GenSA objects and the direct-call failure separately.
as_legacy_optim_display <- function(res) {
  stopifnot(is.finite(res$total_loglikelihood),
            abs(res$optim_result$value+res$total_loglikelihood)<1e-10)
  res$inputs$use_optimx<-FALSE
  res$optim_result$value<-res$total_loglikelihood
  res
}
dec_display<-as_legacy_optim_display(dec)
decj_display<-as_legacy_optim_display(decj)
color_checks<-lapply(c(TRUE,FALSE),function(include) {
  states<-get('rcpp_areas_list_to_states_list',ns)(c('A','B','C'),maxareas=3,include_null_range=include)
  colors<-get_colors_for_states_list_0based(c('A','B','C'),states,max_range_size=3,
                                           plot_null_range=include)
  stopifnot(length(colors)==length(states),length(colors)==if(include)8L else 7L)
  list(include_null_range=include,state_count=length(states),color_count=length(colors))
})
cl<-new.env();load(file.path(original,'RES_clado_events_tables.Rdata'),envir=cl)
an<-new.env();load(file.path(original,'RES_ana_events_tables.Rdata'),envir=an)
BSM_results<-list(master_table_cladogenetic_events=cl$RES_clado_events_tables[[1]],
                  table_w_anagenetic_events=an$RES_ana_events_tables[[1]])
save(BSM_results,file=file.path(inputdir,'one_saved_map.Rdata'))
# A second engineering display uses the same nonempty history with a seven-state
# numbering convention. Drop only the zero-mass empty column; reindex all saved
# state numbers without changing event labels/times or any nonempty probability.
decj_no_empty<-decj
decj_no_empty$inputs$include_null_range<-FALSE
for(n in c('ML_marginal_prob_each_state_at_branch_top_AT_node',
            'ML_marginal_prob_each_state_at_branch_bottom_below_node')) {
  stopifnot(ncol(decj_no_empty[[n]])==8L,all(decj_no_empty[[n]][,1]==0))
  decj_no_empty[[n]]<-decj_no_empty[[n]][,-1,drop=FALSE]
}
stopifnot(is.null(decj_no_empty$inputs$states_list),is.null(decj_no_empty$inputs$states_list_0based))
shift_nonempty<-function(x) {
  stopifnot(all(is.na(x)|(x>=2 & x<=8)))
  x-1L
}
no_empty_map<-BSM_results
for(n in c('sampled_states_AT_nodes','sampled_states_AT_brbots','samp_LEFT_dcorner','samp_RIGHT_dcorner'))
  no_empty_map$master_table_cladogenetic_events[[n]]<-
    shift_nonempty(no_empty_map$master_table_cladogenetic_events[[n]])
for(n in c('current_rangenum_1based','new_rangenum_1based'))
  no_empty_map$table_w_anagenetic_events[[n]]<-shift_nonempty(no_empty_map$table_w_anagenetic_events[[n]])
texts<-no_empty_map$master_table_cladogenetic_events$anagenetic_events_txt_below_node
for(i in seq_along(texts)) if(!is.na(texts[i]) && texts[i]!='none') {
  ev<-events_txt_into_events_table(texts[i])
  for(n in c('current_rangenum_1based','new_rangenum_1based')) ev[[n]]<-shift_nonempty(as.numeric(ev[[n]]))
  texts[i]<-events_table_into_txt(ev)
}
no_empty_map$master_table_cladogenetic_events$anagenetic_events_txt_below_node<-texts
BSM_results_original<-BSM_results;BSM_results<-no_empty_map
save(BSM_results,file=file.path(inputdir,'one_saved_map_no_empty_display.Rdata'))
BSM_results<-BSM_results_original

# Deterministic engineering summary rows exercise every drawing panel in the
# three SSE-summary helpers. These are neither simulations nor study results.
tab<-data.frame(params_txt=rep(c('fixture_A','fixture_B'),each=3),
                d=rep(c(.05,.1),each=3),e=.02,j=rep(c(0,.03),each=3),
                d_DECinf=c(.04,.05,.06,.09,.1,.11),
                d_DECJinf=c(.045,.05,.055,.095,.1,.105),
                e_DECinf=c(.01,.02,.03,.015,.02,.025),
                e_DECJinf=c(.015,.02,.025,.015,.02,.025),
                j_DECinf=0,j_DECJinf=c(0,.005,.01,.025,.03,.035),
                LnL_DEC=c(-20,-21,-22,-23,-24,-25),
                LnL_DECJ=c(-19.8,-20.7,-21.5,-20,-21,-22),
                brate=.3,drate=.3,b_exp=1,d_exp=-1)
accuracy<-matrix(rep(c(0,1,1,0,1,1),147),nrow=6,ncol=147)
accuracy[,1]<-c(10,12,14,16,18,20)
event_names<-c('num_fossil_tips',outer(c('y','s','v','j','d','e','a'),
                                    c('actual','observed'),paste,sep='_'))
events<-as.data.frame(matrix(rep(1:6,length(event_names)),nrow=6))
names(events)<-event_names
plotenv$uniq_params_num<-2L
plotenv$siminf_stats_subset<-tab
plotenv$ancstate_accuracies_DEC_subset<-accuracy
plotenv$ancstate_accuracies_DECJ_subset<-accuracy
plot_inputs<-list(titletxt1='Engineering plotting fixture; no study inference',
                 doPDF=TRUE,siminf_stats_good=tab,simtype='ENGINEERING FIXTURE',
                 brate_plot=.3,drate_plot=.3,b_exp_plot=1,d_exp_plot=-1,
                 uniq_params_txt=c('fixture_A','fixture_B'),segwidth_truth=.05,
                 ancstate_accuracies_DEC_good=accuracy,
                 ancstate_accuracies_DECJ_good=accuracy,simevent_counts_good=events,
                 sim_tipstates_good=matrix(rep(c(2,3,4,5,6,7),50),nrow=6),
                 treeheights_incl_fossil_good=c(20,22,24,26,28,30))
saveRDS(list(plot_inputs=plot_inputs,implicit_globals=list(uniq_params_num=2L,
            siminf_stats_subset=tab,ancstate_accuracies_DEC_subset=accuracy,
            ancstate_accuracies_DECJ_subset=accuracy)),file.path(inputdir,'engineering_summary_inputs.rds'))
write.csv(tab,file.path(inputdir,'engineering_summary_rows.csv'),row.names=FALSE)

sentinel<-file.path(out,'viewer_sentinel.R')
writeLines(c('a<-commandArgs(trailingOnly=TRUE)','writeLines("viewer-called",a[1])',
             'quit(save="no",status=99L)'),sentinel)
calls<-list();marker<-NULL
plotenv$system<-function(command,...) {
  status<-suppressWarnings(system2(rscript,c('--vanilla',shQuote(sentinel),shQuote(marker))))
  calls[[length(calls)+1L]]<<-list(command=command,status=as.integer(status))
  stopifnot(status==99L,file.exists(marker))
  invisible(status)
}
pdf_info<-function(p) {
  stopifnot(file.exists(p),file.info(p)$size>100L)
  bytes<-readBin(p,'raw',n=file.info(p)$size)
  stopifnot(rawToChar(bytes[1:5])=='%PDF-',
            grepl('%%EOF',rawToChar(tail(bytes,50)),fixed=TRUE))
  list(path=p,size_bytes=unname(file.info(p)$size),sha256=digest::digest(file=p,algo='sha256'))
}
invoke<-function(name,folder,include_null_range=TRUE) {
  dir.create(folder)
  if (name=='plot_simDEC_DECJ_inferences') {
    simdir<-file.path(folder,'ps0001_sim001');dir.create(simdir)
    res<-dec_display;save(res,file=file.path(simdir,'DEC_inf.Rdata'))
    res<-decj_display;save(res,file=file.path(simdir,'DECJ_inf.Rdata'))
    SSEsim_results_processed<-list(engineering_metadata_only=TRUE)
    save(SSEsim_results_processed,file=file.path(simdir,'SSEsim_results_processed.Rdata'))
    dej_params<-data.frame(d=.05,e=.02,j=.03,brate=.3,drate=.3,b_exp=1,d_exp=-1)
    params<-file.path(folder,'dej_params.Rdata');save(dej_params,file=params)
    value<-withVisible(plotenv[[name]](simdir=simdir,sim_params_Rdata_table=params,
                                      makePDFs=TRUE,openPDFs=TRUE))
    pdf<-file.path(simdir,'ps0001_sim001_DEC_vs_DECj_SSE.pdf')
    stopifnot(all(c('SSEsim_results_processed','res1','res2','tests') %in% names(value$value)))
    saveRDS(value$value,file.path(folder,'returned_tables.rds'))
  } else if (name=='plot_BSM_fns') {
    pdf<-file.path(folder,'saved_bsm.pdf')
    value<-withVisible(plotenv[[name]](res1=if(include_null_range)decj else decj_no_empty,
                                      BSM_fns_res1=file.path(inputdir,if(include_null_range)'one_saved_map.Rdata' else 'one_saved_map_no_empty_display.Rdata'),
                                      res1name='saved DEC+J map; engineering display',pdffn=pdf,
                                      include_null_range=include_null_range))
    stopifnot(identical(value$value,pdf))
  } else {
    pi<-plot_inputs;pdf<-file.path(folder,paste0(name,'.pdf'));pi$pdffn<-pdf
    value<-withVisible(plotenv[[name]](pi,xwidth_manual=4))
  }
  list(pdf=pdf_info(pdf),returned_class=class(value$value),visible=value$visible,
       return_scalar=if(is.numeric(value$value)&&length(value$value)==1L)value$value else NULL)
}
results<-list()
cases<-c(targets,'plot_BSM_fns_no_empty_display')
for(case_name in cases) {
  name<-if(case_name=='plot_BSM_fns_no_empty_display')'plot_BSM_fns' else case_name
  include_null_range<-case_name!='plot_BSM_fns_no_empty_display'
  cat('\nSTART ',name,'\n',sep='');flush.console()
  saved_cwd<-getwd();warnings<-character()
  result<-tryCatch(withCallingHandlers({
    stopifnot(identical(body(plotenv[[name]]),original_defs[[name]]$body),
              identical(formals(plotenv[[name]]),original_defs[[name]]$formals))
    case<-file.path(out,case_name);dir.create(case)
    marker<-file.path(case,'positive_viewer_called.txt');calls<-list()
    plotenv$interactive<-function() TRUE
    positive<-invoke(name,file.path(case,'positive'),include_null_range)
    stopifnot(length(calls)==1L,calls[[1]]$status==99L,file.exists(marker))
    positive$viewer_calls<-calls
    positive$pdf$current_available<-TRUE
    stopifnot(!file.exists(file.path(case,'noninteractive')))
    plotenv$interactive<-base::interactive
    marker<-file.path(case,'negative_viewer_called.txt');calls<-list()
    negative<-invoke(name,file.path(case,'noninteractive'),include_null_range)
    stopifnot(length(calls)==0L,!file.exists(marker),identical(dev.cur(),c('null device'=1L)))
    if(name %in% targets[2:4])stopifnot(identical(negative$return_scalar,0L),!negative$visible)
    list(status='passed',positive_control=positive,noninteractive=negative,
         source_body_and_formals_unchanged=TRUE,viewer_calls_noninteractive=0L)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w))}),error=function(e)
    list(status='failed_retained',error=conditionMessage(e)))
  result$warnings<-warnings;results[[case_name]]<-result
  while(dev.cur()>1L)dev.off()
  setwd(saved_cwd)
  cat('END ',name,' ',result$status,'\n',sep='');flush.console()
}
manifest<-list(source_branch=branch,source_files=setNames(lapply(source_files,
  function(p)digest::digest(file=p,algo='sha256')),source_files),
  r_version=R.version.string,installed_package=find.package('RustBEARS'),
  input_source=original,scientific_fits_or_simulations_executed=FALSE,
  fixture_scope='Saved three-tip fitted objects and one saved BSM; deterministic six-row engineering summaries. No general plotting/numerical correctness claim.',
  implicit_globals=c('uniq_params_num','siminf_stats_subset','ancstate_accuracies_DEC_subset','ancstate_accuracies_DECJ_subset'),
  legacy_ape_bindings=ape_bindings,
  sse_display_adapter='Only copied optim container: use_optimx=FALSE and value=stored total_loglikelihood. No refitting, probability modification or score recalculation.',
  color_count_checks=color_checks,
  no_empty_display_adapter='Remove only an asserted zero-probability empty column and subtract one from all saved nonempty state numbers; retain event times, labels and nonempty probabilities. Engineering display, not a new fitted model or simulated history.',
  cases=results,passed=sum(vapply(results,function(x)x$status=='passed',logical(1))))
jsonlite::write_json(manifest,file.path(out,'results.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
saveRDS(manifest,file.path(out,'results.rds'))
stopifnot(manifest$passed==6L)
