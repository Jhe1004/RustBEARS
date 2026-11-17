# GPL >= 3. Plug-in conditional model averaging, not a joint Bayesian posterior.
.bgb_u_clades <- function(res) {
    phy<-ape::read.tree(text=paste(res$input_snapshot$tree_text,collapse='\n'))
    nt<-length(phy$tip.label);nn<-nt+phy$Nnode
    members<-vector('list',nn)
    members[seq_len(nt)]<-lapply(phy$tip.label,as.character)
    visit<-function(node) {
        if(length(members[[node]]))return(members[[node]])
        children<-phy$edge[phy$edge[,1]==node,2]
        members[[node]]<<-sort(unlist(lapply(children,visit),use.names=FALSE))
        members[[node]]
    }
    for(i in seq_len(nn))visit(i)
    vapply(members,function(x)digest::digest(x,algo='sha256'),character(1))
}
.bgb_u_model_label <- function(res) {
    p<-res$outputs@params_table
    name<-if(p['s','est']>0)'DEC' else if(p['v','est']>0)'DIVALIKE' else 'BAYAREALIKE'
    if(p['j','type']=='free'||p['j','est']>0)name<-paste0(name,'+J')
    if(p['x','type']=='free'||p['x','est']!=0)name<-paste0(name,'+x')
    if(p['w','type']=='free'||p['w','est']!=1)name<-paste0(name,'+w')
    if('n' %in% rownames(p)&&(p['n','type']=='free'||p['n','est']!=0))name<-paste0(name,'+n')
    if('u' %in% rownames(p)&&(p['u','type']=='free'||p['u','est']!=0))name<-paste0(name,'+u')
    name
}
.bgb_u_allocate_maps <- function(weights,total) {
    if(length(total)!=1L||!is.finite(total)||total<1||total!=floor(total))stop('total maps must be a positive integer')
    active<-which(weights>0)
    if(total<2L*length(active))stop('at least two maps per positive-weight model are required')
    n<-integer(length(weights));n[active]<-2L
    raw<-(total-sum(n))*weights/sum(weights);n<-n+floor(raw)
    left<-total-sum(n)
    if(left>0)n[head(order(-(raw-floor(raw)),seq_along(raw)),left)]<-n[head(order(-(raw-floor(raw)),seq_along(raw)),left)]+1L
    n
}
.bgb_u_combine_maps <- function(parts,weights,tree_ids=NULL,model_ids=NULL) {
    if(!length(parts))return(NULL)
    clado<-ana<-strata<-metadata<-vector('list',length(parts));offset<-0L
    for(i in seq_along(parts)) {
        p<-parts[[i]];clado[[i]]<-p$clado_events_df;ana[[i]]<-p$ana_events_df
        clado[[i]]$map<-clado[[i]]$map+offset;ana[[i]]$map<-ana[[i]]$map+offset
        metadata[[i]]<-p$map_metadata_df
        if(is.null(metadata[[i]]$map))stop('BSM metadata is missing global map indices')
        metadata[[i]]$map<-as.double(metadata[[i]]$map+offset)
        metadata[[i]]$stratum_weight<-metadata[[i]]$stratum_weight*weights[i]
        if(!is.null(tree_ids))metadata[[i]]$tree_id<-tree_ids[i]
        if(!is.null(model_ids))metadata[[i]]$model_id<-model_ids[i]
        strata[i]<-list(p$strata)
        if(is.null(strata[[i]])) {
            strata[[i]]<-metadata[[i]][1L,setdiff(names(metadata[[i]]),c('map','map_id')),drop=FALSE]
            strata[[i]]$map_ids<-list(metadata[[i]]$map)
        } else {
            strata[[i]]$stratum_weight<-strata[[i]]$stratum_weight*weights[i]
            strata[[i]]$map_ids<-lapply(strata[[i]]$map_ids,function(ids)as.double(ids+offset))
            if(!is.null(tree_ids))strata[[i]]$tree_id<-tree_ids[i]
            if(!is.null(model_ids))strata[[i]]$model_id<-model_ids[i]
        }
        for(key in names(metadata[[i]])) {
            clado[[i]][[key]]<-metadata[[i]][[key]][match(clado[[i]]$map,metadata[[i]]$map)]
            ana[[i]][[key]]<-metadata[[i]][[key]][match(ana[[i]]$map,metadata[[i]]$map)]
        }
        offset<-offset+p$num_maps
    }
    structure(list(clado_events_df=do.call(rbind,clado),ana_events_df=do.call(rbind,ana),
        num_maps=offset,areas=parts[[1]]$areas,strata=do.call(rbind,strata),map_metadata_df=do.call(rbind,metadata),
        tree_topology_hash=parts[[1]]$tree_topology_hash,
        likelihood_conditioning=parts[[1]]$likelihood_conditioning,precision_mode=parts[[1]]$precision_mode,
        stop_reason='fixed_count'),class='BSMResult')
}
.bgb_u_attach_rng <- function(result,task) {
    for(field in c('map_metadata_df','clado_events_df','ana_events_df')) {
        result[[field]]$rng_seed_hex<-rep(sub('^0x','',task$seed),nrow(result[[field]]))
        result[[field]]$rng_task_sha256<-rep(task$fingerprint,nrow(result[[field]]))
    }
    result
}
model_average <- function(res_list,criterion=c('AICc','AIC'),num_maps=0L,
        num_threads=getOption('RustBEARS.threads',1L),seedval=12345L) {
    criterion<-match.arg(criterion);.bgb_u_compatible(res_list)
    if(any(!vapply(res_list,.bgb_u_converged,logical(1))))stop('model averaging requires converged fits')
    start<-proc.time()[['elapsed']]
    first<-res_list[[1]];clades<-.bgb_u_clades(first)
    nt<-length(ape::read.tree(text=paste(first$input_snapshot$tree_text,collapse='\n'))$tip.label)
    labels<-vapply(res_list,.bgb_u_model_label,character(1))
    if(anyDuplicated(labels)) {
        if(is.null(names(res_list))||any(!nzchar(names(res_list)))||anyDuplicated(names(res_list)))
            stop('models with the same family label require distinct explicit names')
        labels<-names(res_list)
    }
    scores<-vapply(res_list,function(x) {
        p<-x$outputs@params_table;k<-sum(p$type=='free'&p$min<p$max);aic<-2*k-2*x$total_loglikelihood
        if(criterion=='AIC')aic else if(nt<=k+1)Inf else aic+2*k*(k+1)/(nt-k-1)
    },numeric(1))
    if(any(!is.finite(scores)))stop('all candidate scores must be finite; undefined AICc requires an explicit choice of AIC')
    relative<-exp(-.5*(scores-min(scores)));weights<-relative/sum(relative)
    top<-bottom<-first$ML_marginal_prob_each_state_at_branch_top_AT_node*0
    bottom<-first$ML_marginal_prob_each_state_at_branch_bottom_below_node*0
    for(i in seq_along(res_list)) {
        f<-res_list[[i]];idx<-match(clades,.bgb_u_clades(f))
        if(anyNA(idx))stop('models must share the same rooted clades')
        if(weights[i]>0) {
            top<-top+weights[i]*f$ML_marginal_prob_each_state_at_branch_top_AT_node[idx,,drop=FALSE]
            bottom<-bottom+weights[i]*f$ML_marginal_prob_each_state_at_branch_bottom_below_node[idx,,drop=FALSE]
        }
    }
    aligned<-lapply(res_list,function(f)f$ML_marginal_prob_each_state_at_branch_top_AT_node[match(clades,.bgb_u_clades(f)),,drop=FALSE])
    if(any(!is.finite(top))||any(abs(rowSums(top)-1)>1e-10))stop('invalid ancestral probability arrays')
    model_var<-Reduce(`+`,Map(function(a,w)w*(a-top)^2,aligned,weights))
    minimum<-Reduce(pmin,aligned);maximum<-Reduce(pmax,aligned)
    sensitivity<-expand.grid(node=seq_len(nrow(top)),state=seq_len(ncol(top))-1L)
    sensitivity$probability<-as.vector(top);sensitivity$model_sd<-as.vector(sqrt(model_var))
    sensitivity$min_model_probability<-as.vector(minimum);sensitivity$max_model_probability<-as.vector(maximum)
    bsm<-NULL
    if(num_maps>0) {
        allocation<-.bgb_u_allocate_maps(weights,num_maps);active<-which(allocation>0)
        tasks<-lapply(active,function(i).bgb_u_task_seed(seedval,'model_average_bsm',as.character(i),'1',role='draw'))
        if(anyDuplicated(vapply(tasks,`[[`,character(1),'seed')))stop('duplicate task seeds')
        parts<-lapply(seq_along(active),function(k) {
            i<-active[k];s<-prepare_tree_sampler(res_list[[i]])
            on.exit(close_prepared_sampler(s),add=TRUE)
            .bgb_u_attach_rng(draw_prepared_bsm(s,allocation[i],seed=tasks[[k]]$seed,num_threads=num_threads),tasks[[k]])
        })
        bsm<-.bgb_u_combine_maps(parts,weights[active],model_ids=active)
    }
    structure(list(criterion=criterion,likelihood_conditioning=first$likelihood_conditioning,
        precision_mode=first$precision_mode,weights=stats::setNames(weights,labels),
        model_table=data.frame(model=labels,n_params=vapply(res_list,function(x){p<-x$outputs@params_table;sum(p$type=='free'&p$min<p$max)},integer(1)),
            lnL=vapply(res_list,`[[`,numeric(1),'total_loglikelihood'),criterion_value=scores,weight=weights),
        ancestral_probs_branch_top=top,ancestral_probs_branch_bottom=bottom,node_sensitivity_df=sensitivity,
        joint_bsm=bsm,interpretation='plug-in conditional likelihood-model mixture; not a joint Bayesian posterior',
        elapsed_sec=proc.time()[['elapsed']]-start),class='rustbears_model_average')
}
plot.rustbears_model_average <- function(x,...) {
    graphics::barplot(x$weights,ylab=paste(x$criterion,'weight'),...);invisible(x)
}
.bgb_u_legacy_runBSM <- runBSM
runBSM <- function(res,stochastic_mapping_inputs_list,maxnum_maps_to_try=1,nummaps_goal=maxnum_maps_to_try,
        maxtries_per_branch=40000,save_after_every_try=TRUE,savedir=getwd(),seedval=as.numeric(Sys.time()),
        wait_before_save=.1,master_nodenum_toPrint=0,...,bsm_mode='legacy',num_threads=NULL,
        bsm_engine=c('legacy_r_bridge','native_whole_tree'),target_rmcse=NULL,min_maps=200L,
        batch_size=200L,max_maps=10000L,min_mean_count=1) {
    bsm_engine<-match.arg(bsm_engine)
    if(bsm_engine=='legacy_r_bridge') {
        if(!is.null(target_rmcse))stop('adaptive precision stopping requires bsm_engine="native_whole_tree"')
        return(.bgb_u_legacy_runBSM(res,stochastic_mapping_inputs_list,maxnum_maps_to_try,nummaps_goal,
            maxtries_per_branch,save_after_every_try,savedir,seedval,wait_before_save,master_nodenum_toPrint,
            ...,bsm_mode=bsm_mode,num_threads=num_threads))
    }
    if(is.null(num_threads))num_threads<-getOption('RustBEARS.threads',1L)
    if(missing(seedval))seedval<-floor(seedval)
    sampler<-prepare_tree_sampler(res);on.exit(close_prepared_sampler(sampler),add=TRUE)
    if(is.null(target_rmcse))return(draw_prepared_bsm(sampler,nummaps_goal,seed=seedval,num_threads=num_threads,
        rejection=identical(bsm_mode,'rejection')))
    if(min_maps<2||batch_size<1||max_maps<min_maps||any(c(min_maps,batch_size,max_maps)!=floor(c(min_maps,batch_size,max_maps))))
        stop('invalid adaptive sampling limits')
    out<-NULL;done<-0L
    while(done<max_maps) {
        n<-min(batch_size,max_maps-done)
        piece<-draw_prepared_bsm(sampler,n,start=done+1L,seed=seedval,num_threads=num_threads,
            rejection=identical(bsm_mode,'rejection'))
        if(is.null(out))out<-piece else {
            out$clado_events_df<-rbind(out$clado_events_df,piece$clado_events_df)
            out$ana_events_df<-rbind(out$ana_events_df,piece$ana_events_df)
            out$map_metadata_df<-rbind(out$map_metadata_df,piece$map_metadata_df)
            out$num_maps<-out$num_maps+n
        }
        done<-done+n
        if(done>=min_maps) {
            diagnostic<-bsm_diagnostics(out,target_rmcse=target_rmcse,min_mean_count=min_mean_count)
            if(diagnostic$precision_target_met)break
        }
    }
    out$stop_reason<-if(diagnostic$precision_target_met)'target_precision_met' else 'max_maps_reached'
    out$diagnostics<-bsm_diagnostics(out,target_rmcse=target_rmcse,min_mean_count=min_mean_count)
    out
}
