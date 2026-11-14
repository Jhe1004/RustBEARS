# GPL >= 3. Conditional mapping precision and explicitly selected tree targets.
.bgb_u_map_matrix <- function(bsm) {
    n<-bsm$num_maps
    if(length(n)!=1L||!is.finite(n)||n<2||n!=floor(n))stop('at least two maps are required')
    types<-strsplit('deajysvt','')[[1]]
    events<-matrix(0,n,length(types),dimnames=list(NULL,types))
    routes<-matrix(0,n,length(bsm$areas)^2)
    names<-as.vector(t(outer(bsm$areas,bsm$areas,paste,sep=' -> ')))
    colnames(routes)<-names
    for(table in list(bsm$clado_events_df,bsm$ana_events_df)) {
        if(!nrow(table))next
        if(any(!is.finite(table$map))||any(table$map!=floor(table$map)))stop('invalid map indices')
        ids<-as.integer(table$map)
        if(any(ids<1L|ids>n))stop('map indices must cover a one-based complete batch')
        for(k in seq_len(nrow(table))) {
            type<-match(table$event_type[k],types)
            if(is.na(type))stop('unknown event type')
            events[ids[k],type]<-events[ids[k],type]+1
            if(table$event_type[k] %in% c('d','a','j')&&!is.na(table$source_area[k])&&!is.na(table$destination_area[k])) {
                route<-match(paste(table$source_area[k],table$destination_area[k],sep=' -> '),names)
                if(!is.na(route))routes[ids[k],route]<-routes[ids[k],route]+1
            }
        }
    }
    list(events=events,routes=routes)
}
.bgb_u_ips <- function(x) {
    n<-length(x)
    if(n<20L)stop('MCMC tree diagnostics require at least 20 ordered trees')
    z<-x-mean(x);g0<-sum(z*z)/n
    if(g0==0)return(list(tau=1,ess=NA_real_))
    gamma<-vapply(0:(n-1L),function(k)sum(z[seq_len(n-k)]*z[k+seq_len(n-k)])/n,numeric(1))
    pairs<-gamma[seq.int(1L,n-1L,by=2L)]+gamma[seq.int(2L,n,by=2L)]
    first_bad<-which(pairs<=0)[1]
    if(!is.na(first_bad))pairs<-if(first_bad==1L)numeric() else pairs[seq_len(first_bad-1L)]
    pairs<-cummin(pairs)
    tau<-max(1,(-g0+2*sum(pairs))/g0)
    list(tau=tau,ess=n/tau)
}
.bgb_u_weighted_quantile <- function(x,w,p) {
    ord<-order(x);x<-x[ord];w<-w[ord]/sum(w)
    vapply(p,function(q)x[which(cumsum(w)>=q)[1]],numeric(1))
}
.bgb_u_rank_status <- function(lower,upper,margin) {
    if(!all(is.finite(c(lower,upper,margin)))||lower>upper||margin<=0)stop('invalid paired interval')
    if(lower>=-margin&&upper<=margin)'near_tie' else
        if(lower>margin||upper< -margin)'resolved' else 'unresolved'
}
.bgb_u_scalar_summary <- function(x,strata,tree_target,tree_dependence,alpha) {
    strata<-strata[strata$stratum_weight>0,,drop=FALSE]
    groups<-lapply(strata$map_ids,function(ids)x[ids])
    ns<-lengths(groups);means<-vapply(groups,mean,numeric(1))
    variances<-vapply(groups,function(y)if(length(y)>1L)stats::var(y) else NA_real_,numeric(1))
    w<-strata$stratum_weight/sum(strata$stratum_weight)
    estimate<-sum(w*means);within<-sum(w*w*variances/ns)
    tree_ids<-unique(strata$tree_id)
    tree_w<-vapply(tree_ids,function(id)sum(w[strata$tree_id==id]),numeric(1))
    tree_means<-vapply(tree_ids,function(id){keep<-strata$tree_id==id;sum(w[keep]*means[keep])/sum(w[keep])},numeric(1))
    tree_mapping_var<-vapply(tree_ids,function(id){keep<-strata$tree_id==id;ww<-w[keep]/sum(w[keep]);sum(ww^2*variances[keep]/ns[keep])},numeric(1))
    denominator<-1-sum(w*w)
    between<-if(denominator>1e-15)max(0,sum(w*(means-estimate)^2-w*(1-w)*variances/ns)/denominator) else 0
    if(length(tree_ids)>1L)between<-max(0,sum(tree_w*(tree_means-estimate)^2-
        tree_w*(1-tree_w)*tree_mapping_var)/(1-sum(tree_w^2)))
    ess<-NA_real_;se<-sqrt(within)
    if(tree_target=='posterior_population') {
        k<-length(tree_ids)
        if(k<2L)stop('tree-population precision requires at least two trees')
        if(max(abs(tree_w-1/k))>1e-12)stop('tree-population precision currently requires equal tree weights')
        if(is.null(tree_dependence))stop('posterior_population requires explicit tree_dependence')
        ess<-if(tree_dependence=='mcmc_chain').bgb_u_ips(tree_means)$ess else k
        se<-if(stats::var(tree_means)==0)0 else sqrt(stats::var(tree_means)/ess)
    }
    map_weights<-numeric(length(x))
    for(i in seq_along(groups))map_weights[strata$map_ids[[i]]]<-w[i]/ns[i]
    q<-.bgb_u_weighted_quantile(x,map_weights,c(.025,.975))
    z<-stats::qnorm(1-alpha/2)
    stratum_q<-.bgb_u_weighted_quantile(means,w,c(.025,.975))
    list(mean=estimate,rpi_025=q[1],rpi_975=q[2],se_target=se,mapping_mcse=sqrt(within),
         mean_mc_ci_lower=estimate-z*se,mean_mc_ci_upper=estimate+z*se,
         var_within=sum(w*variances),var_between_trunc=between,stratum_mean_q025=stratum_q[1],
         stratum_mean_q975=stratum_q[2],effective_n_trees=ess)
}
bsm_diagnostics <- function(bsm_output,res=NULL,target_rmcse=.05,tie_margin=.25,
        tree_target=c('fixed_sample','posterior_population'),tree_dependence=NULL,
        top_k_routes=5L,min_mean_count=1,alpha=.05) {
    if(!inherits(bsm_output,'BSMResult'))stop('diagnostics require native BSMResult output')
    tree_target<-match.arg(tree_target)
    if(!is.null(tree_dependence))tree_dependence<-match.arg(tree_dependence,c('independent','mcmc_chain'))
    if(tree_target=='fixed_sample'&&!is.null(tree_dependence))stop('tree_dependence applies only to posterior_population')
    if(tree_target=='posterior_population'&&!is.null(bsm_output$n_trees_missing)&&bsm_output$n_trees_missing>0)
        stop('posterior_population precision is unavailable after incomplete-tree selection')
    if(!is.finite(target_rmcse)||target_rmcse<=0||!is.finite(tie_margin)||tie_margin<=0||
       !is.finite(alpha)||alpha<=0||alpha>=1||!is.finite(min_mean_count)||min_mean_count<0||
       !is.finite(top_k_routes)||top_k_routes<1||top_k_routes!=floor(top_k_routes))stop('invalid precision or interval settings')
    matrices<-.bgb_u_map_matrix(bsm_output)
    strata<-bsm_output$strata
    if(is.null(strata)) {
        strata<-data.frame(tree_id='1',model_id='1',stratum_weight=1)
        strata$map_ids<-list(seq_len(bsm_output$num_maps))
    }
    if(!is.null(bsm_output$tree_order)) {
        ids_order<-as.character(bsm_output$tree_order)
        if(anyDuplicated(ids_order)||!setequal(ids_order,as.character(strata$tree_id)))stop('tree_order must enumerate every tree once')
        strata<-strata[order(match(as.character(strata$tree_id),ids_order)),,drop=FALSE]
    }
    ids<-unlist(strata$map_ids,use.names=FALSE)
    if(any(!is.finite(ids))||any(ids!=floor(ids))||!identical(sort(as.integer(ids)),seq_len(bsm_output$num_maps))||
       any(!is.finite(strata$stratum_weight))||any(strata$stratum_weight<=0)||abs(sum(strata$stratum_weight)-1)>1e-12||
       any(lengths(strata$map_ids)<2L))
        stop('positive normalized strata must partition all maps, with at least two maps per stratum')
    summarize<-function(matrix) {
        rows<-lapply(seq_len(ncol(matrix)),function(i) {
            x<-.bgb_u_scalar_summary(matrix[,i],strata,tree_target,tree_dependence,alpha)
            x$statistic<-colnames(matrix)[i]
            x$assessment<-if(x$mean>=min_mean_count&&x$mean>0&&is.finite(x$se_target))'assessed' else 'not_assessed'
            if(x$assessment=='not_assessed')x$mean_mc_ci_lower<-x$mean_mc_ci_upper<-NA_real_
            x$rmcse<-if(x$mean>0)x$se_target/x$mean else NA_real_
            as.data.frame(x,stringsAsFactors=FALSE)
        })
        if(!length(rows))return(data.frame())
        do.call(rbind,rows)
    }
    diagonal<-vapply(bsm_output$areas,function(a)paste(a,a,sep=' -> '),character(1))
    matrices$routes<-matrices$routes[,!colnames(matrices$routes) %in% diagonal,drop=FALSE]
    event<-summarize(matrices$events);route<-summarize(matrices$routes)
    rank<-if(nrow(route))order(-route$mean,route$statistic) else integer()
    rank<-head(rank[route$assessment[rank]=='assessed'],top_k_routes)
    pairs<-list()
    if(length(rank)>1L)for(i in seq_len(length(rank)-1L)) {
        a<-rank[i];b<-rank[i+1L]
        x<-.bgb_u_scalar_summary(matrices$routes[,a]-matrices$routes[,b],strata,tree_target,tree_dependence,alpha)
        assessed<-route$mean[a]>=min_mean_count&&route$mean[b]>=min_mean_count&&is.finite(x$se_target)
        lo<-x$mean_mc_ci_lower;hi<-x$mean_mc_ci_upper
        status<-if(!assessed)'unresolved' else .bgb_u_rank_status(lo,hi,tie_margin)
        pairs[[i]]<-data.frame(route_a=route$statistic[a],route_b=route$statistic[b],mean_diff=x$mean,
            se_diff=x$se_target,ci_lower=lo,ci_upper=hi,tie_margin=tie_margin,zero_excluded=lo>0||hi<0,
            status=status,assessment=if(assessed)'assessed' else 'not_assessed',effective_n_trees=x$effective_n_trees)
    }
    active<-c(event$rmcse[event$assessment=='assessed'],route$rmcse[route$assessment=='assessed'])
    maximum<-if(length(active))max(active) else Inf
    reason<-bsm_output$stop_reason;if(is.null(reason))reason<-'fixed_count'
    all_events<-rbind(bsm_output$ana_events_df[,c('map','time_bp','event_type')],
                      bsm_output$clado_events_df[,c('map','time_bp','event_type')])
    ett<-data.frame()
    if(nrow(all_events)) {
        height<-max(max(all_events$time_bp),1e-12)
        if(height>0) {
            breaks<-seq(0,height,length.out=21L)
            bins<-pmin(20L,findInterval(all_events$time_bp,breaks,rightmost.closed=TRUE))
            counts<-matrix(0,bsm_output$num_maps,20L)
            for(i in seq_len(nrow(all_events)))counts[all_events$map[i],max(1L,bins[i])]<-counts[all_events$map[i],max(1L,bins[i])]+1
            colnames(counts)<-as.character(seq_len(20L));ett<-summarize(counts)
            ett$time_from<-breaks[-length(breaks)];ett$time_to<-breaks[-1L]
        }
    }
    structure(list(sampling_design=if(length(unique(strata$tree_id))>1L)'two_stage_trees' else
        if(nrow(strata)>1L)'model_stratified' else 'single_tree_iid',stop_reason=reason,
        tree_target=tree_target,tree_dependence=tree_dependence,num_maps_total=bsm_output$num_maps,
        precision_target_met=is.finite(maximum)&&maximum<=target_rmcse,target_rmcse=target_rmcse,
        tie_margin=tie_margin,max_observed_rmcse=maximum,event_summary_df=event,route_summary_df=route,
        route_rank_comparisons_df=if(length(pairs))do.call(rbind,pairs) else data.frame(),ett_df=ett,
        interval_scope='approximate per-comparison exploratory; no simultaneous or optional-stopping guarantee'),
        class='rustbears_bsm_diag')
}
plot.rustbears_bsm_diag <- function(x,...) {
    z<-x$event_summary_df
    locations<-graphics::barplot(stats::setNames(z$mean,z$statistic),ylab='Mean events per history',...)
    graphics::arrows(locations,z$mean_mc_ci_lower,locations,z$mean_mc_ci_upper,angle=90,code=3,length=.05)
    invisible(x)
}
