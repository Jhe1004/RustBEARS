# GPL >= 3. Conditional inference over an external, ordered sample of trees.
.bgb_u_tree_members <- function(tree) {
    n<-length(tree$tip.label)+tree$Nnode;members<-vector('list',n)
    members[seq_along(tree$tip.label)]<-as.list(tree$tip.label)
    visit<-function(node) {
        if(length(members[[node]]))return(members[[node]])
        members[[node]]<<-sort(unlist(lapply(tree$edge[tree$edge[,1]==node,2],visit),use.names=FALSE))
        members[[node]]
    }
    for(i in seq_len(n))visit(i)
    members
}
.bgb_u_tree_key <- function(tree) {
    members<-.bgb_u_tree_members(tree)
    keys<-vapply(members,function(x)digest::digest(sort(x),algo='sha256'),character(1))
    list(members=members,keys=keys,topology_hash=.bgb_u_topology_hash(tree,members))
}
.bgb_u_read_trees <- function(trees) {
    if(inherits(trees,'phylo'))return(list(trees))
    if(inherits(trees,'multiPhylo'))return(unclass(trees))
    if(is.character(trees)&&length(trees)==1L&&file.exists(trees)) {
        first<-readLines(trees,n=1L,warn=FALSE)
        x<-if(grepl('NEXUS',first,ignore.case=TRUE))ape::read.nexus(trees) else ape::read.tree(trees)
        return(if(inherits(x,'phylo'))list(x) else unclass(x))
    }
    if(is.character(trees))return(lapply(trees,function(x)ape::read.tree(text=x)))
    if(is.list(trees)&&all(vapply(trees,inherits,logical(1),'phylo')))return(trees)
    stop('trees must be phylo objects, Newick strings or a Newick/Nexus file')
}
.bgb_u_model_definition <- function(template,name) {
    if(methods::is(name,'BioGeoBEARS_model'))return(name)
    if(!is.character(name)||length(name)!=1L)stop('invalid model definition')
    base<-toupper(strsplit(name,'+',fixed=TRUE)[[1]][1]);mods<-tolower(strsplit(name,'+',fixed=TRUE)[[1]][-1L])
    if(!base %in% c('DEC','DIVALIKE','BAYAREALIKE')||any(!mods %in% c('j','x','w'))||anyDuplicated(mods))
        stop('UnsupportedModelContractError: unknown standard model')
    p<-template@params_table
    p['j','type']<-if('j' %in% mods)'free' else 'fixed'
    p['j',c('init','est')]<-if('j' %in% mods).01 else 0
    categories<-switch(base,DEC=3,DIVALIKE=2,BAYAREALIKE=1)
    p['j','min']<-1e-5;p['j','max']<-categories-1e-5
    p[c('d','e'),'type']<-'free';p[c('d','e'),'min']<-1e-12
    p['ysv','type']<-paste0(categories,'-j')
    p['ys','type']<-if(base=='DEC')'ysv*2/3' else if(base=='DIVALIKE')'ysv*1/2' else 'ysv'
    for(k in c('y','s','v')) {
        is_zero<-(k=='s'&&base!='DEC')||(k=='v'&&base=='BAYAREALIKE')
        p[k,'type']<-if(is_zero)'fixed' else if(base=='DEC')'ysv*1/3' else if(base=='DIVALIKE')'ysv*1/2' else 'ysv'
        p[k,c('init','est')]<-if(is_zero)0 else 1
    }
    p['mx01v','type']<-if(base=='DIVALIKE')'fixed' else 'mx01'
    p['mx01v',c('init','est')]<-if(base=='DIVALIKE').5 else p['mx01','est']
    p['mx01y','type']<-if(base=='BAYAREALIKE')'fixed' else 'mx01'
    p['mx01y',c('init','est')]<-if(base=='BAYAREALIKE').9999 else p['mx01','est']
    for(k in c('x','w')) {
        p[k,'type']<-if(k %in% mods)'free' else 'fixed'
        p[k,c('init','est')]<-if(k=='x')0 else 1
    }
    template@params_table<-p;calc_linked_params_BioGeoBEARS_model_object(template,update_init=TRUE)
}
.bgb_u_fit_tree <- function(template,tree,definition,directory,start=NULL,optimize=TRUE,num_threads=1L) {
    x<-template
    x$trfn<-file.path(directory,'tree.newick');ape::write.tree(tree,file=x$trfn,digits=17)
    x$tmpwd<-directory;x$wd<-directory;x$BioGeoBEARS_model_object<-definition
    if(length(start))for(k in intersect(names(start),rownames(definition@params_table)))
        if(definition@params_table[k,'type']=='free')
            x$BioGeoBEARS_model_object@params_table[k,c('init','est')]<-max(definition@params_table[k,'min'],min(definition@params_table[k,'max'],start[[k]]))
    x$num_cores_to_use<-as.integer(num_threads);x$printlevel<--1L;x$print_optim<-FALSE;x$calc_ancprobs<-TRUE;x$use_optimx<-FALSE
    x<-readfiles_BioGeoBEARS_run(x);result<-NULL
    invisible(utils::capture.output(result<-bears_optim_run(x,
        skip_optim=!optimize||!any(definition@params_table$type=='free'),skip_optim_option='return_all')))
    result
}
fit_over_trees <- function(BioGeoBEARS_run_object,trees,summary_tree=NULL,models=NULL,
        likelihood_conditioning=c('legacy_unconditioned','tip_survival','unconditional_complete'),
        root_prior_spec=NULL,match_mode=c('monophyletic','mrca'),
        init_strategy=c('warm_hybrid','full_cold_start','fixed_ref_params'),cold_start_audit_fraction=.05,
        model_averaging_within_tree=TRUE,max_tree_failure_rate=.02,
        on_excess_failure=c('abort','conservative_bounds'),bsm_maps_per_tree=0L,max_trees=NULL,
        num_threads=getOption('RustBEARS.threads',1L),seedval=12345L,checkpoint_file=NULL,
        criterion=c('AICc','AIC')) {
    started<-proc.time()[['elapsed']]
    likelihood_conditioning<-match.arg(likelihood_conditioning);match_mode<-match.arg(match_mode)
    init_strategy<-match.arg(init_strategy);on_excess_failure<-match.arg(on_excess_failure);criterion<-match.arg(criterion)
    if(length(num_threads)!=1L||!is.finite(num_threads)||num_threads<1||num_threads>256||num_threads!=floor(num_threads))stop('invalid num_threads')
    if(length(bsm_maps_per_tree)!=1L||!is.finite(bsm_maps_per_tree)||bsm_maps_per_tree<0||bsm_maps_per_tree!=floor(bsm_maps_per_tree))stop('invalid bsm_maps_per_tree')
    if(!is.finite(cold_start_audit_fraction)||cold_start_audit_fraction<0||cold_start_audit_fraction>1||
       !is.finite(max_tree_failure_rate)||max_tree_failure_rate<0||max_tree_failure_rate>1)stop('invalid audit or failure fraction')
    sample<-.bgb_u_read_trees(trees)
    if(!is.null(max_trees)) {
        if(length(max_trees)!=1L||!is.finite(max_trees)||max_trees<1||max_trees!=floor(max_trees))stop('invalid max_trees')
        sample<-head(sample,max_trees)
    }
    if(!length(sample))stop('empty tree sample')
    summary<-if(is.null(summary_tree))sample[[1]] else .bgb_u_read_trees(summary_tree)[[1]]
    taxa<-sort(summary$tip.label)
    if(anyDuplicated(taxa)||any(!vapply(sample,function(x)identical(sort(x$tip.label),taxa)&&!anyDuplicated(x$tip.label),logical(1))))
        stop('all trees must have the same unique taxa')
    if(any(!vapply(c(sample,list(summary)),function(t)ape::is.rooted(t)&&ape::is.binary(t)&&
        all(is.finite(t$edge.length))&&all(t$edge.length>=0),logical(1))))stop('only rooted binary trees with nonnegative branch lengths are supported')
    template<-BioGeoBEARS_run_object
    template$likelihood_conditioning<-likelihood_conditioning;template$root_prior<-root_prior_spec
    template$allow_null_tips<-likelihood_conditioning=='unconditional_complete'
    if(is.null(models))models<-list(template$BioGeoBEARS_model_object)
    if(is.character(models))models<-as.list(models)
    definitions<-lapply(models,function(x).bgb_u_model_definition(template$BioGeoBEARS_model_object,x))
    if(!length(definitions))stop('empty candidate model set')
    labels<-if(!is.null(names(models))&&all(nzchar(names(models))))names(models) else
        vapply(seq_along(models),function(i)if(is.character(models[[i]]))models[[i]] else paste0('model_',i),character(1))
    if(anyDuplicated(labels)||anyDuplicated(vapply(definitions,function(x)digest::digest(x@params_table,algo='sha256'),character(1))))
        stop('candidate model names and definitions must be unique')
    geography<-readLines(template$geogfn,warn=FALSE)
    keys<-lapply(sample,.bgb_u_tree_key);wanted<-.bgb_u_tree_key(summary)
    ntree<-length(sample);nn<-length(wanted$keys);nmodel<-length(models)
    mapping<-lapply(keys,function(x) {
        idx<-match(wanted$keys,x$keys);mono<-!is.na(idx)
        if(match_mode=='mrca')for(i in which(!mono)) {
            candidates<-which(vapply(x$members,function(m)all(wanted$members[[i]] %in% m),logical(1)))
            idx[i]<-candidates[which.min(lengths(x$members)[candidates])]
        }
        excess<-numeric(nn);keep<-!is.na(idx);excess[keep]<-lengths(x$members)[idx[keep]]-lengths(wanted$members)[keep]
        list(index=idx,monophyletic=mono,excess=excess)
    })
    distance<-vapply(keys,function(x)sum(!x$keys %in% wanted$keys),integer(1))
    audit<-if(cold_start_audit_fraction>0)head(order(-distance,seq_len(ntree)),max(1L,ceiling(ntree*cold_start_audit_fraction))) else integer()
    checkpoint<-.bgb_u_checkpoint(checkpoint_file,list(kind='multitree',software=as.character(utils::packageVersion('RustBEARS')),
        ffi_abi=2L,precision_mode=.bgb_precision(),trees=sample,summary_tree=summary,geography=geography,
        template=template,definitions=definitions,labels=labels,likelihood_conditioning=likelihood_conditioning,
        root_prior_spec=root_prior_spec,match_mode=match_mode,init_strategy=init_strategy,
        cold_start_audit_fraction=cold_start_audit_fraction,model_averaging_within_tree=model_averaging_within_tree,
        max_tree_failure_rate=max_tree_failure_rate,on_excess_failure=on_excess_failure,
        bsm_maps_per_tree=bsm_maps_per_tree,max_trees=max_trees,seedval=seedval,criterion=criterion))
    directory<-tempfile('rustbears-trees-');dir.create(directory);on.exit(unlink(directory,recursive=TRUE),add=TRUE)
    # Frozen geography remains valid even if the caller later changes its file.
    template$geogfn<-file.path(directory,'geography.data');writeLines(geography,template$geogfn)
    reference<-vector('list',nmodel);records<-checkpoint$records
    seed_tasks<-lapply(seq_len(ntree),function(t)lapply(seq_len(nmodel),function(m)
        .bgb_u_task_seed(seedval,'multitree_bsm',as.character(m),as.character(t),role='draw')))
    all_seeds<-unlist(lapply(seed_tasks,function(x)vapply(x,`[[`,character(1),'seed')))
    if(anyDuplicated(all_seeds))stop('task seed collision')
    for(t in seq_len(ntree)) {
        key<-as.character(t)
        if(!is.null(records[[key]])) {
            if(t==1L)reference<-records[[key]]$reference
            next
        }
        fitted<-vector('list',nmodel);rows<-vector('list',nmodel);errors<-character(nmodel);audit_delta<-numeric()
        for(m in seq_len(nmodel)) {
            initial<-if(init_strategy=='full_cold_start')NULL else reference[[m]]
            attempt<-if(init_strategy=='fixed_ref_params'&&t>1L&&!length(reference[[m]]))simpleError('reference-tree parameter fit is unavailable') else
                tryCatch(.bgb_u_fit_tree(template,sample[[t]],definitions[[m]],directory,initial,
                optimize=init_strategy!='fixed_ref_params'||t==1L,num_threads=num_threads),error=function(e)e)
            valid<-!inherits(attempt,'error')&&if(init_strategy=='fixed_ref_params'&&t>1L)
                is.finite(attempt$total_loglikelihood) else .bgb_u_converged(attempt)
            retried<-FALSE
            if(!valid&&init_strategy!='fixed_ref_params') {
                retried<-TRUE
                attempt<-tryCatch(.bgb_u_fit_tree(template,sample[[t]],definitions[[m]],directory,
                    optimize=TRUE,num_threads=num_threads),error=function(e)e)
                valid<-!inherits(attempt,'error')&&.bgb_u_converged(attempt)
            }
            if(valid) {
                if(t==1L)reference[[m]]<-.bgb_u_params(attempt)
                if(init_strategy=='warm_hybrid'&&t %in% audit&&t>1L) {
                    cold<-tryCatch(.bgb_u_fit_tree(template,sample[[t]],definitions[[m]],directory,
                        optimize=TRUE,num_threads=num_threads),error=function(e)e)
                    if(!inherits(cold,'error')&&.bgb_u_converged(cold)) {
                        delta<-cold$total_loglikelihood-attempt$total_loglikelihood;audit_delta<-c(audit_delta,delta)
                        if(delta>1e-3)warning('cold-start audit found a better optimum; consider init_strategy="full_cold_start"')
                        if(delta>0)attempt<-cold
                    } else warning('cold-start audit did not converge')
                }
                fitted[[m]]<-attempt
            } else errors[m]<-if(inherits(attempt,'error'))conditionMessage(attempt) else 'optimizer did not converge'
            rows[[m]]<-data.frame(tree_id=as.character(t),tree_topology_hash=keys[[t]]$topology_hash,model_id=labels[m],
                obj=if(valid)attempt$total_loglikelihood else NA_real_,converged=valid,retried=retried,error=errors[m],
                precision_mode=if(.bgb_precision()==1L)'strict_f64' else 'legacy_f32',likelihood_conditioning=likelihood_conditioning)
        }
        complete<-all(vapply(fitted,function(f)!is.null(f),logical(1)))
        record<-list(status=if(complete)'complete' else 'incomplete_model_set',rows=do.call(rbind,rows),reference=reference,
            audit_delta=audit_delta,fits=fitted,top=NULL,bottom=NULL,bsm=NULL,weights=NULL,params=NULL)
        if(complete) {
            names(fitted)<-labels
            average<-model_average(fitted,criterion=criterion)
            w<-average$weights
            if(!isTRUE(model_averaging_within_tree)) {w[]<-0;w[which.min(average$model_table$criterion_value)]<-1}
            top<-average$ancestral_probs_branch_top*0;bottom<-average$ancestral_probs_branch_bottom*0
            for(m in seq_len(nmodel))if(w[m]>0) {
                top<-top+w[m]*fitted[[m]]$ML_marginal_prob_each_state_at_branch_top_AT_node
                bottom<-bottom+w[m]*fitted[[m]]$ML_marginal_prob_each_state_at_branch_bottom_below_node
            }
            idx<-mapping[[t]]$index
            record$top<-top[idx,,drop=FALSE];record$bottom<-bottom[idx,,drop=FALSE]
            record$weights<-w;record$rows$criterion_value<-average$model_table$criterion_value;record$rows$weight<-w
            record$params<-lapply(fitted,.bgb_u_params)
            if(bsm_maps_per_tree>0) {
                allocation<-.bgb_u_allocate_maps(w,bsm_maps_per_tree);active<-which(allocation>0)
                parts<-lapply(active,function(m) {
                    sampler<-prepare_tree_sampler(fitted[[m]]);on.exit(close_prepared_sampler(sampler),add=TRUE)
                    .bgb_u_attach_rng(draw_prepared_bsm(sampler,allocation[m],seed=seed_tasks[[t]][[m]]$seed,num_threads=num_threads),seed_tasks[[t]][[m]])
                })
                record$bsm<-.bgb_u_combine_maps(parts,w[active],tree_ids=rep(t,length(active)),model_ids=active)
                record$bsm$tree_topology_hash<-keys[[t]]$topology_hash
            }
        } else {record$rows$criterion_value<-NA_real_;record$rows$weight<-NA_real_}
        record$rows$tree_status<-record$status
        records[[key]]<-record;checkpoint$records<-records;.bgb_u_save_checkpoint(checkpoint)
    }
    records<-records[as.character(seq_len(ntree))]
    ok<-vapply(records,function(x)x$status=='complete',logical(1));nok<-sum(ok)
    narea<-ncol(getranges_from_LagrangePHYLIP(template$geogfn)@df)
    nstate<-if(nok>0L)ncol(records[[which(ok)[1]]]$top) else
        sum(choose(narea,seq_len(min(narea,template$max_range_size))))+as.integer(isTRUE(template$include_null_range))
    summarize<-function(field) {
        value<-array(NA_real_,c(nn,nstate,ntree))
        for(t in which(ok))value[,,t]<-records[[t]][[field]]
        eligible<-sapply(mapping,function(x)if(match_mode=='monophyletic')x$monophyletic else rep(TRUE,nn))
        count<-rowSums(eligible[,ok,drop=FALSE]);denominator<-rowSums(eligible)
        total<-apply(value,c(1,2),sum,na.rm=TRUE)
        mean<-total/count;mean[count==0,]<-NA_real_
        low<-total/denominator;high<-(total+denominator-count)/denominator
        low[denominator==0,]<-NA_real_;high[denominator==0,]<-NA_real_
        sd<-apply(value,c(1,2),stats::sd,na.rm=TRUE)
        if(field=='bottom') {
            # A root has no incoming branch. Such MRCA mappings are undefined.
            undefined<-vapply(seq_len(nn),function(u)any(vapply(seq_len(ntree),function(t)
                !is.na(mapping[[t]]$index[u])&&length(keys[[t]]$members[[mapping[[t]]$index[u]]])==length(taxa),logical(1))),logical(1))
            mean[undefined,]<-NA_real_;low[undefined,]<-NA_real_;high[undefined,]<-NA_real_;sd[undefined,]<-NA_real_
        }
        list(mean=mean,low=low,high=high,sd=sd,count=count)
    }
    top<-summarize('top');bottom<-summarize('bottom')
    joint<-NULL
    if(bsm_maps_per_tree>0&&nok>0) {
        parts<-lapply(records[ok],`[[`,'bsm');joint<-.bgb_u_combine_maps(parts,rep(1/nok,nok),tree_ids=which(ok))
        joint$summary_scope<-'complete-case trees only; event counts are unbounded'
        joint$n_trees_missing<-ntree-nok
        joint$tree_order<-as.character(which(ok))
    }
    winner<-vapply(records[ok],function(x)labels[which.max(x$weights)],character(1))
    sensitivity<-list()
    if(nok>0)for(m in seq_len(nmodel))for(p in rownames(definitions[[m]]@params_table)[definitions[[m]]@params_table$type=='free']) {
        vals<-vapply(records[ok],function(x)x$params[[m]][p],numeric(1));q<-stats::quantile(vals,c(.025,.5,.975),names=FALSE)
        sensitivity[[length(sensitivity)+1L]]<-data.frame(model_id=labels[m],param=p,q025=q[1],median=q[2],q975=q[3],
            interpretation='cross-tree conditional point-estimate sensitivity; not a confidence interval')
    }
    audit_values<-unlist(lapply(records,`[[`,'audit_delta'))
    result<-structure(list(summary_tree=summary,n_trees_planned=ntree,n_trees_completed=nok,n_trees_incomplete=ntree-nok,
        models=labels,likelihood_conditioning=likelihood_conditioning,precision_mode=records[[1]]$rows$precision_mode[1],
        match_mode=match_mode,cold_start_audit_max_delta=if(length(audit_values))max(audit_values) else NULL,
        clade_support_all_trees=rowMeans(sapply(mapping,`[[`,'monophyletic')),clade_complete_trees_count=top$count,
        mean_mrca_excess_tips=rowMeans(sapply(mapping,`[[`,'excess')),complete_case_probs_top=top$mean,
        complete_case_probs_bottom=bottom$mean,prob_failure_bounds_top=list(lower=top$low,upper=top$high),
        prob_failure_bounds_bottom=list(lower=bottom$low,upper=bottom$high),cross_tree_prob_sd_top=top$sd,
        per_tree_df=do.call(rbind,lapply(records,`[[`,'rows')),model_winner_freq=if(nok) table(factor(winner,levels=labels))/nok else stats::setNames(rep(NA_real_,length(labels)),labels),
        param_tree_sensitivity_df=if(length(sensitivity))do.call(rbind,sensitivity) else data.frame(),joint_bsm=joint,
        elapsed_sec=proc.time()[['elapsed']]-started,task_fingerprint=checkpoint$fingerprint,
        interpretation='conditional ML/Akaike summaries over external trees; complete-case means with deterministic missing-fit bounds'),class='rustbears_multitree')
    if((ntree-nok)/ntree>max_tree_failure_rate&&on_excess_failure=='abort')
        stop(structure(list(message='TreeFailureError: incomplete tree fraction exceeds max_tree_failure_rate',call=NULL,result=result),
            class=c('rustbears_tree_failure','error','condition')))
    result
}
plot.rustbears_multitree <- function(x,min_clade_support=.5,...) {
    if(!is.finite(min_clade_support)||min_clade_support<0||min_clade_support>1)stop('invalid support threshold')
    z<-x$complete_case_probs_top;z[x$clade_support_all_trees<min_clade_support,]<-NA_real_
    graphics::image(seq_len(nrow(z)),seq_len(ncol(z)),z,xlab='Summary-tree node',ylab='State',zlim=c(0,1),...)
    invisible(x)
}
