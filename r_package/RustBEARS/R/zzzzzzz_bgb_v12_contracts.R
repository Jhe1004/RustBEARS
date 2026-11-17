# GPL >= 3. Explicit probability contracts; historical function formals survive.
.bgb_v12_runtime <- new.env(parent=emptyenv())
.bgb_v12_runtime$contract <- 0L
.bgb_v12_runtime$root_prior <- NULL
.bgb_v12_runtime$last_root <- NULL
.bgb_v12_runtime$last_log_survival <- NULL
.bgb_v12_contract_names <- c('legacy_unconditioned','tip_survival','unconditional_complete')
.bgb_v12_standard_clado <- function(w) {
    if(length(w)!=4L||any(!is.finite(w))||any(w<0))return(FALSE)
    eq<-function(a,b)abs(a-b)<1e-12
    (eq(w[1],w[2])&&eq(w[1],w[4])) ||
    (eq(w[1],0)&&eq(w[2],w[4])) || (eq(w[1],0)&&eq(w[2],0))
}

define_root_prior <- function(mode=c('uniform_nonnull','custom_weights','legacy_flat_weights'),custom_weights=NULL) {
    mode<-match.arg(mode)
    if(mode=='custom_weights') {
        if(!is.numeric(custom_weights)||!length(custom_weights)||any(!is.finite(custom_weights))||
           any(custom_weights<0)||!is.finite(sum(custom_weights))||sum(custom_weights)<=0)
            stop('root weights must be finite, nonnegative and have positive finite total')
        if(!is.null(names(custom_weights))&&anyDuplicated(names(custom_weights)))stop('duplicate root state labels')
    } else if(!is.null(custom_weights))stop('custom_weights requires mode="custom_weights"')
    structure(list(mode=mode,custom_weights=custom_weights),class='RootPriorSpec')
}

.bgb_v12_state_label <- function(state) {
    if(!length(state)||anyNA(state)||any(state<0))return('_')
    paste(sort(as.integer(state)),collapse=',')
}

.bgb_v12_resolve_root <- function(spec,states,contract,allowed=NULL) {
    if(is.null(spec))spec<-define_root_prior(if(contract==0L)'legacy_flat_weights' else 'uniform_nonnull')
    if(!inherits(spec,'RootPriorSpec'))stop('root_prior must be created by define_root_prior()')
    labels<-vapply(states,.bgb_v12_state_label,character(1))
    if(!length(labels)||anyDuplicated(labels))stop('invalid or duplicate state labels')
    if(!is.null(allowed)) {
        allowed<-as.logical(allowed)
        if(length(allowed)!=length(states)||anyNA(allowed))stop('invalid root allowed mask')
    }
    if(spec$mode=='custom_weights') {
        w<-spec$custom_weights
        if(is.null(names(w))) {
            if(length(w)!=length(states))stop('root weight length does not match states')
        } else {
            if(!setequal(names(w),labels))stop('root state labels must match canonical zero-based range labels')
            w<-w[labels]
        }
        w<-as.numeric(w)
    } else w<-as.numeric(labels!='_')
    if(!is.null(allowed)) {
        if(spec$mode=='custom_weights'&&any(w[!allowed]>0))
            stop('custom root weights assign positive mass to states disallowed in the root stratum')
        w[!allowed]<-0
    }
    if(contract==1L)w[labels=='_']<-0
    total<-sum(w)
    if(any(!is.finite(w))||any(w<0)||!is.finite(total)||total<=0)stop('root prior has no valid probability mass')
    probabilities<-w/total
    weights<-if(contract!=0L||spec$mode=='uniform_nonnull')probabilities else w
    list(likelihood_weights=weights,sampling_probs=probabilities,
         log_weight_scale=log(sum(weights)),state_labels=labels,spec=spec)
}

.bgb_v12_all_states <- function(x) {
    if(!is.null(x$all_geog_states_list_usually_inferred_from_areas_maxareas))
        return(x$all_geog_states_list_usually_inferred_from_areas_maxareas)
    tipranges<-getranges_from_LagrangePHYLIP(x$geogfn)
    areas<-getareas_from_tipranges_object(tipranges)
    rcpp_areas_list_to_states_list(areas=areas,maxareas=x$max_range_size,include_null_range=isTRUE(x$include_null_range))
}

.bgb_v12_layer_mask <- function(inputs,states,i,include_null) {
    allowed<-adjacency<-manual<-rep(TRUE,length(states))
    if(!is.null(inputs$list_of_areas_allowed_mats))
        allowed<-vapply(states,check_if_state_is_allowed,logical(1),inputs$list_of_areas_allowed_mats[[i]])
    if(!is.null(inputs$list_of_areas_adjacency_mats))
        adjacency<-vapply(states,check_if_state_is_allowed_by_adjacency,logical(1),inputs$list_of_areas_adjacency_mats[[i]])
    if(!is.null(inputs$lists_of_states_lists_0based))
        manual<-states %in% inputs$lists_of_states_lists_0based[[i]]
    if(include_null)allowed[1]<-adjacency[1]<-manual[1]<-TRUE
    allowed&adjacency&manual
}

.bgb_v12_validate_scope <- function(x) {
    if(!is.null(x$trait)||!is.null(x$trait_model)||
       length(x$fixnode)||isTRUE(x$use_detection_model)||isTRUE(x$useAmbiguities)||isTRUE(x$rootedge))
        stop('UnsupportedModelContractError: probability contracts exclude trait, detection, root-edge, ambiguity and fixed-node models')
    if(!isTRUE(x$include_null_range))
        stop('UnsupportedModelContractError: probability contracts require include_null_range=TRUE')
    p<-x$BioGeoBEARS_model_object@params_table
    for(key in c('a','b'))if(key %in% rownames(p)) {
        target<-if(key=='b')1 else 0
        if(p[key,'type']!='fixed'||p[key,'init']!=target||p[key,'est']!=target)
            stop('UnsupportedModelContractError: nonstandard rate modifiers a and b await validation')
    }
    for(field in c('list_of_areas_allowed_mats','list_of_areas_adjacency_mats')) {
        mats<-x[[field]]
        if(!is.null(mats))for(m in mats) {
            vals<-as.numeric(as.matrix(m))
            if(any(!is.finite(vals))||any(!vals %in% c(0,1)))
                stop('UnsupportedModelContractError: ',field,' must contain binary 0/1 entries')
        }
    }
    check_pos<-function(key,mats,diag_ignore=FALSE) {
        if(!key %in% rownames(p)||p[key,'type']!='free'||p[key,'min']>=0||is.null(mats))return(invisible(NULL))
        for(m in mats) {
            arr<-as.matrix(m);storage.mode(arr)<-'double'
            if(diag_ignore&&nrow(arr)==ncol(arr))diag(arr)<-1
            if(any(!is.finite(arr))||any(arr<=0))
                stop('free ',key,' with negative bounds requires strictly positive modifier entries')
        }
    }
    check_pos('x',x$list_of_distances_mats,TRUE)
    check_pos('n',x$list_of_envdistances_mats,TRUE)
    check_pos('w',x$list_of_dispersal_multipliers_mats,FALSE)
    if('u' %in% rownames(p)&&p['u','type']=='free'&&p['u','min']<0&&!is.null(x$list_of_area_of_areas)) {
        for(a in x$list_of_area_of_areas) {
            vals<-as.numeric(unlist(a,use.names=FALSE))
            if(any(!is.finite(vals))||any(vals<=0))
                stop('free u with negative bounds requires strictly positive area_of_areas entries')
        }
    }
    if(is.numeric(x$timeperiods)) {
        if(length(x$timeperiods)<2L||is.null(x$tree_sections_list)||is.null(x$master_table))
            stop('UnsupportedModelContractError: stratified probability contracts require section_the_tree(..., make_master_table=TRUE) over >=2 timeperiods')
        if(any(isTRUE(x$master_table$fossils)))
            stop('UnsupportedModelContractError: probability contracts exclude fossil tips')
        orig_rows<-which(x$master_table$piececlass=='orig_tip')
        if(length(orig_rows)&&any(x$master_table$time_bp[orig_rows]>=x$timeperiods[1L]))
            stop('UnsupportedModelContractError: probability contracts require all tips in the present stratum')
        states<-.bgb_v12_all_states(x)
        k<-length(x$timeperiods)
        masks<-lapply(seq_len(k),function(i).bgb_v12_layer_mask(x,states,i,TRUE))
        if(any(vapply(masks,function(m)sum(m)<=1L,logical(1))))
            stop('UnsupportedModelContractError: every stratum must retain at least one occupied geographic state')
        for(i in seq_len(k-1L))if(any(masks[[i+1L]]&!masks[[i]]))
            stop('UnsupportedModelContractError: probability contracts require monotone state-space emergence across time strata')
    }
    invisible(NULL)
}

.bgb_v12_validate_run <- function(x,contract) {
    if(contract==0L&&is.null(x$root_prior))return(invisible(NULL))
    .bgb_v12_validate_scope(x)
    if(contract==2L&&!isTRUE(x$allow_null_tips))stop('unconditional_complete requires allow_null_tips=TRUE')
    if(contract==1L&&isTRUE(x$allow_null_tips))stop('tip_survival requires allow_null_tips=FALSE')
    if(!isTRUE(x$include_null_range))stop('probability contracts require include_null_range=TRUE')
    invisible(NULL)
}

.bgb_v12_original_bears_optim_run <- bears_optim_run
bears_optim_run <- .bgb_v12_original_bears_optim_run
body(bears_optim_run) <- quote({
    # The upstream optimx formatter assumes at least one estimated parameter.
    # A fully fixed model is an evaluation and uses its scalar optim shell.
    if(!any(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type=='free')) {
        skip_optim<-TRUE
        BioGeoBEARS_run_object$use_optimx<-FALSE
    }
    name<-BioGeoBEARS_run_object$likelihood_conditioning
    if(is.null(name))name<-'legacy_unconditioned'
    if(length(name)!=1L||!name %in% .bgb_v12_contract_names)stop('unknown likelihood_conditioning')
    code<-as.integer(match(name,.bgb_v12_contract_names)-1L)
    .bgb_v12_validate_run(BioGeoBEARS_run_object,code)
    old<-as.list(.bgb_v12_runtime)
    on.exit({for(k in names(old)).bgb_v12_runtime[[k]]<-old[[k]]},add=TRUE)
    .bgb_v12_runtime$contract<-code
    .bgb_v12_runtime$root_prior<-BioGeoBEARS_run_object$root_prior
    .bgb_v12_runtime$last_root<-NULL
    .bgb_v12_runtime$last_log_survival<-NULL
    snapshot<-list(tree_text=readLines(BioGeoBEARS_run_object$trfn,warn=FALSE),
                   geography_text=if(!isTRUE(BioGeoBEARS_run_object$use_detection_model))readLines(BioGeoBEARS_run_object$geogfn,warn=FALSE) else NULL)
    result<-.bgb_v12_original_bears_optim_run(BioGeoBEARS_run_object,skip_optim,skip_optim_option)
    if(!identical(snapshot$tree_text,readLines(BioGeoBEARS_run_object$trfn,warn=FALSE))||
       (!is.null(snapshot$geography_text)&&!identical(snapshot$geography_text,readLines(BioGeoBEARS_run_object$geogfn,warn=FALSE))))
        stop('input files changed during fitting')
    if(is.list(result)) {
        result$likelihood_conditioning<-name
        result$root_prior<-BioGeoBEARS_run_object$root_prior
        result$resolved_root_prior<-.bgb_v12_runtime$last_root
        result$log_tip_survival_prob<-.bgb_v12_runtime$last_log_survival
        result$precision_mode<-c('legacy_f32','strict_f64','legacy_fused')[.bgb_precision()+1L]
        result$input_snapshot<-snapshot
    }
    result
})

.bgb_v12_original_define_run <- define_BioGeoBEARS_run
define_BioGeoBEARS_run <- .bgb_v12_original_define_run
body(define_BioGeoBEARS_run) <- substitute({
    value<-do.call(original,mget(names(formals(original)),envir=environment(),inherits=FALSE))
    value$likelihood_conditioning<-'legacy_unconditioned'
    value$root_prior<-NULL
    value$allow_null_tips<-FALSE
    value
},list(original=as.name('.bgb_v12_original_define_run')))

close_prepared_sampler <- function(sampler) {
    if(!inherits(sampler,'PreparedTreeSampler'))stop('expected PreparedTreeSampler')
    .Call(.bgb_c_bgb_sampler_release,sampler$handle)
    invisible(NULL)
}

prepare_tree_sampler <- function(fit_or_model,conditioning='observed',root_distribution='model_prior') {
    if(!root_distribution %in% c('model_prior','empirical_posterior'))stop('unknown root_distribution')
    if(is.null(fit_or_model$inputs)) {
        fit_or_model<-bears_optim_run(fit_or_model,skip_optim=TRUE,skip_optim_option='return_all')
    }
    fit<-fit_or_model;x<-fit$inputs
    if(is.null(x)||is.null(fit$outputs))stop('expected a fitted result or BioGeoBEARS run object')
    .bgb_v12_validate_scope(x)
    if(is.null(fit$input_snapshot))stop('prepared sampling requires a v1.2 fit with an input snapshot; reevaluate the run object')
    if(!identical(fit$input_snapshot$tree_text,readLines(x$trfn,warn=FALSE))||
       !identical(fit$input_snapshot$geography_text,readLines(x$geogfn,warn=FALSE)))
        stop('fitted input files have changed; reevaluate before preparing a sampler')
    if(isTRUE(x$force_sparse))stop('prepared sampling currently requires a dense Q matrix; set force_sparse=FALSE')
    old_precision<-options(RustBEARS.precision=if(fit$precision_mode=='strict_f64')'StrictF64' else 'LegacyF32',
                           RustBEARS.legacy_fma=identical(fit$precision_mode,'legacy_fused'))
    on.exit(options(old_precision),add=TRUE)
    contract<-fit$likelihood_conditioning
    if(is.null(contract))contract<-'legacy_unconditioned'
    if(conditioning=='from_model') {
        if(contract=='legacy_unconditioned')stop('Contract-C requires explicit simulation conditioning')
        conditioning<-contract
    }
    code<-match(conditioning,c('observed','tip_survival','unconditional_complete'))-1L
    if(length(code)!=1L||is.na(code))stop('unsupported simulation conditioning')
    phy<-ape::reorder.phylo(ape::read.tree(x$trfn),'pruningwise')
    if(!is.null(phy$root.edge)&&any(phy$root.edge!=0))stop('UnsupportedModelContractError: root edges await validation')
    nt<-length(phy$tip.label);nn<-nt+phy$Nnode
    if(any(phy$edge.length<x$min_branchlength)||any(table(phy$edge[,1])!=2L))
        stop('UnsupportedModelContractError: prepared sampling requires a binary tree without fossil hooks')
    roots<-setdiff(phy$edge[,1],phy$edge[,2]);roots<-unique(roots)
    if(length(roots)!=1L)stop('invalid root')
    tipranges<-getranges_from_LagrangePHYLIP(x$geogfn)
    matrix_geog<-as.matrix(tipranges@df)
    matrix_geog<-matrix_geog[phy$tip.label,,drop=FALSE]
    area_names<-colnames(matrix_geog);na<-length(area_names);areas_0<-seq_len(na)-1L
    if(code==0L&&contract=='legacy_unconditioned'&&any(rowSums(matrix_geog=='0')==ncol(matrix_geog)))
        stop('UnsupportedModelContractError: Contract-C observed mapping with empty tips is undefined; use unconditional_complete')
    stratified<-is.numeric(x$timeperiods)&&length(x$timeperiods)>=2L
    all_states<-.bgb_v12_all_states(x)
    if(stratified) {
        states<-all_states
        k_layers<-length(x$timeperiods)
        layer_masks<-lapply(seq_len(k_layers),function(i).bgb_v12_layer_mask(x,states,i,isTRUE(x$include_null_range)))
        root_allowed<-layer_masks[[k_layers]]
    } else {
        states<-if(!is.null(x$states_list))x$states_list else {
            keep<-.bgb_v12_layer_mask(x,all_states,1L,isTRUE(x$include_null_range))
            all_states[keep]
        }
        k_layers<-1L
        layer_masks<-list(rep(TRUE,length(states)))
        root_allowed<-NULL
    }
    if(length(states)&&!is.na(states[[1]][1])&&identical(states[[1]],'_'))states[[1]]<-NA
    n<-length(states)
    null<-as.integer(.bgb_v12_state_label(states[[1]])=='_')
    ranges<-vapply(states,function(st)if(!length(st)||anyNA(st)||any(st<0))'_' else paste(area_names[st+1L],collapse=''),character(1))
    p_tab<-fit$outputs@params_table
    if(!.bgb_v12_standard_clado(p_tab[c('s','v','j','y'),'est']))
        stop('UnsupportedModelContractError: sampling supports the six standard cladogenesis models')
    prior<-.bgb_v12_resolve_root(x$root_prior,states,match(contract,.bgb_v12_contract_names)-1L,allowed=root_allowed)
    weights<-if(code==0L)prior$likelihood_weights else prior$sampling_probs
    if(code==1L) {
        weights[ranges=='_']<-0
        if(sum(weights)<=0)stop('tip-survival simulation requires nonempty root mass')
        weights<-weights/sum(weights)
    }
    if(root_distribution=='empirical_posterior') {
        if(code==0L)stop('empirical_posterior is a simulation variant, not an observed BSM root prior')
        weights<-fit$ML_marginal_prob_each_state_at_branch_top_AT_node[roots,]
        if(any(!is.finite(weights))||sum(weights)<=0)stop('invalid empirical root posterior')
        weights<-weights/sum(weights)
        if(code==1L)code<-4L
    }
    tips<-matrix(0,nt,n)
    for(i in seq_len(nt)) {
        row<-as.character(matrix_geog[i,])
        if(any(!row %in% c('0','1','?')))stop('invalid tip geography')
        if(any(row=='?')&&!isTRUE(x$useAmbiguities))stop('ambiguous tips require useAmbiguities=TRUE')
        for(s in seq_len(n)) {
            present<-rep(FALSE,ncol(matrix_geog))
            state<-states[[s]]
            if(!anyNA(state)&&all(state>=0))present[state+1L]<-TRUE
            tips[i,s]<-as.numeric(all(row=='?'|row==as.character(as.integer(present))))
        }
        if(code==0L&&!any(tips[i,layer_masks[[1L]]]>0))
            stop('tip geography for ',phy$tip.label[i],' is outside the allowed present-stratum state space')
    }
    depths<-ape::node.depth.edgelength(phy);ages<-max(depths[seq_len(nt)])-depths
    if(!is.null(x$fossils_older_than))f_cut<-x$fossils_older_than else f_cut<-0.6
    if(any(ages[seq_len(nt)]>f_cut))
        stop('UnsupportedModelContractError: prepared sampling excludes non-contemporaneous fossil tips')
    parent<-rep(NA_integer_,nn);branch<-numeric(nn)
    parent[phy$edge[,2]]<-phy$edge[,1];branch[phy$edge[,2]]<-phy$edge.length
    b<-p_tab['b','est'];branch<-branch^b
    first<-seq.int(1L,nrow(phy$edge),by=2L)
    steps<-cbind(phy$edge[first,1]-1L,phy$edge[first,2]-1L,phy$edge[first+1L,2]-1L)
    child<-which(!is.na(parent))
    if(stratified) {
        times<-as.numeric(x$timeperiods)
        if(times[length(times)]<=max(ages))times[length(times)]<-max(ages)+1
        layer_at_node<-findInterval(ages,times)
        seg_list<-lapply(child,function(ch) {
            old<-ages[parent[ch]];young<-ages[ch]
            inner<-rev(times[times>young&times<old])
            cuts<-c(old,inner,young)
            mids<-0.5*(cuts[-length(cuts)]+cuts[-1L])
            cbind(ch-1L,findInterval(mids,times),cuts[-length(cuts)]-cuts[-1L],cuts[-length(cuts)])
        })
        segments<-do.call(rbind,seg_list)
    } else {
        layer_at_node<-rep(0L,nn)
        segments<-cbind(child-1L,0,branch[child],ages[parent[child]])
    }
    ones<-matrix(1,na,na)
    get_mat<-function(field,idx,diag_one=FALSE) {
        val<-x[[field]]
        m<-if(is.null(val))ones else as.matrix(val[[idx]])
        storage.mode(m)<-'double'
        if(diag_one)diag(m)[is.na(diag(m))|diag(m)==0]<-1
        m
    }
    layer_buffers<-list()
    for(i in seq_len(k_layers)) {
        mask<-layer_masks[[i]]
        cols<-which(mask)
        local_states<-states[mask]
        dist_m<-get_mat('list_of_distances_mats',i,TRUE)
        env_m<-get_mat('list_of_envdistances_mats',i,TRUE)
        man_m<-get_mat('list_of_dispersal_multipliers_mats',i,FALSE)
        area_v<-if(is.null(x$list_of_area_of_areas))rep(1,na) else as.numeric(unlist(x$list_of_area_of_areas[[i]],use.names=FALSE))
        dispersal<-(dist_m^p_tab['x','est'])*(env_m^p_tab['n','est'])*(man_m^p_tab['w','est'])
        diag(dispersal)[is.na(diag(dispersal))]<-1
        elist<-(area_v^p_tab['u','est'])*p_tab['e','est']
        dmat_d<-dispersal*p_tab['d','est']
        amat_a<-dispersal*p_tab['a','est']
        if(any(!is.finite(dispersal))||any(dispersal<0)||any(!is.finite(elist))||any(elist<0))
            stop('spatial/environmental modifier exponents produced invalid rate multipliers')
        old_defer<-.bgb_m2_runtime$defer_model;.bgb_m2_runtime$defer_model<-FALSE
        q_local<-tryCatch(rcpp_states_list_to_DEmat(areas_list=areas_0,states_list=local_states,dmat=dmat_d,elist=elist,
            amat=amat_a,include_null_range=isTRUE(x$include_null_range),normalize_TF=TRUE,makeCOO_TF=FALSE),
            finally={.bgb_m2_runtime$defer_model<-old_defer})
        q<-matrix(0,n,n)
        q[cols,cols]<-q_local
        clado_states<-if(null)local_states[-1L] else local_states
        inp<-get_spPmat_inputs_from_BGB(x,local_states,dispersal)
        inp$s<-p_tab['s','est'];inp$v<-p_tab['v','est'];inp$j<-p_tab['j','est'];inp$y<-p_tab['y','est']
        cached<-.bgb_m2_priors(max(lengths(clado_states)),inp)
        ones_c<-rep(1,length(clado_states))
        coo<-rcpp_calc_anclikes_sp_COOweights_faster(ones_c,ones_c,clado_states,inp$s,inp$v,inp$j,inp$y,
            dispersal,cached$priors[[1]],cached$priors[[2]],cached$priors[[3]],cached$priors[[4]],cached$maxmin,printmat=FALSE)
        idx_map<-cols[(null+1L):length(cols)]-1L
        events<-cbind(idx_map[coo[[1]]+1L],idx_map[coo[[2]]+1L],idx_map[coo[[3]]+1L],coo[[4]])
        layer_buffers<-c(layer_buffers,list(q,as.numeric(mask),as.numeric(t(events)),dispersal))
    }
    buffers<-c(list(c(nt,roots-1L)),.bgb_state_buffers(states),
        list(tips,as.numeric(t(steps)),ages,as.numeric(layer_at_node),rep(0,nn),matrix(1,nn,n),
             as.numeric(t(segments))),layer_buffers)
    handle<-.Call(.bgb_c_bgb_sampler_prepare,buffers,weights,as.integer(code))
    structure(list(handle=handle,states=states,ranges=ranges,areas=area_names,
        tip_labels=phy$tip.label,num_nodes=nn,conditioning=conditioning,
        tree_topology_hash=.bgb_u_tree_key(phy)$topology_hash,model_id=1L,
        likelihood_conditioning=contract,root_distribution=root_distribution,
        root_sampling_probs=weights/sum(weights),precision_mode=fit$precision_mode),class='PreparedTreeSampler')
}

simulate_datasets <- function(fit_or_model,n_reps=100L,conditioning='from_model',
        root_distribution='model_prior',num_threads=1L,seed=12345,start=1,return_node_states=FALSE) {
    owned<-!inherits(fit_or_model,'PreparedTreeSampler')
    sampler<-if(owned)prepare_tree_sampler(fit_or_model,conditioning,root_distribution) else fit_or_model
    if(owned)on.exit(close_prepared_sampler(sampler),add=TRUE)
    if(sampler$conditioning=='observed')stop('tip simulation requires a forward sampling condition')
    raw<-.Call(.bgb_c_bgb_sampler_simulate,sampler$handle,start,n_reps,seed,num_threads)
    tips<-matrix(raw[[1]],nrow=n_reps,byrow=TRUE,dimnames=list(NULL,sampler$tip_labels))
    structure(list(tip_states=tips,node_states=if(return_node_states)matrix(raw[[2]],nrow=n_reps,byrow=TRUE) else NULL,
        states=sampler$states,ranges=sampler$ranges,areas=sampler$areas,tip_labels=sampler$tip_labels,
        conditioning=sampler$conditioning,root_distribution=sampler$root_distribution,
        root_sampling_probs=sampler$root_sampling_probs,seed=seed,start=start),class='SimulatedBatch')
}

.bgb_v12_event_frame <- function(raw,columns,sampler) {
    out<-as.data.frame(matrix(raw,ncol=length(columns),byrow=TRUE));names(out)<-columns
    out$event_type<-strsplit('deajysvt','')[[1]][out$event_type+1L]
    for(key in intersect(c('node','left_node','right_node','stratum'),names(out)))out[[key]]<-out[[key]]+1L
    for(key in intersect(c('source_area','destination_area'),names(out)))out[[key]]<-c(NA_character_,sampler$areas)[out[[key]]+2L]
    for(key in grep('_state$',names(out),value=TRUE))out[[sub('_state$','_range',key)]]<-sampler$ranges[out[[key]]+1L]
    out
}

draw_prepared_bsm <- function(sampler,num_maps=100L,start=1,seed=12345,num_threads=1L,rejection=FALSE) {
    if(!inherits(sampler,'PreparedTreeSampler')||sampler$conditioning!='observed')stop('expected an observed PreparedTreeSampler')
    raw<-.Call(.bgb_c_bgb_sampler_draw,sampler$handle,start,num_maps,seed,num_threads,as.integer(rejection))
    clado<-.bgb_v12_event_frame(raw[[1]],c('map','node','left_node','right_node','ancestral_state','left_state','right_state',
        'time_bp','stratum','event_type','source_area','destination_area'),sampler)
    ana<-.bgb_v12_event_frame(raw[[2]],c('map','node','stratum','from_state','to_state','time_bp','segment_time',
        'event_type','source_area','destination_area'),sampler)
    metadata<-data.frame(tree_id=1L,tree_topology_hash=sampler$tree_topology_hash,model_id=sampler$model_id,
        map=as.double(start+seq_len(num_maps)-1L),map_id=as.double(start+seq_len(num_maps)-1L),stratum_weight=1,precision_mode=sampler$precision_mode,
        likelihood_conditioning=sampler$likelihood_conditioning,stringsAsFactors=FALSE)
    annotate<-function(events) {
        index<-match(events$map,metadata$map_id)
        for(key in names(metadata))events[[key]]<-metadata[[key]][index]
        events
    }
    structure(list(clado_events_df=annotate(clado),ana_events_df=annotate(ana),num_maps=num_maps,start=start,seed=seed,
        map_metadata_df=metadata,tree_topology_hash=sampler$tree_topology_hash,
        likelihood_conditioning=sampler$likelihood_conditioning,areas=sampler$areas,
        precision_mode=sampler$precision_mode),class='BSMResult')
}
