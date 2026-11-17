# GPL >= 3. Uncertainty workflows over the native probability contracts.
.bgb_u_params <- function(res) stats::setNames(res$outputs@params_table$est,rownames(res$outputs@params_table))
.bgb_u_free <- function(res) rownames(res$outputs@params_table)[res$outputs@params_table$type=='free']
.bgb_u_converged <- function(res) {
    if(!is.finite(res$total_loglikelihood))return(FALSE)
    if(!length(.bgb_u_free(res)))return(TRUE)
    o<-res$optim_result
    if(!is.null(o$convergence))return(identical(as.integer(o$convergence),0L))
    if(!is.null(o$convcode)) {
        if(is.null(o$value)||!any(is.finite(o$value)))return(FALSE)
        winner<-which.max(ifelse(is.finite(o$value),o$value,-Inf))
        return(isTRUE(o$convcode[winner]==0))
    }
    FALSE
}
.bgb_u_assert_fit <- function(res) {
    if(!is.list(res)||is.null(res$outputs)||is.null(res$inputs)||is.null(res$likelihood_conditioning)||
       !is.finite(res$total_loglikelihood))stop('expected a finite v1.2 fitted result')
    if(is.null(res$input_snapshot))stop('fit has no frozen input snapshot')
    invisible(res)
}
.bgb_u_root <- function(res) {
    if(!is.null(res$resolved_root_prior))return(res$resolved_root_prior$likelihood_weights)
    n<-ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
    if(is.null(n))stop('fit needs ancestral-state output')
    w<-c(if(isTRUE(res$inputs$include_null_range))0,rep(1,n-as.integer(isTRUE(res$inputs$include_null_range))))
    if(is.numeric(res$inputs$timeperiods)&&length(res$inputs$timeperiods)>1L) {
        states<-.bgb_v12_all_states(res$inputs)
        if(length(states)==n) {
            mask<-.bgb_v12_layer_mask(res$inputs,states,length(res$inputs$timeperiods),isTRUE(res$inputs$include_null_range))
            w[!mask]<-0
        }
    }
    w
}
.bgb_u_compatible <- function(fits) {
    if(!length(fits))stop('at least one fit is required')
    lapply(fits,.bgb_u_assert_fit);first<-fits[[1]]
    for(f in fits) {
        if(!identical(.bgb_u_data_signature(f),.bgb_u_data_signature(first))||
           !identical(f$likelihood_conditioning,first$likelihood_conditioning)||
           !identical(f$precision_mode,first$precision_mode)||
           !identical(.bgb_u_root(f),.bgb_u_root(first)))
            stop('ContractMismatchError: fits must share observations, root weights, precision and likelihood contract')
    }
    invisible(fits)
}
.bgb_u_data_signature <- function(res) {
    phy<-ape::read.tree(text=paste(res$input_snapshot$tree_text,collapse='\n'))
    nt<-length(phy$tip.label);nn<-nt+phy$Nnode;clades<-vector('list',nn)
    for(i in seq_len(nt))clades[[i]]<-phy$tip.label[i]
    descendants<-function(node) {
        if(!is.null(clades[[node]]))return(clades[[node]])
        result<-sort(unlist(lapply(phy$edge[phy$edge[,1]==node,2],descendants),use.names=FALSE))
        clades[[node]]<<-result;result
    }
    keys<-vapply(seq_len(nn),function(node)paste(vapply(descendants(node),function(label)
        paste0(nchar(label,type='bytes'),':',label),character(1)),collapse='|'),character(1))
    edge<-lapply(seq_len(nrow(phy$edge)),function(i)list(parent=keys[phy$edge[i,1]],child=keys[phy$edge[i,2]],length=phy$edge.length[i]))
    edge<-edge[order(vapply(edge,function(x)x$child,character(1)))]
    # Source files may permute data rows without changing labelled observations.
    lines<-trimws(res$input_snapshot$geography_text);lines<-lines[nzchar(lines)]
    observations<-lapply(lines[-1L],function(line)strsplit(line,'[[:space:]]+')[[1]])
    observations<-observations[order(vapply(observations,`[[`,character(1),1L))]
    list(tips=sort(phy$tip.label),edges=edge,root_edge=phy$root.edge,geography_header=lines[1L],observations=observations)
}
.bgb_u_refit <- function(res,fixed=list(),start=NULL,geogfn=NULL,contract=NULL,optimize=TRUE,num_threads=1L,calc_ancprobs=TRUE) {
    x<-res$inputs;x$BioGeoBEARS_model_object<-res$outputs
    p<-x$BioGeoBEARS_model_object@params_table
    p$init<-p$est
    if(length(start))for(k in intersect(names(start),rownames(p)))if(p[k,'type']=='free')
        p[k,c('init','est')]<-max(p[k,'min'],min(p[k,'max'],start[[k]]))
    if(length(fixed))for(k in names(fixed)) {
        if(!k %in% rownames(p)||!is.finite(fixed[[k]]))stop('invalid fixed profile parameter')
        p[k,'type']<-'fixed';p[k,c('init','est')]<-fixed[[k]]
    }
    x$BioGeoBEARS_model_object@params_table<-p
    if(!is.null(geogfn))x$geogfn<-geogfn
    if(!is.null(contract)) {
        x$likelihood_conditioning<-contract
        x$allow_null_tips<-contract=='unconditional_complete'
        if(is.null(x$root_prior))x$root_prior<-define_root_prior(if(contract=='legacy_unconditioned')'legacy_flat_weights' else 'uniform_nonnull')
    }
    x$num_cores_to_use<-as.integer(num_threads);x$printlevel<--1L;x$print_optim<-FALSE
    x$calc_ancprobs<-isTRUE(calc_ancprobs);x$use_optimx<-FALSE
    old<-options(RustBEARS.precision=if(res$precision_mode=='strict_f64')'StrictF64' else 'LegacyF32')
    on.exit(options(old),add=TRUE)
    result<-NULL
    invisible(utils::capture.output(result<-bears_optim_run(x,
        skip_optim=!optimize||!any(p$type=='free'),skip_optim_option='return_all')))
    # The upstream scalar-only branch wraps its score as the first unnamed list
    # element when adding inputs/outputs. Normalize that shell for callers.
    if(!isTRUE(calc_ancprobs)&&is.list(result)&&is.null(result$total_loglikelihood)&&
       length(result[[1]])==1L&&is.numeric(result[[1]]))result$total_loglikelihood<-as.numeric(result[[1]])
    result
}
.bgb_u_prepare_comparison <- function(a,b,target='match_fits',refit_observed=FALSE) {
    .bgb_u_compatible(list(a,b));changed<-FALSE
    if(target=='match_fits')target<-a$likelihood_conditioning
    if(!target %in% .bgb_v12_contract_names)stop('unknown target_conditioning')
    if(target!=a$likelihood_conditioning) {
        if(!isTRUE(refit_observed))stop('ContractMismatchError: changing likelihood contract requires refit_observed=TRUE')
        a<-.bgb_u_refit(a,contract=target);b<-.bgb_u_refit(b,contract=target);changed<-TRUE
        .bgb_u_compatible(list(a,b))
    }
    list(null=a,alt=b,recomputed_observed_under_new_contract=changed)
}
.bgb_u_json_string <- function(x) {
    if(length(x)!=1L||is.na(x)||!is.character(x))stop('task string must be a scalar')
    esc<-vapply(utf8ToInt(enc2utf8(x)),function(k) {
        if(k==34L)return('\\"')
        if(k==92L)return('\\\\')
        simple<-c('8'='\\b','9'='\\t','10'='\\n','12'='\\f','13'='\\r')
        if(as.character(k) %in% names(simple))return(simple[[as.character(k)]])
        if(k<32L||k>126L) {
            if(k<=65535L)return(sprintf('\\u%04x',k))
            z<-k-65536L;return(paste0(sprintf('\\u%04x',55296L+z%/%1024L),sprintf('\\u%04x',56320L+z%%1024L)))
        }
        intToUtf8(k)
    },character(1))
    paste0('"',paste0(esc,collapse=''),'"')
}
.bgb_u_json_utf8_string <- function(x) {
    esc<-vapply(utf8ToInt(enc2utf8(x)),function(k) {
        if(k==34L)return('\\"')
        if(k==92L)return('\\\\')
        simple<-c('8'='\\b','9'='\\t','10'='\\n','12'='\\f','13'='\\r')
        if(as.character(k) %in% names(simple))return(simple[[as.character(k)]])
        if(k<32L)return(sprintf('\\u%04x',k))
        intToUtf8(k)
    },character(1))
    paste0('"',paste0(esc,collapse=''),'"')
}
.bgb_u_topology_hash <- function(tree,members=NULL) {
    if(is.null(members))members<-.bgb_u_tree_members(tree)
    members<-lapply(members,function(x)sort(enc2utf8(x),method='radix'))
    columns<-lapply(seq_len(max(lengths(members))),function(k)vapply(members,function(x)
        if(length(x)>=k)x[k] else '',character(1)))
    ordered<-do.call(order,c(columns,list(method='radix')))
    arrays<-vapply(members[ordered],function(x)paste0('[',paste0(vapply(x,.bgb_u_json_utf8_string,character(1)),collapse=','),']'),character(1))
    text<-paste0('[',paste0(arrays,collapse=','),']')
    digest::digest(charToRaw(enc2utf8(text)),algo='sha256',serialize=FALSE)
}
.bgb_u_task_seed <- function(seed,purpose,model_id,tree_id,outer_rep=0,inner_rep=0,role='',layer_id=0,step_id=0) {
    for(value in list(purpose,model_id,tree_id))if(length(value)!=1L||!is.character(value)||is.na(value)||!nzchar(value))
        stop('purpose, model_id and tree_id must be nonempty stable string identifiers')
    task<-list(inner_rep=inner_rep,layer_id=layer_id,model_id=model_id,outer_rep=outer_rep,
               purpose=purpose,role=role,step_id=step_id,tree_id=tree_id)
    fields<-vapply(names(task),function(k) {
        value<-task[[k]]
        if(k %in% c('inner_rep','outer_rep','layer_id','step_id')) {
            # The bridge accepts exact decimal/hex strings for all u64 values;
            # doubles remain restricted to their exact integer range.
            encoded<-.Call(.bgb_c_bgb_seed_decimal,value)
        } else encoded<-.bgb_u_json_string(value)
        paste0('"',k,'":',encoded)
    },character(1))
    payload<-paste0('{"master_seed":',.Call(.bgb_c_bgb_seed_decimal,seed),',"task":{',paste0(fields,collapse=','),'}}')
    bytes<-c(charToRaw('RustBEARS/uncertainty/task-seed/v1'),as.raw(0),charToRaw(payload))
    hash<-digest::digest(bytes,algo='sha256',serialize=FALSE)
    list(seed=paste0('0x',substr(hash,1,16)),fingerprint=hash,payload=payload)
}
.bgb_u_checkpoint <- function(path,manifest) {
    native<-getLoadedDLLs()[['RustBEARS']][['path']]
    manifest$native_sha256<-digest::digest(file=native,algo='sha256')
    manifest$r_version<-as.character(getRversion())
    code_names<-c('profile_ci','profile_surface','bootstrap_lrt','model_average','fit_over_trees','bsm_diagnostics',
        'simulate_datasets','prepare_tree_sampler','.bgb_u_refit','.bgb_u_task_seed','.bgb_s_fit','.bgb_s_point','.bgb_s_pair','.bgb_s_fit_snapshot',
        '.bgb_u_combine_maps','.bgb_u_scalar_summary','.bgb_u_model_definition','.bgb_u_topology_hash')
    manifest$r_interface_sha256<-digest::digest(stats::setNames(lapply(code_names,function(name) {
        fn<-get0(name,mode='function',inherits=TRUE)
        if(is.null(fn))NULL else list(formals=deparse(formals(fn)),body=deparse(body(fn)))
    }),code_names),algo='sha256')
    fingerprint<-digest::digest(manifest,algo='sha256')
    if(is.null(path))return(list(path=NULL,fingerprint=fingerprint,manifest=manifest,records=list()))
    if(file.exists(path)) {
        saved<-readRDS(path)
        if(!identical(saved$fingerprint,fingerprint))stop('CheckpointMismatchError: task configuration or native library changed')
        return(saved)
    }
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    state<-list(path=path,fingerprint=fingerprint,manifest=manifest,records=list())
    .bgb_u_save_checkpoint(state);state
}
.bgb_u_save_checkpoint <- function(state) {
    if(is.null(state$path))return(invisible(NULL))
    part<-tempfile('.checkpoint-',tmpdir=dirname(state$path))
    on.exit(unlink(part),add=TRUE);saveRDS(state,part,version=3)
    if(!file.rename(part,state$path))stop('atomic checkpoint replacement failed')
    manifest_path<-paste0(state$path,'.manifest.json')
    manifest_part<-tempfile('.manifest-',tmpdir=dirname(state$path))
    on.exit(unlink(manifest_part),add=TRUE)
    header<-list(format='RustBEARS.R.checkpoint.v1',task_fingerprint=state$fingerprint,
        payload_sha256=digest::digest(file=state$path,algo='sha256'),manifest=state$manifest)
    writeLines(.bgb_u_manifest_json(header),manifest_part,useBytes=TRUE)
    if(!file.rename(manifest_part,manifest_path))stop('atomic checkpoint manifest replacement failed')
    invisible(NULL)
}
.bgb_u_manifest_json <- function(value) {
    if(is.null(value))return('null')
    if(isS4(value))return(.bgb_u_manifest_json(list(s4_class=class(value),slots=stats::setNames(
        lapply(methods::slotNames(value),function(k)methods::slot(value,k)),methods::slotNames(value)))))
    if(is.data.frame(value))return(.bgb_u_manifest_json(list(data_frame=as.list(value),row_names=rownames(value))))
    if(is.list(value)) {
        fields<-vapply(value,.bgb_u_manifest_json,character(1))
        if(is.null(names(value)))return(paste0('[',paste0(fields,collapse=','),']'))
        return(paste0('{',paste0(vapply(names(value),.bgb_u_json_string,character(1)),':',fields,collapse=','),'}'))
    }
    # Hexadecimal double strings preserve every finite binary64 value and avoid
    # non-standard JSON NaN/Infinity tokens. RDS remains the complete payload.
    vals<-if(is.double(value))sprintf('%a',value) else if(is.raw(value))sprintf('%02x',as.integer(value)) else as.character(value)
    vals[is.na(vals)]<-'NA'
    encoded<-vapply(vals,.bgb_u_json_string,character(1))
    paste0('{"type":',.bgb_u_json_string(typeof(value)),',"values":[',paste0(encoded,collapse=','),'],"names":',
        if(is.null(names(value)))'null' else paste0('[',paste0(vapply(names(value),.bgb_u_json_string,character(1)),collapse=','),']'),
        ',"dim":',if(is.null(dim(value)))'null' else paste0('[',paste0(dim(value),collapse=','),']'),'}')
}
.bgb_u_write_simulated <- function(batch,rep,filename) {
    rows<-vapply(batch$tip_states[rep,],function(index) {
        state<-batch$states[[index+1L]];bits<-rep('0',length(batch$areas))
        if(!anyNA(state)&&all(state>=0))bits[state+1L]<-'1'
        paste0(bits,collapse='')
    },character(1))
    writeLines(c(sprintf('%d %d (%s)',length(rows),length(batch$areas),paste(batch$areas,collapse=' ')),
        paste(batch$tip_labels,rows)),filename)
    filename
}
simulate_biogeog_datasets <- function(res,n_reps=100L,conditioning=c('from_model','tip_survival','unconditional_complete'),
        root_distribution=c('model_prior','empirical_posterior'),seedval=12345L,num_threads=getOption('RustBEARS.threads',1L)) {
    simulate_datasets(res,n_reps=n_reps,conditioning=match.arg(conditioning),root_distribution=match.arg(root_distribution),
        seed=seedval,num_threads=num_threads,return_node_states=TRUE)
}
