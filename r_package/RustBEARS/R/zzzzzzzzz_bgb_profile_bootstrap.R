# GPL >= 3. Native R profile and parametric-bootstrap workflows; no Python bridge.
.bgb_s_positive <- function(x,name) {
    if(length(x)!=1L||!is.numeric(x)||!is.finite(x)||x<1||x!=floor(x))stop(name,' must be a positive integer')
    as.integer(x)
}
.bgb_s_name <- function(res) {
    p<-res$outputs@params_table
    base<-if(p['v','type']=='fixed'&&p['v','est']==0)'BAYAREALIKE' else
        if(p['s','type']=='fixed'&&p['s','est']==0)'DIVALIKE' else 'DEC'
    paste0(base,if(p['j','type']=='free')'+J' else '',if(p['x','type']=='free')'+x' else '',
           if(p['w','type']=='free')'+w' else '')
}
.bgb_s_check <- function(res) {
    .bgb_u_assert_fit(res);.bgb_v12_validate_scope(res$inputs)
    if(!identical(res$input_snapshot$tree_text,readLines(res$inputs$trfn,warn=FALSE))||
       !identical(res$input_snapshot$geography_text,readLines(res$inputs$geogfn,warn=FALSE)))
        stop('fitted input files changed; refit before uncertainty analysis')
    invisible(res)
}
.bgb_s_manifest <- function(res) {
    files<-lapply(res$inputs[intersect(names(res$inputs),c('trfn','geogfn','distsfn','dispersal_multipliers_fn',
        'timesfn','envdistsfn','area_of_areas_fn','areas_allowed_fn','areas_adjacency_fn'))],function(path) {
            if(length(path)==1L&&!is.na(path)&&nzchar(path)&&file.exists(path))readBin(path,'raw',n=file.info(path)$size) else path
        })
    list(version=as.character(utils::packageVersion('RustBEARS')),abi=1L,inputs=res$inputs,
         files=files,snapshot=res$input_snapshot,parameters=res$outputs@params_table,
         root=.bgb_u_root(res),contract=res$likelihood_conditioning,precision=res$precision_mode,
         objective=res$total_loglikelihood)
}
.bgb_s_store <- function(path,manifest,resume) {
    if(isTRUE(resume)&&(is.null(path)||!file.exists(path)))stop('resume requires an existing checkpoint_file')
    if(!isTRUE(resume)&&!is.null(path)&&file.exists(path))stop('checkpoint exists; use resume=TRUE')
    .bgb_u_checkpoint(path,manifest)
}
.bgb_s_condition <- function(res,conditioning) {
    if(conditioning=='from_model') {
        if(res$likelihood_conditioning=='legacy_unconditioned')stop('Contract-C requires explicit simulation_conditioning')
        conditioning<-res$likelihood_conditioning
    }
    if(!conditioning %in% c('tip_survival','unconditional_complete'))stop('unsupported simulation_conditioning')
    if(res$likelihood_conditioning!='legacy_unconditioned'&&conditioning!=res$likelihood_conditioning)
        stop('ContractMismatchError: simulation must match the fitted probability contract')
    conditioning
}
.bgb_s_tasks <- function(seed,tasks) {
    result<-lapply(tasks,function(t)do.call(.bgb_u_task_seed,c(list(seed=seed),t)))
    hashes<-vapply(result,`[[`,character(1),'fingerprint');seeds<-vapply(result,`[[`,character(1),'seed')
    if(anyDuplicated(hashes))stop('duplicate random task identity')
    if(anyDuplicated(seeds))stop('64-bit task seed collision')
    Map(function(value,identity){value$identity<-identity;value},result,tasks)
}
.bgb_s_modified <- function(res,fixed=list(),start=NULL,geogfn=NULL) {
    p<-res$outputs@params_table;p$init<-p$est
    for(k in intersect(names(start),rownames(p)))if(p[k,'type']=='free')
        p[k,c('init','est')]<-max(p[k,'min'],min(p[k,'max'],start[[k]]))
    for(k in names(fixed)) {
        if(!k %in% rownames(p)||length(fixed[[k]])!=1L||!is.finite(fixed[[k]]))stop('invalid fixed parameter')
        p[k,'type']<-'fixed';p[k,c('init','est')]<-fixed[[k]]
    }
    res$outputs@params_table<-p
    res$outputs<-calc_linked_params_BioGeoBEARS_model_object(res$outputs,update_init=TRUE)
    res$inputs$BioGeoBEARS_model_object<-res$outputs
    if(!is.null(geogfn)) {
        res$inputs$geogfn<-geogfn
        if(res$likelihood_conditioning!='tip_survival')res$inputs$allow_null_tips<-TRUE
    }
    res
}
.bgb_s_score_factory <- function(res,num_threads=1L) {
    x<-res$inputs;x$BioGeoBEARS_model_object<-res$outputs;x$printlevel<--1L;x$print_optim<-FALSE
    phy<-ape::reorder.phylo(ape::read.tree(x$trfn),'pruningwise')
    mats<-get_Qmat_COOmat_from_BioGeoBEARS_run_object(x,BioGeoBEARS_model_object=res$outputs,
        max_range_size=x$max_range_size,include_null_range=x$include_null_range)
    states<-mats$states_list;areas<-seq_len(length(getareas_from_tipranges_object(getranges_from_LagrangePHYLIP(x$geogfn))))-1L
    tips<-tipranges_to_tip_condlikes_of_data_on_each_state(getranges_from_LagrangePHYLIP(x$geogfn),phy,
        states_list=states,maxareas=x$max_range_size,include_null_range=x$include_null_range,
        useAmbiguities=FALSE,trait_as_tip_condlikes=NULL,allow_null_tips=isTRUE(x$allow_null_tips))
    function(params) {
        old<-as.list(.bgb_v12_runtime);old_threads<-.bgb_m2_runtime$threads
        old_precision<-options(RustBEARS.precision=if(res$precision_mode=='strict_f64')'StrictF64' else 'LegacyF32')
        on.exit({for(k in names(old)).bgb_v12_runtime[[k]]<-old[[k]]
            .bgb_m2_runtime$threads<-old_threads;options(old_precision)},add=TRUE)
        .bgb_v12_runtime$contract<-match(res$likelihood_conditioning,.bgb_v12_contract_names)-1L
        .bgb_v12_runtime$root_prior<-if(is.null(x$root_prior))
            define_root_prior(if(res$likelihood_conditioning=='legacy_unconditioned')'legacy_flat_weights' else 'uniform_nonnull') else x$root_prior
        .bgb_m2_runtime$threads<-num_threads
        calc_loglike_for_optim(params,x,phy,tips,print_optim=FALSE,areas_list=areas,states_list=states,
            force_sparse=isTRUE(x$force_sparse),cluster_already_open=NULL,return_what='loglike',calc_ancprobs=FALSE)
    }
}
.bgb_s_fit <- function(res,fixed=list(),start=.bgb_u_params(res),geogfn=NULL,num_threads=1L,maxit=500L,ancestors=TRUE) {
    modified<-.bgb_s_modified(res,fixed,start,geogfn);p<-modified$outputs@params_table
    free<-rownames(p)[p$type=='free'];score<-.bgb_s_score_factory(modified,num_threads)
    theta<-stats::setNames(p[free,'est'],free)
    if(length(free)) {
        search<-.bgb_composite_search(theta,function(params,...)score(params),
            list(BioGeoBEARS_run_object=modified$inputs,phy=ape::read.tree(modified$inputs$trfn)),
            p[free,'min'],p[free,'max'],maxit=maxit)
        if(search$convergence!=0L||!is.finite(search$value)) {
            raw<-modified
            raw$outputs<-calc_linked_params_BioGeoBEARS_model_object(
                params_into_BioGeoBEARS_model_object(modified$outputs,search$par))
            snapshot<-list(obj=search$value,converged=FALSE,params=.bgb_u_params(raw),
                n_evals=search$count,message=paste('constrained optimizer did not converge; code',search$convergence))
            stop(structure(list(message='constrained optimizer did not converge',call=NULL,fit_snapshot=snapshot),
                class=c('rustbears_fit_failure','error','condition')))
        }
        theta<-search$par
    } else search<-list(value=score(numeric()),convergence=0L,count=1L)
    if(!is.finite(search$value))stop('constrained optimum is not finite')
    if(ancestors)out<-.bgb_u_refit(modified,start=theta,optimize=FALSE,num_threads=num_threads,calc_ancprobs=TRUE)
    else {
        out<-modified;out$outputs<-calc_linked_params_BioGeoBEARS_model_object(
            params_into_BioGeoBEARS_model_object(modified$outputs,theta))
        out$inputs$BioGeoBEARS_model_object<-out$outputs;out$total_loglikelihood<-score(theta)
        out$input_snapshot<-list(tree_text=readLines(out$inputs$trfn,warn=FALSE),geography_text=readLines(out$inputs$geogfn,warn=FALSE))
        out[grep('^ML_marginal|^relative_probs|^computed_likelihoods|^condlikes',names(out))]<-NULL
        out$log_tip_survival_prob<-NULL
    }
    out$optim_result<-list(convergence=search$convergence,counts=stats::setNames(search$count,'function'),par=theta,value=out$total_loglikelihood)
    out$n_evals_uncertainty<-search$count
    out
}
.bgb_s_point <- function(res,fixed,num_threads=1L,maxit=500L) {
    modified<-.bgb_s_modified(res,fixed);p<-modified$outputs@params_table;free<-rownames(p)[p$type=='free']
    interior<-stats::setNames((p[free,'min']+p[free,'max'])/2,free)
    probe<-.bgb_s_score_factory(modified,num_threads)(interior)
    if(identical(as.numeric(probe),-Inf))return(list(profile_obj=-Inf,converged=TRUE,structural_zero=TRUE,
        params=c(.bgb_u_params(modified)[setdiff(rownames(p),free)],interior),n_evals=1L,fit=NULL))
    if(!is.finite(probe))stop('nonfinite numerical support probe')
    fit<-tryCatch(.bgb_s_fit(modified,num_threads=num_threads,maxit=maxit,ancestors=FALSE),error=function(e)
        .bgb_s_fit(modified,start=interior,num_threads=num_threads,maxit=2L*maxit,ancestors=FALSE))
    list(profile_obj=fit$total_loglikelihood,converged=TRUE,structural_zero=FALSE,
         params=.bgb_u_params(fit),n_evals=fit$n_evals_uncertainty+1L,fit=fit)
}
.bgb_s_failure_intervals <- function(r,failed,planned,confidence=.95) {
    if(planned<1||r<0||failed<0||r+failed>planned)stop('invalid failure counts')
    alpha<-1-confidence;upper<-r+failed
    list(failure_bounds=c((r+1)/(planned+1),(upper+1)/(planned+1)),
         mc_ci95=c(if(r==0)0 else stats::qbeta(alpha/2,r,planned-r+1),
                   if(upper==planned)1 else stats::qbeta(1-alpha/2,upper+1,planned-upper)))
}
.bgb_s_nested <- function(null,alt) {
    n<-null$outputs@params_table;a<-alt$outputs@params_table
    if(!identical(rownames(n),rownames(a)))stop('models have different parameter sets')
    extra<-character();defaults<-c(d=1e-12,e=1e-12,a=1e-12,j=1e-5)
    for(k in rownames(n)) {
        nt<-n[k,'type'];at<-a[k,'type']
        if(!nt %in% c('fixed','free')||!at %in% c('fixed','free')) {
            if(nt!=at)stop('non-nested parameter link: ',k)
        } else if(at=='fixed') {
            if(nt!='fixed'||n[k,'est']!=a[k,'est'])stop('alternative fixes null parameter: ',k)
        } else if(nt=='free') {
            if(n[k,'min']<a[k,'min']||n[k,'max']>a[k,'max'])stop('alternative bounds do not contain null bounds: ',k)
        } else {
            value<-n[k,'est']
            if(value==0&&k %in% names(defaults)&&a[k,'min']>0) {
                if(a[k,'min']!=defaults[[k]])stop('custom positive bound excludes null physical zero: ',k)
                a[k,'min']<-0
            } else if(value<a[k,'min']||value>a[k,'max'])stop('null value outside alternative box: ',k)
            if(a[k,'min']<a[k,'max'])extra<-c(extra,k)
        }
    }
    if(!length(extra))stop('alternative must add a free parameter')
    alt$outputs@params_table<-a;alt$inputs$BioGeoBEARS_model_object<-alt$outputs
    list(alt=alt,extra=extra)
}
.bgb_s_fit_snapshot <- function(fit) {
    list(obj=fit$total_loglikelihood,converged=.bgb_u_converged(fit),params=.bgb_u_params(fit),
        n_evals=if(is.null(fit$n_evals_uncertainty))NA_integer_ else fit$n_evals_uncertainty,
        message=if(is.null(fit$optim_result$message))'' else as.character(fit$optim_result$message))
}
.bgb_s_pair <- function(null,alt,num_threads,maxit=500L,geogfn=NULL,audit=NULL) {
    if(is.null(audit))audit<-new.env(parent=emptyenv())
    audit$raw_obj_null<-audit$raw_obj_alt<-audit$raw_d_stat<-NA_real_
    audit$alternative_fits<-list()
    nf<-tryCatch(.bgb_s_fit(null,geogfn=geogfn,num_threads=num_threads,maxit=maxit),error=function(e) {
        audit$null_fit<-if(is.null(e$fit_snapshot))list(error=conditionMessage(e)) else e$fit_snapshot
        if(!is.null(e$fit_snapshot))audit$raw_obj_null<-e$fit_snapshot$obj
        stop(e)
    })
    audit$null_fit<-.bgb_s_fit_snapshot(nf);audit$raw_obj_null<-nf$total_loglikelihood
    attempts<-list();errors<-character()
    starts<-list(alternative_start=.bgb_u_params(alt),injected_null_start=.bgb_u_params(nf))
    for(source in names(starts)) {
        start<-starts[[source]]
        value<-tryCatch(.bgb_s_fit(alt,start=start,geogfn=geogfn,num_threads=num_threads,maxit=maxit),error=function(e)e)
        snapshot<-if(inherits(value,'error')) {
            out<-if(is.null(value$fit_snapshot))list() else value$fit_snapshot
            out$error<-conditionMessage(value);out
        } else .bgb_s_fit_snapshot(value)
        audit$alternative_fits[[length(audit$alternative_fits)+1L]]<-c(list(start_source=source,initial_params=start),snapshot)
        if(inherits(value,'error'))errors<-c(errors,conditionMessage(value)) else attempts[[length(attempts)+1L]]<-value
    }
    if(!length(attempts)) {
        raw<-vapply(audit$alternative_fits,function(f)if(is.null(f$obj))NA_real_ else f$obj,numeric(1))
        if(any(!is.na(raw))) {
            audit$raw_obj_alt<-max(raw,na.rm=TRUE)
            audit$raw_d_stat<-2*(audit$raw_obj_alt-audit$raw_obj_null)
        }
        stop('alternative fit failed: ',paste(errors,collapse='; '))
    }
    af<-attempts[[which.max(vapply(attempts,`[[`,numeric(1),'total_loglikelihood'))]]
    audit$raw_obj_alt<-af$total_loglikelihood;audit$raw_d_stat<-2*(af$total_loglikelihood-nf$total_loglikelihood)
    if(2*(af$total_loglikelihood-nf$total_loglikelihood) < -1e-5)stop('negative nested statistic after null-start injection')
    list(null=nf,alt=af)
}

bootstrap_lrt <- function(res_null,res_alt,n_boot=399L,target_conditioning=c('match_fits','tip_survival','unconditional_complete','legacy_unconditioned'),
    refit_observed=FALSE,root_distribution=c('model_prior','empirical_posterior'),calc_surrogate_power=FALSE,n_alt_sim=199L,
    alpha=.05,num_threads=getOption('RustBEARS.threads',1L),seedval=12345L,checkpoint_file=NULL,
    simulation_conditioning='from_model',resume=FALSE,maxit=500L) {
    began<-proc.time()[3];n_boot<-.bgb_s_positive(n_boot,'n_boot');num_threads<-.bgb_s_positive(num_threads,'num_threads')
    maxit<-.bgb_s_positive(maxit,'maxit');target_conditioning<-match.arg(target_conditioning);root_distribution<-match.arg(root_distribution)
    if(length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)stop('alpha must be in (0,1)')
    if(calc_surrogate_power)n_alt_sim<-.bgb_s_positive(n_alt_sim,'n_alt_sim') else n_alt_sim<-0L
    .bgb_s_check(res_null);.bgb_s_check(res_alt)
    pair<-.bgb_u_prepare_comparison(res_null,res_alt,target_conditioning,refit_observed)
    nested<-.bgb_s_nested(pair$null,pair$alt)
    observed<-.bgb_s_pair(pair$null,nested$alt,num_threads,maxit)
    nf<-observed$null;af<-observed$alt
    if(.bgb_u_converged(pair$null)&&pair$null$total_loglikelihood>nf$total_loglikelihood)nf<-pair$null
    if(.bgb_u_converged(nested$alt)&&nested$alt$total_loglikelihood>af$total_loglikelihood)af<-nested$alt
    if(2*(af$total_loglikelihood-nf$total_loglikelihood) < -1e-5)stop('observed alternative remains below null')
    d_obs<-max(0,2*(af$total_loglikelihood-nf$total_loglikelihood));resolved<-.bgb_s_condition(nf,simulation_conditioning)
    tasks<-list()
    for(phase in c('null','alt'))for(rep in seq_len(if(phase=='null')n_boot else n_alt_sim))for(role in c('simulate','fit_null','fit_alt')) {
        id<-paste(phase,rep,role,sep=':');model<-if(role=='fit_alt'||role=='simulate'&&phase=='alt')af else nf
        tasks[[id]]<-list(purpose=paste0('bootstrap_',phase),model_id=.bgb_s_name(model),tree_id='single_tree',outer_rep=rep,role=role)
    }
    rng<-.bgb_s_tasks(seedval,tasks)
    state<-.bgb_s_store(checkpoint_file,list(kind='bootstrap_lrt.R.v2',null=.bgb_s_manifest(nf),alt=.bgb_s_manifest(af),
        n_boot=n_boot,n_alt_sim=n_alt_sim,alpha=alpha,seed=seedval,simulation=resolved,root_distribution=root_distribution,
        maxit=maxit,random_tasks=tasks,recomputed=pair$recomputed_observed_under_new_contract),resume)
    temporary<-tempfile('rustbears-bootstrap-',fileext='.data');on.exit(unlink(temporary),add=TRUE)
    for(phase in c('null','alt')) {
        count<-if(phase=='null')n_boot else n_alt_sim;if(!count)next
        generator<-if(phase=='null')nf else af
        sampler<-prepare_tree_sampler(generator,resolved,root_distribution)
        tryCatch(for(rep in seq_len(count)) {
            id<-paste(phase,rep,sep=':');if(!is.null(state$records[[id]]))next
            keys<-paste(phase,rep,c('simulate','fit_null','fit_alt'),sep=':');task_rng<-rng[keys]
            errors<-character();fits<-NULL;batch<-NULL;audits<-list()
            for(attempt in 1:2) {
                audit<-new.env(parent=emptyenv())
                audit$attempt<-attempt;audit$maxiter_multiplier<-attempt
                audit$initial_null<-.bgb_u_params(nf);audit$initial_alt<-.bgb_u_params(af)
                audit$raw_obj_null<-audit$raw_obj_alt<-audit$raw_d_stat<-NA_real_
                answer<-tryCatch({
                    if(is.null(batch))batch<-simulate_datasets(sampler,n_reps=1L,start=rep,seed=task_rng[[1]]$seed,num_threads=num_threads,return_node_states=TRUE)
                    .bgb_u_write_simulated(batch,1L,temporary)
                    .bgb_s_pair(nf,af,num_threads,maxit*attempt,temporary,audit=audit)
                },error=function(e)e)
                audit$status<-if(inherits(answer,'error'))'failed' else 'converged'
                audit$error<-if(inherits(answer,'error'))conditionMessage(answer) else ''
                audits[[attempt]]<-as.list(audit)
                if(inherits(answer,'error'))errors<-c(errors,conditionMessage(answer)) else {fits<-answer;break}
            }
            row<-list(id=id,rep_id=rep,phase=phase,status=if(is.null(fits))'failed' else 'converged',retried=length(errors)>0,
                obj_null=NA_real_,obj_alt=NA_real_,d_stat=NA_real_,error=paste(errors,collapse=' | '),rng=task_rng,
                simulated_dataset=batch,fits=fits,attempts=audits,
                raw_obj_null=audit$raw_obj_null,raw_obj_alt=audit$raw_obj_alt,raw_d_stat=audit$raw_d_stat)
            if(!is.null(fits)) {
                row$obj_null<-fits$null$total_loglikelihood;row$obj_alt<-fits$alt$total_loglikelihood
                row$d_stat<-max(0,2*(row$obj_alt-row$obj_null));row$params_null<-.bgb_u_params(fits$null);row$params_alt<-.bgb_u_params(fits$alt)
            }
            state$records[[id]]<-row;.bgb_u_save_checkpoint(state)
        },finally=close_prepared_sampler(sampler))
    }
    frame_for<-function(phase,count)do.call(rbind,lapply(seq_len(count),function(rep) {
        row<-state$records[[paste(phase,rep,sep=':')]]
        as.data.frame(row[c('rep_id','phase','status','retried','obj_null','obj_alt','d_stat',
            'raw_obj_null','raw_obj_alt','raw_d_stat','error')],stringsAsFactors=FALSE)
    }))
    frame<-frame_for('null',n_boot);ok<-frame$status=='converged';failed<-sum(!ok);r<-sum(frame$d_stat[ok]>=d_obs)
    bounds<-.bgb_s_failure_intervals(r,failed,n_boot);power<-critical<-alt_frame<-power_bounds<-power_ci<-NULL
    if(n_alt_sim) {
        reference<-c(frame$d_stat[ok],rep(Inf,failed));critical<-sort(reference)[ceiling((1-alpha)*(n_boot-1))+1L]
        alt_frame<-frame_for('alt',n_alt_sim);aok<-alt_frame$status=='converged';hits<-sum(alt_frame$d_stat[aok]>critical)
        power<-hits/n_alt_sim;power_bounds<-c(power,(hits+sum(!aok))/n_alt_sim)
        power_ci<-.bgb_s_failure_intervals(hits,sum(!aok),n_alt_sim)$mc_ci95
    }
    structure(list(null_model=.bgb_s_name(nf),alt_model=.bgb_s_name(af),likelihood_conditioning=nf$likelihood_conditioning,
        recomputed_observed_under_new_contract=pair$recomputed_observed_under_new_contract,root_distribution=root_distribution,
        precision_mode=nf$precision_mode,d_obs=d_obs,df_nominal=length(nested$extra),pval_asymptotic=stats::pchisq(d_obs,length(nested$extra),lower.tail=FALSE),
        n_planned=n_boot,n_converged=sum(ok),n_retried_recovered=sum(ok&frame$retried),n_failed=failed,n_exceedances=r,
        pval_mc_conservative=bounds$failure_bounds[2],pval_failure_bounds=bounds$failure_bounds,pval_mc_ci95=bounds$mc_ci95,
        replicates_df=frame,surrogate_power_fixed_threshold=power,critical_value=critical,surrogate_replicates_df=alt_frame,
        surrogate_power_failure_bounds=power_bounds,surrogate_power_mc_ci95=power_ci,simulation_conditioning=resolved,
        observed_null_fit=nf,observed_alt_fit=af,replicate_artifacts=state$records,checkpoint_fingerprint=state$fingerprint,
        elapsed_sec=unname(proc.time()[3]-began)),class='rustbears_bootstrap_lrt')
}

.bgb_s_grid <- function(mle,bound,n) {
    if(bound==mle)return(mle)
    if(bound==0&&mle>0)g<-c(mle,exp(seq(log(mle*.8),log(max(mle*1e-8,1e-15)),length.out=max(2,n-2L))),0)
    else if(bound>mle&&mle>=0&&bound>100*max(mle,1e-8)) {
        scale<-max(mle,1e-6);g<-mle+expm1(seq(0,log1p((bound-mle)/scale),length.out=n))*scale
    } else g<-seq(mle,bound,length.out=n)
    g[1]<-mle;g[length(g)]<-bound;sort(unique(g),decreasing=bound<mle)
}
.bgb_s_endpoint <- function(point,grid,objective,critical,tol,physical=FALSE) {
    rows<-lapply(grid,point);drop<-vapply(rows,function(row)2*(objective-row$profile_obj)-critical,numeric(1))
    outside<-which(drop>0)
    if(!length(outside))return(list(value=if(physical&&tail(grid,1)==0)0 else NA_real_,
        status=if(physical&&tail(grid,1)==0)'at_physical_zero' else if(tail(grid,1)<grid[1])'censored_at_box_lower' else 'censored_at_box_upper',
        residual=0,disconnected=FALSE))
    i<-outside[1];if(i==1L)stop('profile at fitted maximum is below the reference objective')
    a<-grid[i-1L];b<-grid[i];fa<-drop[i-1L];fb<-drop[i];root<-NA_real_;residual<-Inf
    for(iter in seq_len(250L)) {
        mid<-(a+b)/2;if(mid==a||mid==b)break
        fm<-2*(objective-point(mid)$profile_obj)-critical
        if(is.na(fm))stop('profile endpoint has an invalid objective')
        if(is.finite(fm)&&abs(fm)<tol) {root<-mid;residual<-abs(fm);break}
        if(fm>0){b<-mid;fb<-fm}else{a<-mid;fa<-fm}
    }
    if(!is.finite(root)||residual>=tol)stop('ProfileOptimizationError: endpoint deviance residual exceeds tolerance')
    list(value=root,status=if(physical&&isTRUE(tail(rows,1)[[1]]$structural_zero))'infeasible_zero_approached' else 'interior_bracketed',
        residual=residual,disconnected=any(drop[seq_along(drop)>i]<=0))
}
.bgb_s_hessian <- function(res,num_threads) {
    p<-res$outputs@params_table;names<-rownames(p)[p$type=='free'];theta<-stats::setNames(p[names,'est'],names)
    if(!length(names))return(list(se=NULL,corr=NULL))
    h<-pmax(abs(theta)*1e-3,1e-5)
    if(any(theta-h<=p[names,'min'])||any(theta+h>=p[names,'max']))return(list(se=NULL,corr=NULL))
    score<-.bgb_s_score_factory(res,num_threads);objective<-function(v)-score(v);f0<-objective(theta)
    H<-matrix(0,length(names),length(names),dimnames=list(names,names))
    for(i in seq_along(names)) {
        ei<-theta*0;ei[i]<-h[i];H[i,i]<-(objective(theta+ei)-2*f0+objective(theta-ei))/h[i]^2
        if(i>1L)for(j in seq_len(i-1L)) {
            ej<-theta*0;ej[j]<-h[j]
            H[i,j]<-H[j,i]<-(objective(theta+ei+ej)-objective(theta+ei-ej)-objective(theta-ei+ej)+objective(theta-ei-ej))/(4*h[i]*h[j])
        }
    }
    if(any(!is.finite(H))||min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)<=0)return(list(se=NULL,corr=NULL))
    covariance<-solve(H);list(se=stats::setNames(sqrt(diag(covariance)),names),corr=stats::cov2cor(covariance))
}
.bgb_s_profile_plugin <- function(res,param,n_boot,level,conditioning,rng,state,commit,num_threads,maxit) {
    generation<-digest::digest(list(params=.bgb_u_params(res),objective=res$total_loglikelihood),algo='sha256')
    sampler<-prepare_tree_sampler(res,conditioning,'model_prior');on.exit(close_prepared_sampler(sampler),add=TRUE)
    temporary<-tempfile('rustbears-plugin-',fileext='.data');on.exit(unlink(temporary),add=TRUE)
    rows<-list()
    for(rep in seq_len(n_boot)) {
        id<-paste('plugin',param,generation,rep,sep=':')
        row<-state()$records[[id]]
        if(is.null(row)) {
            keys<-paste(param,rep,c('simulate','full','constrained'),sep=':');task_rng<-rng[keys]
            batch<-simulate_datasets(sampler,n_reps=1L,start=rep,seed=task_rng[[1]]$seed,num_threads=num_threads,return_node_states=TRUE)
            .bgb_u_write_simulated(batch,1L,temporary)
            full<-.bgb_s_fit(res,geogfn=temporary,num_threads=num_threads,maxit=maxit)
            point<-.bgb_s_point(full,stats::setNames(list(.bgb_u_params(res)[[param]]),param),num_threads,maxit)
            if(point$profile_obj>full$total_loglikelihood+1e-6)
                full<-.bgb_s_fit(full,start=point$params,num_threads=num_threads,maxit=2L*maxit)
            if(!is.finite(point$profile_obj)||point$profile_obj>full$total_loglikelihood+5e-6)
                stop('bootstrap-plugin unrestricted fit did not dominate its constraint')
            row<-list(id=id,kind='plugin',param=param,rep_id=rep,generation=generation,
                d_stat=max(0,2*(full$total_loglikelihood-point$profile_obj)),obj_full=full$total_loglikelihood,
                obj_constrained=point$profile_obj,params_full=.bgb_u_params(full),params_constrained=point$params,
                full_fit=full,constrained_fit=point$fit,simulated_dataset=batch,rng=task_rng,
                n_evals=full$n_evals_uncertainty+point$n_evals)
            commit(id,row)
        }
        rows[[rep]]<-row
    }
    frame<-do.call(rbind,lapply(rows,function(row)as.data.frame(row[c('rep_id','d_stat','obj_full','obj_constrained','n_evals')])))
    list(critical=sort(frame$d_stat)[ceiling(level*(n_boot-1))+1L],frame=frame)
}

profile_ci <- function(res,params=NULL,level=.95,method=c('asymptotic_chi2','bootstrap_plugin'),n_boot_plugin=199L,
    include_onesided_zero_ref=FALSE,eval_exact_zero=TRUE,multi_basin_guard=TRUE,n_grid=16L,tol=1e-4,
    num_threads=getOption('RustBEARS.threads',1L),seedval=12345L,simulation_conditioning='from_model',
    checkpoint_file=NULL,resume=FALSE,maxit=500L) {
    began<-proc.time()[3];.bgb_s_check(res);method<-match.arg(method)
    if(length(level)!=1L||!is.finite(level)||level<=0||level>=1)stop('level must lie in (0,1)')
    if(include_onesided_zero_ref&&level<=.5)stop('one-sided zero reference requires level>.5')
    n_grid<-.bgb_s_positive(n_grid,'n_grid');if(n_grid<3L)stop('n_grid must be at least 3')
    num_threads<-.bgb_s_positive(num_threads,'num_threads');maxit<-.bgb_s_positive(maxit,'maxit')
    if(length(tol)!=1L||!is.finite(tol)||tol<=0)stop('tol must be positive and finite')
    if(is.null(params))params<-.bgb_u_free(res)
    if(!length(params)||anyDuplicated(params)||!all(params %in% .bgb_u_free(res)))stop('params must be distinct free parameter names')
    resolved<-NULL;tasks<-list();rng<-list()
    if(method=='bootstrap_plugin') {
        n_boot_plugin<-.bgb_s_positive(n_boot_plugin,'n_boot_plugin');resolved<-.bgb_s_condition(res,simulation_conditioning)
        for(param in params)for(rep in seq_len(n_boot_plugin))for(role in c('simulate','full','constrained'))
            tasks[[paste(param,rep,role,sep=':')]]<-list(purpose='profile_plugin',model_id=.bgb_s_name(res),
                tree_id='single_tree',outer_rep=rep,role=paste0(param,'_',role))
        rng<-.bgb_s_tasks(seedval,tasks)
    }
    state<-.bgb_s_store(checkpoint_file,list(kind='profile_ci.R.v1',fit=.bgb_s_manifest(res),params=params,level=level,
        method=method,n_boot_plugin=n_boot_plugin,include_onesided_zero_ref=include_onesided_zero_ref,eval_exact_zero=eval_exact_zero,
        multi_basin_guard=multi_basin_guard,n_grid=n_grid,tol=tol,seed=seedval,simulation_conditioning=resolved,maxit=maxit,random_tasks=tasks),resume)
    commit<-function(id,row){state$records[[id]]<<-row;.bgb_u_save_checkpoint(state)}
    kind<-if(res$likelihood_conditioning=='legacy_unconditioned')'profile_support_interval' else
        if(method=='bootstrap_plugin')'bootstrap_plugin_interval' else 'asymptotic_confidence_interval'
    original<-res;rebased<-FALSE;intervals<-list();finished<-FALSE
    for(restart in 0:3) {
        intervals<-list();candidate<-NULL
        for(param in params) {
            used<-new.env(parent=emptyenv())
            generation<-digest::digest(list(box=res$outputs@params_table[,c('type','min','max'),drop=FALSE],
                params=.bgb_u_params(res),objective=res$total_loglikelihood),algo='sha256')
            point<-function(value) {
                id<-paste('point',generation,param,sprintf('%a',value),sep=':')
                row<-state$records[[id]]
                if(is.null(row)) {
                    row<-.bgb_s_point(res,stats::setNames(list(value),param),num_threads,maxit)
                    row$param<-param;row$value<-value;row$kind<-'point';row$id<-id;commit(id,row)
                }
                used[[id]]<-row;row
            }
            p<-res$outputs@params_table;mle<-p[param,'est'];physical<-eval_exact_zero&&param %in% c('d','e','a','j')&&p[param,'min']>=0
            left<-.bgb_s_grid(mle,if(physical)0 else p[param,'min'],n_grid);right<-.bgb_s_grid(mle,p[param,'max'],n_grid)
            invisible(lapply(unique(c(left,right)),point))
            best<-function(){values<-as.list(used);values[[which.max(vapply(values,`[[`,numeric(1),'profile_obj'))]]}
            if(best()$profile_obj>res$total_loglikelihood+1e-6){candidate<-best();break}
            critical<-stats::qchisq(level,1);calibration<-NULL
            if(method=='bootstrap_plugin') {
                plugin<-.bgb_s_profile_plugin(res,param,n_boot_plugin,level,resolved,rng,function()state,commit,num_threads,maxit)
                critical<-plugin$critical;calibration<-plugin$frame
            }
            low<-.bgb_s_endpoint(point,left,res$total_loglikelihood,critical,tol,physical)
            high<-.bgb_s_endpoint(point,right,res$total_loglikelihood,critical,tol)
            one<-NA_real_
            if(include_onesided_zero_ref&&physical&&!point(0)$structural_zero)
                one<-.bgb_s_endpoint(point,right,res$total_loglikelihood,stats::qchisq(2*level-1,1),tol)$value
            if(best()$profile_obj>res$total_loglikelihood+1e-6){candidate<-best();break}
            rows<-as.list(used);rows<-rows[order(vapply(rows,`[[`,numeric(1),'value'))]
            frame<-do.call(rbind,lapply(rows,function(row) {
                out<-data.frame(value=row$value,profile_obj=row$profile_obj,delta_deviance=max(0,2*(res$total_loglikelihood-row$profile_obj)),
                    converged=row$converged,structural_zero=row$structural_zero,n_evals=row$n_evals)
                for(k in setdiff(names(row$params),param))out[[paste0('nuisance_',k)]]<-row$params[[k]]
                out
            }))
            intervals[[param]]<-list(param=param,interval_kind=kind,mle=mle,lower=low$value,upper=high$value,
                status_lower=low$status,status_upper=high$status,onesided_upper_from_zero=one,se_hessian_local=NA_real_,
                endpoint_deviance_residual_max=max(low$residual,high$residual),n_evals=sum(frame$n_evals),profile_curve_df=frame,
                critical_value=critical,disconnected_support_detected=low$disconnected||high$disconnected,bootstrap_plugin_df=calibration)
        }
        if(is.null(candidate)){finished<-TRUE;break}
        if(!multi_basin_guard)stop('ProfileOptimizationError: profile exceeds maximum; enable multi_basin_guard')
        p<-res$outputs@params_table
        for(k in .bgb_u_free(res))if(k %in% c('d','e','a','j')&&candidate$params[[k]]==0&&p[k,'min']>0)p[k,'min']<-0
        res$outputs@params_table<-p
        improved<-.bgb_s_fit(res,start=candidate$params,num_threads=num_threads,maxit=2L*maxit)
        if(improved$total_loglikelihood<candidate$profile_obj-1e-6)stop('ProfileOptimizationError: global rebase failed')
        res<-improved;rebased<-TRUE
    }
    if(!finished)stop('ProfileOptimizationError: maximum continued changing after three rebases')
    hessian<-.bgb_s_hessian(res,num_threads)
    for(param in names(intervals))if(param %in% names(hessian$se))intervals[[param]]$se_hessian_local<-hessian$se[[param]]
    structure(list(level=level,method=method,interval_kind=kind,likelihood_conditioning=res$likelihood_conditioning,
        precision_mode=res$precision_mode,obj_mle=res$total_loglikelihood,rebased_mle=rebased,updated_fit=if(rebased)res else NULL,
        intervals=intervals,hessian_corr=hessian$corr,total_evals=sum(vapply(state$records,function(r)r$n_evals,numeric(1))),
        elapsed_sec=unname(proc.time()[3]-began),checkpoint_fingerprint=state$fingerprint,
        interpretation='Connected profile component containing the fitted maximum; bootstrap_plugin is local calibration at that estimate.'),
        class='rustbears_profile_ci')
}

profile_surface <- function(res,param_x='d',param_y='j',grid_size=25L,reoptimize_nuisance=TRUE,
    num_threads=getOption('RustBEARS.threads',1L),grid_x=NULL,grid_y=NULL,maxit=500L) {
    .bgb_s_check(res)
    if(param_x==param_y||!all(c(param_x,param_y) %in% .bgb_u_free(res)))stop('axes must be distinct free parameters')
    grid_size<-.bgb_s_positive(grid_size,'grid_size');num_threads<-.bgb_s_positive(num_threads,'num_threads')
    grid<-function(values,name) {
        p<-res$outputs@params_table
        if(is.null(values))return(seq(p[name,'min'],p[name,'max'],length.out=grid_size))
        if(!is.numeric(values)||!length(values)||any(!is.finite(values))||any(values<p[name,'min'])||any(values>p[name,'max']))stop('surface grid is outside the parameter box')
        values
    }
    rows<-list()
    for(x in grid(grid_x,param_x))for(y in grid(grid_y,param_y)) {
        fixed<-stats::setNames(list(x,y),c(param_x,param_y))
        if(!reoptimize_nuisance)fixed<-c(fixed,as.list(.bgb_u_params(res)[setdiff(.bgb_u_free(res),c(param_x,param_y))]))
        point<-.bgb_s_point(res,fixed,num_threads,maxit)
        row<-data.frame(x=x,y=y,profile_obj=point$profile_obj,delta_deviance=2*(res$total_loglikelihood-point$profile_obj),
            converged=point$converged,structural_zero=point$structural_zero)
        names(row)[1:2]<-c(param_x,param_y)
        for(k in setdiff(names(point$params),c(param_x,param_y)))row[[paste0('nuisance_',k)]]<-point$params[[k]]
        rows[[length(rows)+1L]]<-row
    }
    out<-do.call(rbind,rows);class(out)<-c('rustbears_profile_surface','data.frame');out
}

as.data.frame.rustbears_profile_ci <- function(x,...)do.call(rbind,lapply(x$intervals,function(i)
    as.data.frame(i[c('param','interval_kind','mle','lower','upper','status_lower','status_upper',
        'onesided_upper_from_zero','se_hessian_local','endpoint_deviance_residual_max','critical_value','disconnected_support_detected','n_evals')],stringsAsFactors=FALSE)))
print.rustbears_profile_ci <- function(x,...) {print(as.data.frame(x),...);invisible(x)}
plot.rustbears_profile_ci <- function(x,...) {
    old<-graphics::par(mfrow=c(ceiling(length(x$intervals)/2),min(2,length(x$intervals))));on.exit(graphics::par(old),add=TRUE)
    for(i in x$intervals) {graphics::plot(i$profile_curve_df$value,i$profile_curve_df$delta_deviance,type='l',xlab=i$param,ylab='Twice profile score drop',...)
        graphics::abline(h=i$critical_value,lty=2,col='red')}
    invisible(x)
}
print.rustbears_bootstrap_lrt <- function(x,...) {
    print(data.frame(d_obs=x$d_obs,n_planned=x$n_planned,n_failed=x$n_failed,pval_mc_conservative=x$pval_mc_conservative,
        failure_lower=x$pval_failure_bounds[1],failure_upper=x$pval_failure_bounds[2],mc_lower=x$pval_mc_ci95[1],mc_upper=x$pval_mc_ci95[2]),...)
    invisible(x)
}
plot.rustbears_bootstrap_lrt <- function(x,...) {
    values<-x$replicates_df$d_stat[x$replicates_df$status=='converged']
    if(!length(values))stop('no successful bootstrap statistics to plot')
    graphics::hist(values,xlab='Bootstrap score difference',main='Fixed fitted-model bootstrap',...)
    graphics::abline(v=x$d_obs,col='red');invisible(x)
}
