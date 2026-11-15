# Scientific assertions for inherited examples. Historical optimizer coordinates
# and one RNG realization are not reference probability laws. See
# docs/LEGACY_TEST_REMEDIATION.md in the source repository for retained failures.
.bgb_test_assert_three_epoch_maps <- function(counts, ana, clado, areas) {
    testthat::expect_length(ana, 50L)
    testthat::expect_length(clado, 50L)
    kinds <- c('founder','a','d','e','subset','vicariance','sympatry')
    per_map <- matrix(0, length(ana), length(kinds), dimnames=list(NULL,kinds))
    directions <- matrix(0, length(areas), length(areas), dimnames=list(areas,areas))
    for (i in seq_along(ana)) {
        a <- ana[[i]]; c <- clado[[i]]
        testthat::expect_true(all(is.finite(a$abs_event_time)))
        testthat::expect_true(all(a$abs_event_time >= a$time_top-1e-12 & a$abs_event_time <= a$time_bot+1e-12))
        testthat::expect_true(all(a$event_type %in% c('d','e')))
        for (column in c('current_rangetxt','new_rangetxt')) {
            middle <- a$abs_event_time > .1+1e-12 & a$abs_event_time < 1.5-1e-12
            older <- a$abs_event_time > 1.5+1e-12
            testthat::expect_false(any(grepl('C',a[[column]][middle])))
            testthat::expect_true(all(a[[column]][older] == 'A'))
        }
        current <- strsplit(a$current_rangetxt,'',fixed=TRUE)
        next_range <- strsplit(a$new_rangetxt,'',fixed=TRUE)
        for (j in seq_len(nrow(a))) {
            if (a$event_type[j] == 'd') {
                testthat::expect_length(setdiff(next_range[[j]],current[[j]]),1L)
                testthat::expect_length(setdiff(current[[j]],next_range[[j]]),0L)
            } else {
                testthat::expect_length(setdiff(current[[j]],next_range[[j]]),1L)
                testthat::expect_length(setdiff(next_range[[j]],current[[j]]),0L)
            }
        }
        ct <- sub('.*\\(([ysvj])\\).*','\\1',c$clado_event_type)
        ct <- ct[!is.na(ct) & ct %in% c('y','s','v','j')]
        testthat::expect_length(ct,2L)
        per_map[i,] <- c(sum(ct=='j'),sum(a$event_type=='a'),sum(a$event_type=='d'),
                        sum(a$event_type=='e'),sum(ct=='s'),sum(ct=='v'),sum(ct=='y'))
        rows <- which(a$event_type %in% c('d','a'))
        for (j in rows) {
            from <- a$ana_dispersal_from[j]; to <- a$dispersal_to[j]
            testthat::expect_true(from %in% areas && to %in% areas && from != to)
            directions[from,to] <- directions[from,to]+1
        }
        rows <- which(grepl('\\(j\\)',c$clado_event_type))
        for (j in rows) {
            from <- c$clado_dispersal_from[j]; to <- c$clado_dispersal_to[j]
            testthat::expect_true(from %in% areas && to %in% areas && from != to)
            directions[from,to] <- directions[from,to]+1
        }
    }
    per_map <- cbind(per_map,ALL_disp=per_map[,'founder']+per_map[,'a']+per_map[,'d'],
        ana_disp=per_map[,'a']+per_map[,'d'],all_ana=rowSums(per_map[,c('a','d','e'),drop=FALSE]),
        all_clado=rowSums(per_map[,c('founder','subset','vicariance','sympatry'),drop=FALSE]),
        total_events=rowSums(per_map))
    reference <- rbind(means=colMeans(per_map),stdevs=apply(per_map,2,sd),sums=colSums(per_map))
    testthat::expect_equal(unname(as.matrix(counts$summary_counts_BSMs[rownames(reference),colnames(reference)])),
                          unname(reference),tolerance=1e-12)
    testthat::expect_equal(unname(as.matrix(counts$all_dispersals_counts_fromto_means)),
                          unname(directions/length(ana)),tolerance=1e-12)
}

.bgb_test_assert_three_epoch_fixed_score <- function(inputs) {
    previous <- options(RustBEARS.precision='LegacyF32')
    on.exit(options(previous),add=TRUE)
    x <- inputs
    for (p in c('d','e','j')) x$BioGeoBEARS_model_object@params_table[p,c('init','est')] <- c(d=.4,e=.2,j=.1)[[p]]
    x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object,update_init=TRUE)
    x$calc_ancprobs <- FALSE; x$print_optim <- FALSE
    observed <- bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_loglike')
    # Independent BioGeoBEARS 1.1.4 @ 876a48c, genuine cladoRcpp/rexpokit,
    # OLD (.1,1.5,100) epochs, noC2 masks, eight full states, raw root measure.
    testthat::expect_lt(abs(observed-(-18.410726319261713)),1e-10)
}

# Full fixed-parameter evaluation preserves trait observations and linked rows.
.bgb_test_fixed <- function(res, candidate) {
    x <- res$inputs
    x$BioGeoBEARS_model_object@params_table <- res$outputs@params_table
    p <- names(candidate)
    x$BioGeoBEARS_model_object@params_table[p,'init'] <- candidate
    x$BioGeoBEARS_model_object@params_table[p,'est'] <- candidate
    # A fixed-point oracle may include an archived value outside today's
    # optimizer box. Preserve the value; only this non-optimizing container's
    # box is expanded so the legacy input checker permits its evaluation.
    x$BioGeoBEARS_model_object@params_table[p,'min'] <- pmin(x$BioGeoBEARS_model_object@params_table[p,'min'],candidate)
    x$BioGeoBEARS_model_object@params_table[p,'max'] <- pmax(x$BioGeoBEARS_model_object@params_table[p,'max'],candidate)
    x$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(x$BioGeoBEARS_model_object,update_init=TRUE)
    x$use_optimx <- FALSE; x$calc_ancprobs <- TRUE
    x$print_optim <- FALSE; x$printlevel <- -1L; x$num_cores_to_use <- 1L
    bears_optim_run(x,skip_optim=TRUE,skip_optim_option='return_all')
}

.bgb_test_plain <- function(x) {
    if (isS4(x) || is.environment(x) || is.function(x) || typeof(x) %in% c('externalptr','weakref','language')) return(FALSE)
    if (is.list(x)) return(all(vapply(x,.bgb_test_plain,logical(1))))
    is.null(x) || is.atomic(x)
}

.bgb_test_oracle <- function(res,candidates,label) {
    stopifnot(length(find.package('BioGeoBEARS',lib.loc=c(getOption('RustBEARS.test_oracle_library'),.libPaths()),quiet=TRUE))==1L)
    allowed <- setdiff(names(formals(define_BioGeoBEARS_run)),c('BioGeoBEARS_model_object','states_list','tmpwd'))
    allowed <- c(allowed,'trait_Pmat_txt')
    config <- res$inputs[intersect(allowed,names(res$inputs))]
    stopifnot(all(vapply(config,.bgb_test_plain,logical(1))))
    # Serialize only data. In particular, never load a RustBEARS S4 class in the
    # subprocess: doing so would silently turn this into a self-comparison.
    payload <- list(run_config=config,params_table=res$outputs@params_table,
                    candidates=candidates,oracle_library=getOption('RustBEARS.test_oracle_library'))
    if (!is.null(res$inputs$trait)) payload$trait_df <- res$inputs$trait@df
    folder <- file.path(getOption('RustBEARS.test_evidence_dir',getwd()),'legacy-oracle-evidence',label)
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    folder <- normalizePath(folder,mustWork=TRUE)
    input <- file.path(folder,'input.rds'); output <- file.path(folder,'output.rds')
    saveRDS(payload,input,version=2)
    runner <- getOption('RustBEARS.test_oracle_runner',system.file('validation','legacy_trait_oracle.R',package='RustBEARS'))
    stopifnot(file.exists(runner))
    # Pass the active test library explicitly to a clean process on every OS.
    old <- Sys.getenv('R_LIBS',unset=NA_character_)
    on.exit(if(is.na(old))Sys.unsetenv('R_LIBS') else Sys.setenv(R_LIBS=old),add=TRUE)
    Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
    status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(runner),shQuote(input),shQuote(output)),
                      stdout=file.path(folder,'stdout.log'),stderr=file.path(folder,'stderr.log'))
    testthat::expect_identical(status,0L,info=paste('independent oracle:',label))
    if(status!=0L)stop('upstream oracle failed; inspect ',folder)
    answer <- readRDS(output)
    testthat::expect_false(answer$metadata$loaded_rustbears)
    answer
}

.bgb_test_assert_archive_fit <- function(res,candidates,label) {
    p <- res$outputs@params_table; free <- rownames(p)[p$type=='free']
    current <- setNames(p[free,'est'],free)
    testthat::expect_true(all(is.finite(current)))
    testthat::expect_true(all(current>=p[free,'min'] & current<=p[free,'max']))
    testthat::expect_true(is.finite(res$total_loglikelihood))
    fixed <- .bgb_test_fixed(res,current)
    testthat::expect_lt(abs(fixed$total_loglikelihood-res$total_loglikelihood),1e-8)
    o <- res$optim_result
    if(!is.null(o$convergence)) {
        testthat::expect_identical(as.integer(o$convergence),0L)
    } else if(!is.null(o$convcode)) {
        testthat::expect_equal(o$convcode[which.max(o$value)],0)
    } else {
        # GenSA has no convergence code and cannot certify a global optimum.
        # Check that an actual search completed; oracle baselines below test its
        # attained objective, not an unsupported global-convergence claim.
        testthat::expect_true(is.matrix(o$trace.mat) && nrow(o$trace.mat)>0L)
        testthat::expect_true(is.numeric(o$counts) && all(is.finite(o$counts)) && all(o$counts>0))
        testthat::expect_lt(abs(-o$value-res$total_loglikelihood),1e-8)
    }
    candidates <- lapply(candidates,function(x) {
        stopifnot(setequal(names(x),free),length(x)==length(free))
        x[free]
    })
    initial <- setNames(res$inputs$BioGeoBEARS_model_object@params_table[free,'init'],free)
    initial_fixed <- .bgb_test_fixed(res,initial)
    testthat::expect_gte(res$total_loglikelihood+1e-8,initial_fixed$total_loglikelihood)
    candidates <- c(list(fitted=current,initial=initial),candidates)
    oracle <- .bgb_test_oracle(res,candidates,label)
    previous <- options(RustBEARS.precision='LegacyF32')
    on.exit(options(previous),add=TRUE)
    native <- lapply(candidates,function(x).bgb_test_fixed(res,x))
    top <- 'relative_probs_of_each_state_at_branch_top_AT_node_UPPASS'
    bottom <- 'relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS'
    for(name in names(candidates)) {
        testthat::expect_true(abs(native[[name]]$total_loglikelihood-oracle$scores[[name]]$lnL)<1e-8,info=paste(label,name,'lnL'))
        for(item in list(c(top,'ancestral_top'),c(bottom,'ancestral_bottom'))) {
            a <- native[[name]][[item[[1]]]]; b <- oracle$scores[[name]][[item[[2]]]]
            testthat::expect_identical(dim(a),dim(b))
            testthat::expect_identical(is.na(a),is.na(b))
            testthat::expect_true(all(is.finite(a[!is.na(a)])))
            testthat::expect_true(all(is.finite(b[!is.na(b)])))
            testthat::expect_true(max(abs(a-b),na.rm=TRUE)<1e-7,info=paste(label,name,item[[2]]))
        }
    }
    eligible <- vapply(oracle$scores,function(x)isTRUE(x$candidate_within_declared_bounds),logical(1))
    testthat::expect_true(eligible[['fitted']])
    # Out-of-box archived candidates remain fixed-score oracle checks, but are
    # not feasible lower bounds for optimization in this declared box.
    baseline <- max(vapply(oracle$scores[eligible],function(x)x$lnL,numeric(1)))
    evidence_path <- file.path(getOption('RustBEARS.test_evidence_dir',getwd()),'legacy-oracle-evidence',label,'comparison.rds')
    saveRDS(list(fit=res,oracle=oracle,native=native,eligible=eligible,baseline=baseline,
                original_search_minus_baseline=native$fitted$total_loglikelihood-baseline,
                separate_bounded_status='not_yet_completed'),evidence_path,version=2)
    # Keep the heuristic fit and its shortfall separate. A random search does
    # not promise the same optimum on every platform. A second, explicitly
    # defined optimization regression starts at each feasible fitted/archive
    # point and selects only an actually converged result.
    starts <- candidates[eligible & names(candidates)!='initial']
    x <- res$inputs
    x$BioGeoBEARS_model_object@params_table <- p
    x$use_optimx <- FALSE; x$print_optim <- FALSE; x$printlevel <- -1L; x$num_cores_to_use <- 1L
    x$calc_ancprobs <- FALSE
    options(RustBEARS.precision='StrictF64')
    objective <- function(value) {
        current_run <- x
        current_run$BioGeoBEARS_model_object@params_table[free,'init'] <- value
        current_run$BioGeoBEARS_model_object@params_table[free,'est'] <- value
        current_run$BioGeoBEARS_model_object <- calc_linked_params_BioGeoBEARS_model_object(current_run$BioGeoBEARS_model_object,update_init=TRUE)
        invisible(capture.output(score <- bears_optim_run(current_run,skip_optim=TRUE,skip_optim_option='return_loglike')))
        -score
    }
    strict_start_scores <- vapply(starts,function(start)-objective(start),numeric(1))
    strict_baseline <- max(strict_start_scores)
    attempts <- lapply(starts,function(start) {
        scale <- pmax(abs(start),1e-2)
        lower <- p[free,'min']/scale; upper <- p[free,'max']/scale
        rhobeg <- min(.1,.2*min(upper-lower))
        # Multiplication can place an exact scaled endpoint one ulp outside
        # its physical box. Map optimizer coordinates back into that same box;
        # the historical candidate inputs themselves are never altered.
        physical <- function(z)pmax(p[free,'min'],pmin(p[free,'max'],z*scale))
        answer <- minqa::bobyqa(par=start/scale,fn=function(z)objective(physical(z)),
            lower=lower,upper=upper,control=list(rhobeg=rhobeg,rhoend=rhobeg*1e-6,maxfun=5000,iprint=0))
        list(par=setNames(physical(answer$par),free),objective=answer$fval,convergence=answer$ierr,
             message=answer$msg,evaluations=answer$feval,raw=answer,scale=scale)
    })
    saveRDS(attempts,file.path(getOption('RustBEARS.test_evidence_dir',getwd()),'legacy-oracle-evidence',label,'bounded_attempts.rds'))
    converged <- vapply(attempts,function(a)a$convergence==0L && is.finite(a$objective),logical(1))
    testthat::expect_true(any(converged),info=paste(label,'independent bounded multi-start convergence'))
    if(!any(converged))stop('no converged independent bounded search: ',label)
    winner <- names(which.min(vapply(attempts[converged],function(a)a$objective,numeric(1))))
    pp <- attempts[[winner]]$par
    testthat::expect_true(all(is.finite(pp) & pp>=p[free,'min'] & pp<=p[free,'max']))
    polished_strict <- .bgb_test_fixed(res,setNames(pp,free))
    testthat::expect_lt(abs(polished_strict$total_loglikelihood+attempts[[winner]]$objective),1e-8)
    options(RustBEARS.precision='LegacyF32')
    polished <- .bgb_test_fixed(res,setNames(pp,free))
    testthat::expect_true(polished_strict$total_loglikelihood+1e-7>=strict_baseline,
                         info=paste(label,'separate bounded optimizer vs StrictF64 feasible baseline'))
    final_oracle <- .bgb_test_oracle(res,list(fitted=setNames(pp,free)),paste0(label,'-bounded'))
    testthat::expect_lt(abs(polished$total_loglikelihood-final_oracle$scores$fitted$lnL),1e-8)
    for(item in list(c(top,'ancestral_top'),c(bottom,'ancestral_bottom'))) {
        a <- polished[[item[[1]]]]; b <- final_oracle$scores$fitted[[item[[2]]]]
        testthat::expect_identical(dim(a),dim(b));testthat::expect_identical(is.na(a),is.na(b))
        testthat::expect_true(all(is.finite(a[!is.na(a)])) && all(is.finite(b[!is.na(b)])))
        testthat::expect_lt(max(abs(a-b),na.rm=TRUE),1e-7)
    }
    deficit <- native$fitted$total_loglikelihood-baseline
    cat('\nLegacy search audit: ',label,', original LegacyF32-minus-archive=',format(deficit,digits=16),
        ', original StrictF64-minus-archive=',format(strict_start_scores[['fitted']]-strict_baseline,digits=16),
        ', separate bounded StrictF64=',format(polished_strict$total_loglikelihood,digits=16),
        ', separate bounded LegacyF32=',format(polished$total_loglikelihood,digits=16),'\n',sep='')
    saveRDS(list(fit=res,oracle=oracle,native=native,eligible=eligible,baseline=baseline,
                original_search_minus_baseline=deficit,original_strict_search_minus_baseline=strict_start_scores[['fitted']]-strict_baseline,strict_start_scores=strict_start_scores,strict_baseline=strict_baseline,final_oracle=final_oracle,separate_bounded_attempts=attempts,
                separate_bounded_winner=winner,separate_bounded_fit_strict=polished_strict,separate_bounded_fit_legacy=polished),
            evidence_path,version=2)
    invisible(res)
}

.bgb_test_trait_archives <- function(fits,prefixes,env,label) {
    for(i in seq_along(fits)) {
        p <- fits[[i]]$outputs@params_table
        free <- rownames(p)[p$type=='free']
        candidates <- lapply(prefixes[[i]],function(prefix) {
            setNames(vapply(free,function(name)get(paste0(prefix,'_',name),envir=env,inherits=FALSE),numeric(1)),free)
        })
        names(candidates) <- prefixes[[i]]
        .bgb_test_assert_archive_fit(fits[[i]],candidates,paste0(label,'-',names(fits)[i]))
    }
    invisible(fits)
}
