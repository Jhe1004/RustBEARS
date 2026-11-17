# GPL >= 2. Prepare immutable R-side topology and call the native evaluator.
.bgb_m3_runtime <- new.env(parent = emptyenv())
.bgb_m3_runtime$contexts <- list()
.bgb_m3_runtime$native_calls <- 0L
.bgb_m3_runtime$fallbacks <- 0L
.bgb_m3_runtime$last_fallback <- NULL

.bgb_m3_clear_contexts <- function() {
    for (entry in .bgb_m3_runtime$contexts)
        .Call(.bgb_c_bgb_tree_release, entry$handle)
    .bgb_m3_runtime$contexts <- list()
    invisible(NULL)
}

.bgb_m3_layer_mask <- function(inputs, states, i, include_null) {
    allowed <- adjacency <- manual <- rep(TRUE, length(states))
    if (is.null(inputs$lists_of_states_lists_0based) &&
        (!is.null(inputs$list_of_areas_allowed_mats) || !is.null(inputs$list_of_areas_adjacency_mats)))
        stop("missing states list for stratified restrictions")
    if (!is.null(inputs$list_of_areas_allowed_mats))
        allowed <- vapply(states, check_if_state_is_allowed, logical(1), inputs$list_of_areas_allowed_mats[[i]])
    if (!is.null(inputs$list_of_areas_adjacency_mats))
        adjacency <- vapply(states, check_if_state_is_allowed_by_adjacency, logical(1), inputs$list_of_areas_adjacency_mats[[i]])
    if (!is.null(inputs$lists_of_states_lists_0based))
        manual <- states %in% inputs$lists_of_states_lists_0based[[i]]
    if (include_null) allowed[1] <- adjacency[1] <- manual[1] <- TRUE
    allowed & adjacency & manual
}

.bgb_m3_fixes <- function(master, rows, fixnode, fixlikes, mask, trait_mask) {
    empty <- list(nodes = numeric(), down = numeric(), up = numeric())
    if (!length(fixnode)) return(empty)
    if (length(fixnode) == 1L) {
        matches <- which(master$node == fixnode & master$SUBnode.type %in% c("root", "internal"))
        stopifnot(length(matches) == 1L)
        if (!matches %in% rows) return(empty)
        return(list(nodes = master$SUBnode[matches] - 1L, down = as.numeric(fixlikes)[mask],
                    up = as.numeric(fixlikes)[trait_mask]))
    }
    requested <- fixnode %in% master$node[rows]
    if (!any(requested)) return(empty)
    nodes <- master$SUBnode[rows][master$node[rows] %in% fixnode]
    stopifnot(length(nodes) == sum(requested))
    order <- order(nodes)
    # Both original loops sort local node numbers after selecting likes in
    # requested master-node order. Preserve that operation order.
    likes <- fixlikes[requested, , drop = FALSE][order, , drop = FALSE]
    list(nodes = nodes[order] - 1L, down = as.numeric(likes[, trait_mask, drop = FALSE]),
         up = as.numeric(likes[, trait_mask, drop = FALSE]))
}

.bgb_m3_branch_spec <- function(master, row, length, minimum, exact_ages = FALSE) {
    if (exact_ages) {
        ok <- is.finite(length) && length > 0 && master$edge.length[row] >= minimum && !isTRUE(master$fossils[row])
        return(c(length, ok, length, ok, length))
    }
    age <- master$time_bp[row]
    top <- master$time_top[row]
    bottom <- master$time_bot[row]
    fossil <- isTRUE(master$fossils[row])
    down <- !(fossil && age > top)
    up <- age >= top
    duration <- length
    if (fossil && age >= top && age < bottom) {
        duration <- length - (age - top)
        down <- up <- TRUE
    }
    if (master$edge.length[row] < minimum) down <- up <- FALSE
    c(length, down, duration, up, duration)
}

.bgb_m3_compile <- function(inputs, phy, areas, states, n, minimum, include_null, fixnode, fixlikes,
                            exact_ages = FALSE, root_prior = numeric()) {
    sections <- inputs$tree_sections_list
    master <- inputs$master_table
    k <- length(inputs$timeperiods)
    nt <- length(phy$tip.label)
    nn <- nt + phy$Nnode
    na <- length(areas)
    nr <- nrow(master)
    trait <- !is.null(inputs$trait)
    ntraits <- if (trait) ncol(inputs$trait_Pmat_txt) else 1L
    stopifnot(k >= 2L, length(sections) == k, is.data.frame(master), nr >= nt,
              n == length(states) * ntraits, is.finite(minimum), na > 0)
    exact_node_ages <- parent_of <- NULL
    if (exact_ages) {
        depths <- ape::node.depth.edgelength(phy)
        exact_node_ages <- max(depths[seq_len(nt)]) - depths
        parent_of <- integer(nn)
        parent_of[phy$edge[, 2]] <- phy$edge[, 1]
    }
    groups <- split(seq_len(nr), paste(master$stratum, master$piecenum, master$piececlass))
    rows_by_node <- split(seq_len(nr), paste(master$stratum, master$node))
    oldest_route <- function(i, node) {
        if (i == k) return(numeric())
        row <- rows_by_node[[paste(i + 1L, node)]]
        stopifnot(length(row) == 1L)
        tip <- if (master$piececlass[row] == "subbranch") 0L else master$SUBnode[row] - 1L
        c(master$piecenum[row] - 1L, tip)
    }
    root <- master$node[which(master$node.type == "root")]
    stopifnot(length(root) == 1L)
    ordered <- function(rows) rows[order(master$node[rows])]
    ttl <- ordered(which(master$SUBnode.type %in% c("internal", "orig_tip", "root")))
    down <- ordered(which(master$piececlass == "subtree" & master$SUBnode.type %in% c("internal", "tip")))
    tip_routes <- vapply(seq_len(nt), function(node) {
        orig <- which(master$node == node & master$piececlass == "orig_tip")
        stopifnot(length(orig) == 1L)
        row <- which(master$node == node & master$piececlass != "orig_tip" & master$time_top == master$time_top[orig])
        stopifnot(length(row) == 1L)
        c(master$stratum[row] - 1L, master$piecenum[row] - 1L,
          if (master$piececlass[row] == "subbranch") 0L else master$SUBnode[row] - 1L)
    }, numeric(3))
    config <- list(c(n, nt, nn, nr, root - 1L, k, include_null), ttl - 1L, down - 1L, as.numeric(tip_routes))
    current_labels <- phy$tip.label
    layers <- vector("list", k)
    ones <- matrix(1, na, na)
    geography <- function(field, index, default) {
        value <- inputs[[field]]
        if (is.null(value)) default else value[[index]]
    }
    empty <- numeric()
    for (i in seq_len(k)) {
        t_top <- if (i == 1L) 0 else inputs$timeperiods[i - 1L]
        t_bot <- inputs$timeperiods[i]
        section <- sections[[i]]
        pieces <- section$return_pieces_list
        mask <- .bgb_m3_layer_mask(inputs, states, i, include_null)
        trait_mask <- rep(mask, ntraits)
        local_states <- states[mask]
        local_n <- sum(trait_mask)
        stopifnot(local_n > as.integer(include_null))
        clado_states <- if (include_null) local_states[-1L] else local_states
        q_data <- if (trait) list(empty, empty) else .bgb_state_buffers(local_states)
        clado_data <- if (trait) list(empty, empty) else .bgb_state_buffers(clado_states)
        config <- c(config, list(c(length(pieces), na), which(trait_mask) - 1L, mask), q_data, clado_data)
        for (j in seq_along(pieces)) {
            piece <- pieces[[j]]
            if (is.numeric(piece)) {
                row <- groups[[paste(i, j, "subbranch")]]
                slot <- which(current_labels == section$return_pieces_basenames[[j]])
                stopifnot(length(row) == 1L, length(slot) == 1L, length(piece) == 1L)
                seg_len <- if (exact_ages) {
                    mnode <- master$node[row]
                    min(t_bot, exact_node_ages[parent_of[mnode]]) - max(t_top, exact_node_ages[mnode])
                } else piece
                spec <- .bgb_m3_branch_spec(master, row, seg_len, minimum, exact_ages = exact_ages)
                config <- c(config, list(c(0, slot - 1L, row - 1L, spec),
                                         oldest_route(i, master$node[row])))
            } else {
                stopifnot(inherits(piece, "phylo"))
                piece <- reorder(piece, "pruningwise")
                piece_tips <- length(piece$tip.label)
                piece_nodes <- piece_tips + piece$Nnode
                rows <- groups[[paste(i, j, "subtree")]]
                rows <- rows[match(seq_len(piece_nodes), master$SUBnode[rows])]
                slots <- which(current_labels %in% piece$tip.label)
                stopifnot(length(slots) == piece_tips, !anyNA(rows),
                          all(master$node.type[rows] %in% c("internal", "tip", "root")))
                piece_root <- piece$edge[nrow(piece$edge), 1]
                master_piece_root <- master$node[rows[piece_root]]
                if (exact_ages) {
                    edge_lengths <- vapply(piece$edge[, 2], function(lc) {
                        mc <- master$node[rows[lc]]
                        exact_node_ages[parent_of[mc]] - max(t_top, exact_node_ages[mc])
                    }, numeric(1))
                    root_length <- if (master_piece_root == root) 0 else t_bot - exact_node_ages[master_piece_root]
                } else {
                    edge_lengths <- piece$edge.length
                    root_length <- if (is.null(piece$root.edge)) 0 else piece$root.edge
                }
                older <- if (root_length > 0) oldest_route(i, master_piece_root) else empty
                fixes <- .bgb_m3_fixes(master, rows, fixnode, fixlikes, mask, trait_mask)
                piece_prior <- if (i == k && master_piece_root == root && length(root_prior))
                    root_prior[which(trait_mask)] else empty
                subtree_config <- list(c(piece_tips, local_n, piece_nodes, minimum, include_null,
                                         TRUE, root_length, 1, 3, 0),
                                       piece$edge - 1L, edge_lengths, empty, fixes$nodes, fixes$down,
                                       empty, piece_prior, empty, empty, empty, empty, c(0, 0))
                config <- c(config, list(1, slots - 1L, rows - 1L, master$node[rows] - 1L,
                                         master$SUBnode.type[rows] == "root", older,
                                         fixes$nodes, fixes$up), subtree_config)
            }
        }
        current_labels <- section$tree_to_chainsaw$tip.label
        num_mat <- function(m, diag_one = FALSE) {
            out <- as.matrix(m); storage.mode(out) <- "double"
            if (diag_one) diag(out)[is.na(diag(out)) | diag(out) == 0] <- 1
            out
        }
        layers[[i]] <- list(states = local_states, width = local_n,
                            distances = num_mat(geography("list_of_distances_mats", i, ones), TRUE),
                            env_down = num_mat(geography("list_of_envdistances_mats", i, ones), TRUE),
                            env_up = num_mat(geography("list_of_envdistances_mats", 1L, ones), TRUE),
                            manual = num_mat(geography("list_of_dispersal_multipliers_mats", i, ones)),
                            area = as.numeric(unlist(geography("list_of_area_of_areas", i, rep(1, na)), use.names = FALSE)),
                            prior_width = max(lengths(clado_states)))
    }
    # The original dense Uppass allocates matrices with for (m in ...).
    # That binding survives to the next layer and scales its jump distances.
    previous_m <- NULL
    for (i in rev(seq_len(k))) {
        layers[[i]]$upper_loop_m <- previous_m
        for (piece in sections[[i]]$return_pieces_list) {
            if (inherits(piece, "phylo")) previous_m <- length(piece$edge.length)
        }
    }
    tip_lookup <- .bgb_m3_compile_tip_lookup(master, groups, rows_by_node)
    list(config = config, layers = layers, nr = nr, nn = nn, nt = nt,
         trait = trait, master = master, tip_lookup = tip_lookup)
}

.bgb_m3_compile_tip_lookup <- function(master, groups, rows_by_node) {
    lookup <- list()
    for (row in which(master$piececlass == "subbranch")) {
        i <- master$stratum[row]
        j <- master$piecenum[row]
        node <- master$node[row]
        key <- paste(i, j)
        if (length(groups[[paste(i, j, "subbranch")]]) != 1L) next
        fossil_here <- master$time_bp[row] >= master$time_top[row] &&
            master$time_bp[row] < master$time_bot[row] && master$fossils[row]
        if (is.na(fossil_here)) next
        if (fossil_here) {
            target <- which(master$stratum == i & master$piecenum == j &
                            master$SUBnode.type == "tip" & master$node == node)
        } else if (i > 1L) {
            target <- rows_by_node[[paste(i - 1L, node)]]
            target <- target[master$piececlass[target] != "orig_tip"]
        } else {
            target <- rows_by_node[[paste(0, node)]]
            target <- target[master$piecenum[target] == j & master$piececlass[target] == "orig_tip"]
        }
        if (length(target) == 1L)
            lookup[[key]] <- list(row = target, down_table = fossil_here || i > 1L)
    }
    lookup
}

.bgb_m3_cached_subbranch_tip <- function(stratum, piecenum, master_table, condlikes, down_table) {
    if (!is.numeric(stratum) || length(stratum) != 1L || is.na(stratum) ||
        !is.numeric(piecenum) || length(piecenum) != 1L || is.na(piecenum))
        return(NULL)
    for (entry in .bgb_m3_runtime$contexts) {
        if (!identical(entry$master, master_table)) next
        route <- entry$tip_lookup[[paste(stratum, piecenum)]]
        if (is.null(route)) return(NULL)
        table <- if (route$down_table) down_table else condlikes
        if (!is.matrix(table) || nrow(table) < route$row) return(NULL)
        values <- table[route$row, ]
        total <- sum(values)
        if (!is.na(total) && total != 0) return(values)
        return(NULL)
    }
    NULL
}

.bgb_m3_jts <- function(text, parameters) {
    if (is.null(text)) return(NULL)
    result <- matrix(0, nrow(text), ncol(text))
    for (i in seq_len(nrow(text))) {
        diag_value <- 1
        for (j in setdiff(seq_len(ncol(text)), i)) {
            result[i, j] <- parameters[text[i, j], "est"]
            diag_value <- 1 - result[i, j]
        }
        result[i, i] <- diag_value
    }
    result
}

.bgb_m3_layer_model <- function(layer, env, inputs, areas, sparse, trait, upper) {
    table <- inputs$BioGeoBEARS_model_object@params_table
    dispersal <- layer$distances ^ table["x", "est"]
    dispersal <- dispersal * env ^ table["n", "est"]
    dispersal <- dispersal * layer$manual ^ table["w", "est"]
    elist <- layer$area ^ table["u", "est"] * table["e", "est"]
    dmat <- dispersal * table["d", "est"]
    amat <- dispersal * table["a", "est"]
    inp <- get_spPmat_inputs_from_BGB(inputs, layer$states, dispersal)
    if (upper && !trait && !sparse && !is.null(layer$upper_loop_m))
        inp$dmat <- layer$upper_loop_m * dispersal
    empty <- numeric()
    q <- rows <- cols <- empty
    m <- jts <- NULL
    if (trait) {
        q_trait <- modify_Qmat_with_trait(Qmat = NULL, inputs, numstates_geogtrait = layer$width,
                                         areas_list = areas, states_list = layer$states,
                                         dispersal_multipliers_matrix = dispersal, elist = elist,
                                         force_sparse = sparse)
        m <- q_trait$m
        if (sparse) {
            rows <- q_trait$Qmat[, "ia"]
            cols <- q_trait$Qmat[, "ja"]
            q <- q_trait$Qmat[, "a"]
        } else {
            q <- q_trait$Qmat
        }
        if (!upper) jts <- .bgb_m3_jts(inputs$jts_txt_matrix, table)
        q_parameters <- list(rep(0, 4), empty, empty, empty)
    } else {
        q_parameters <- list(c(1, 1, sparse, 1e-26), dmat, elist, amat)
    }
    if (!trait) {
        prior <- .bgb_m2_priors(layer$prior_width, inp)
        na <- nrow(dispersal)
        priors <- lapply(prior$priors, function(probability) {
            result <- matrix(0, na, na)
            width <- min(na, nrow(probability))
            result[seq_len(width), seq_len(width)] <- probability[seq_len(width), seq_len(width), drop = FALSE]
            result
        })
        maxmin <- c(prior$maxmin, rep(0, max(0, na - length(prior$maxmin))))
        return(c(list(c(FALSE, sparse, 0, 0), empty, empty, empty, empty, empty,
                      empty, empty, empty, empty), q_parameters,
                 list(c(inp$s, inp$v, inp$j, inp$y), inp$dmat), priors, list(maxmin, na)))
    }
    if (upper) {
        prior <- .bgb_m2_priors(max(lengths(inp$l)), inp)
        ones <- rep(1, length(inp$l))
        coo <- rcpp_calc_anclikes_sp_COOweights_faster(ones, ones, inp$l, inp$s, inp$v, inp$j, inp$y,
                                                       dispersal, prior$priors[[1]], prior$priors[[2]],
                                                       prior$priors[[3]], prior$priors[[4]], prior$maxmin,
                                                       printmat = FALSE, m = m)
        denom <- rcpp_calc_rowsums_for_COOweights_columnar(coo)
        # The original upper trait expansion can emit fewer rows than the Q.
        # Only emitted events enter the Uppass sums.
        width <- layer$width - as.integer(inputs$include_null_range)
        denom <- c(denom, rep(0, max(0, width - length(denom))))
        prep <- list(anc = coo[[1]], left = coo[[2]], right = coo[[3]], weight = coo[[4]],
                     denom = denom)
    } else {
        prep <- .bgb_m2_clado(list(use_cpp = TRUE, cppSpMethod = 3, spPmat_inputs = inp,
                                   include_null_range = inputs$include_null_range, m = m, jts_matrix = jts),
                              layer$width)
    }
    clado_params <- if (is.null(prep$model)) rep(list(empty), 8L) else
        c(list(prep$model$parameters, prep$model$dmat), prep$model$priors, list(prep$model$maxmin, prep$model$na))
    c(list(c(trait && sparse, sparse, 0, length(m)), q, rows, cols,
           prep$anc, prep$left, prep$right, prep$weight, prep$denom, empty),
      q_parameters, clado_params)
}

.bgb_m3_eval <- function(inputs, phy, areas, states, tips, minimum, include_null, fixnode, fixlikes,
                         sparse, return_what, ancestors, return_table, ttl_from_table, rootedge, method) {
    contract <- .bgb_v12_runtime$contract
    root_spec <- .bgb_v12_runtime$root_prior
    contract_native <- (contract != 0L || !is.null(root_spec))
    if (identical(.bgb_stratified_engine(), "R")) {
        if (contract_native)
            stop("UnsupportedModelContractError: stratified probability contracts require native evaluation")
        return(NULL)
    }
    flags <- list(include_null, sparse, ancestors, return_table, ttl_from_table, rootedge)
    valid_flags <- all(vapply(flags, function(x) is.logical(x) && length(x) == 1L && !is.na(x), logical(1))) &&
        is.matrix(tips) && length(inputs$timeperiods) >= 2L && !is.null(inputs$master_table) &&
        identical(inputs$include_null_range, include_null) && identical(as.numeric(method), 3) &&
        !(ancestors && !ttl_from_table) &&
        !(identical(return_what, "all") && !(ttl_from_table || return_table))
    if (!valid_flags || (!contract_native && inputs$printlevel >= 3)) {
        if (contract_native)
            stop("UnsupportedModelContractError: this stratified configuration is outside the validated probability contract")
        return(NULL)
    }
    resolved_root <- NULL
    if (contract_native) {
        if (!is.null(inputs$trait) || length(fixnode) || !isTRUE(include_null) ||
            any(phy$edge.length < minimum) || any(isTRUE(inputs$master_table$fossils)))
            stop("UnsupportedModelContractError: stratified probability contracts reject traits, fixed internal nodes and fossils")
        p_tab <- inputs$BioGeoBEARS_model_object@params_table
        if (!.bgb_v12_standard_clado(p_tab[c("s", "v", "j", "y"), "est"]))
            stop("UnsupportedModelContractError: only the six standard cladogenesis models are validated")
        k_len <- length(inputs$timeperiods)
        masks <- lapply(seq_len(k_len), function(i) .bgb_m3_layer_mask(inputs, states, i, include_null))
        for (i in seq_len(k_len - 1L)) {
            if (any(masks[[i + 1L]] & !masks[[i]]))
                stop("UnsupportedModelContractError: probability contracts require monotone state-space emergence across time strata")
        }
        resolved_root <- .bgb_v12_resolve_root(root_spec, states, contract, allowed = masks[[k_len]])
        .bgb_v12_runtime$last_root <- resolved_root
    }
    runtime <- .bgb_m3_runtime
    precision <- .bgb_precision()
    threads <- as.integer(.bgb_m2_runtime$threads)
    keys <- c("timeperiods", "tree_sections_list", "master_table", "list_of_areas_allowed_mats",
              "list_of_areas_adjacency_mats", "lists_of_states_lists_0based", "list_of_distances_mats",
              "list_of_envdistances_mats", "list_of_dispersal_multipliers_mats", "list_of_area_of_areas",
              "trait", "trait_Pmat_txt", "jts_txt_matrix")
    signature <- list(inputs[keys], phy, areas, states, dim(tips), minimum, include_null,
                      fixnode, fixlikes, sparse, precision, threads, contract_native,
                      if (contract_native) resolved_root$likelihood_weights else NULL)
    evaluated <- tryCatch(withCallingHandlers({
        hit <- which(vapply(runtime$contexts, function(entry) identical(entry$signature, signature), logical(1)))
        if (length(hit)) {
            entry <- runtime$contexts[[hit[1L]]]
            runtime$contexts <- c(list(entry), runtime$contexts[-hit[1L]])
        } else {
            entry <- .bgb_m3_compile(inputs, phy, areas, states, ncol(tips), minimum, include_null, fixnode, fixlikes,
                                     exact_ages = contract_native,
                                     root_prior = if (contract_native) resolved_root$likelihood_weights else numeric())
            entry$signature <- signature
            entry$handle <- .Call(.bgb_c_bgb_strat_handle)
            if (length(runtime$contexts) == 4L)
                .Call(.bgb_c_bgb_tree_release, runtime$contexts[[4L]]$handle)
            runtime$contexts <- c(list(entry), head(runtime$contexts, 3L))
        }
        full <- identical(return_what, "all")
        jump <- inputs$BioGeoBEARS_model_object@params_table["j", "est"]
        models <- lapply(entry$layers, function(layer) {
            down <- .bgb_m3_layer_model(layer, layer$env_down, inputs, seq_along(areas) - 1L, sparse, entry$trait, FALSE)
            different_j <- !sparse && jump != 0 && !is.null(layer$upper_loop_m)
            separate <- !contract_native && ancestors && (entry$trait || different_j || !identical(layer$env_down, layer$env_up))
            if (separate) {
                up <- .bgb_m3_layer_model(layer, layer$env_up, inputs, seq_along(areas) - 1L, sparse, entry$trait, TRUE)
                c(list(1), down, up)
            } else {
                c(list(0), down)
            }
        })
        buffers <- c(list(c(ancestors, full, return_table, ttl_from_table, rootedge), tips),
                     unlist(models, recursive = FALSE))
        native <- if (contract_native)
            .Call(.bgb_c_bgb_tree_contract, entry$handle, entry$config, buffers, precision, threads, contract, numeric())
        else
            .Call(.bgb_c_bgb_tree, entry$handle, entry$config, buffers, precision, threads)
        list(native = native, entry = entry, full = full)
    }, warning = function(w) stop(conditionMessage(w), call. = FALSE)), error = identity)
    if (inherits(evaluated, "error")) {
        if (contract_native) stop(evaluated)
        runtime$fallbacks <- runtime$fallbacks + 1L
        runtime$last_fallback <- conditionMessage(evaluated)
        return(NULL)
    }
    .bgb_v12_runtime$last_log_survival <- if (contract == 1L) tail(evaluated$native[[1]], 1L) else NULL
    runtime$native_calls <- runtime$native_calls + 1L
    native <- evaluated$native
    entry <- evaluated$entry
    lnL <- native[[1]][1]
    if (ancestors) {
        cat("\nUppass started for (STRATIFIED) marginal ancestral states estimation!\n")
        if (!is.null(inputs$list_of_areas_allowed_mats)) {
            for (i in rev(seq_along(inputs$timeperiods))) {
                cat("i=", i, "\n")
                print(inputs$list_of_areas_allowed_mats[[i]])
            }
        }
        cat("\nUppass completed for (STRATIFIED) marginal ancestral states estimation!\n")
    }
    if (!evaluated$full) return(list(value = lnL))
    n <- ncol(tips)
    result <- list(final_all_condlikes_of_each_state = matrix(native[[2]], native[[1]][2], n),
                   condlikes_table = matrix(native[[3]], entry$nr, n))
    if (ancestors) {
        result$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE <- matrix(native[[4]], entry$nr, n)
        result$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS <- matrix(native[[5]], entry$nn, n)
        result$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS <- matrix(native[[6]], entry$nn, n)
        result$ML_marginal_prob_each_state_at_branch_bottom_below_node <- matrix(native[[7]], entry$nn, n)
        result$ML_marginal_prob_each_state_at_branch_top_AT_node <- matrix(native[[8]], entry$nn, n)
    }
    result$grand_total_likelihood <- lnL
    result$total_loglikelihood <- lnL
    if (ancestors) {
        bottom <- matrix(native[[9]], entry$nn, n)
        rownames(bottom) <- replace(rep("", entry$nn), entry$nt + 1L, "root_row")
        result$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS <- bottom
        rownames(result$ML_marginal_prob_each_state_at_branch_bottom_below_node) <- rownames(bottom)
        result$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS <- matrix(native[[10]], entry$nn, n)
    }
    list(value = result)
}
