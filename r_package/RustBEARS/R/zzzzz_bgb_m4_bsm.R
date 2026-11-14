# GPL >= 2. R owns RNG scope, compatibility fallbacks and table formatting.
.bgb_m4_runtime <- new.env(parent = emptyenv())
.bgb_m4_runtime$active <- NULL
.bgb_m4_runtime$legacy_branches <- 0L
.bgb_m4_runtime$legacy_continuations <- 0L

.bgb_m4_bottom_age <- function(trtable, row, stratified, brlen) {
    if (stratified == TRUE) trtable$time_top[row] + trtable$SUBtime_bp[row] + brlen
    else trtable$time_bp[row] + brlen
}

.bgb_m4_legacy_branch <- function(q, n, from, to, duration, maxtries) {
    if (identical(.bgb_bsm_engine(), "R") ||
        !is.matrix(q) || !is.double(q) || !identical(dim(q), c(n, n)) ||
        length(from) != 1L || length(to) != 1L || anyNA(c(from, to)) ||
        from < 1L || from > n || to < 1L || to > n || from != trunc(from) || to != trunc(to) ||
        length(maxtries) != 1L || !is.finite(maxtries) || maxtries < 1 ||
        maxtries > .Machine$integer.max || maxtries != trunc(maxtries) ||
        !exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) ||
        RNGkind()[1L] == "user-supplied") return(NULL)
    if (is.na(q[from, from])) return(NULL)
    if (q[from, from] == 0 && from == to) {
        .bgb_m4_runtime$legacy_branches <- .bgb_m4_runtime$legacy_branches + 1L
        return(list(c(1, 1), numeric()))
    }
    if (q[from, from] >= 0) return(NULL)
    attempts <- 0
    while (attempts < maxtries) {
        result <- .Call(.bgb_c_bgb_bsm_legacy, q, n, from, to, duration, min(maxtries - attempts, 512))
        if (is.null(result)) return(NULL)
        attempts <- attempts + result[[1L]][1L]
        if (result[[1L]][2L] == 1) {
            result[[1L]][1L] <- attempts
            .bgb_m4_runtime$legacy_branches <- .bgb_m4_runtime$legacy_branches + 1L
            return(result)
        }
    }
    .bgb_m4_runtime$legacy_continuations <- .bgb_m4_runtime$legacy_continuations + 1L
    # The last upstream attempt owns its diagnostic, manual synthesis and
    # try-error behaviour. Earlier failed attempts have advanced R's RNG.
    list(c(attempts, 0), numeric())
}

.bgb_m4_event_row <- function(nodenum_at_top_of_branch, trynum, brlen, current_rangenum_1based,
                              new_rangenum_1based, abs_event_time, event_time, states_list,
                              ranges_list, areas) {
if ( length(states_list[[new_rangenum_1based]]) > length(states_list[[current_rangenum_1based]]) )
	{
	event_type = "d"
	TF = states_list[[new_rangenum_1based]] %in% states_list[[current_rangenum_1based]]
	new_area_TF = TF == FALSE
	new_area_num_0based = states_list[[new_rangenum_1based]][new_area_TF]
	new_area_num_1based = new_area_num_0based + 1
	lost_area_num_1based = "-"
	dispersal_to = areas[new_area_num_1based]
	extirpation_from = "-"
	}
if ( length(states_list[[new_rangenum_1based]]) < length(states_list[[current_rangenum_1based]]) )
	{
	event_type = "e"
	TF = states_list[[current_rangenum_1based]] %in% states_list[[new_rangenum_1based]]
	lost_area_TF = TF == FALSE
	lost_area_num_0based = states_list[[current_rangenum_1based]][lost_area_TF]
	lost_area_num_1based = lost_area_num_0based + 1
	new_area_num_1based = "-"
	dispersal_to = "-"
	extirpation_from = areas[lost_area_num_1based]
	}
if ( length(states_list[[new_rangenum_1based]]) == length(states_list[[current_rangenum_1based]]) )
	{
	if (ranges_list[[new_rangenum_1based]] == "_")
		{
		event_type = "e"
		TF = states_list[[current_rangenum_1based]] %in% states_list[[new_rangenum_1based]]
		lost_area_TF = TF == FALSE
		lost_area_num_0based = states_list[[current_rangenum_1based]][lost_area_TF]
		lost_area_num_1based = lost_area_num_0based + 1
		new_area_num_1based = "-"
		dispersal_to = "-"
		extirpation_from = areas[lost_area_num_1based]
		} else {
		event_type = "a" # anagenetic range-switching
		new_area_num_1based = 1+states_list[[new_rangenum_1based]]
		dispersal_to = areas[new_area_num_1based]
		lost_area_num_1based = 1+states_list[[current_rangenum_1based]]
		extirpation_from = areas[lost_area_num_1based]
		}
	} # END if ( length(states_list[[new_rangenum_1based]]) == 
event_txt = paste(ranges_list[[current_rangenum_1based]], "->", ranges_list[[new_rangenum_1based]], sep="")
current_rangetxt = ranges_list[[current_rangenum_1based]]
new_rangetxt = ranges_list[[new_rangenum_1based]]
tmprow = c(nodenum_at_top_of_branch, trynum, brlen, current_rangenum_1based, new_rangenum_1based, current_rangetxt, new_rangetxt, abs_event_time, event_time, event_type, event_txt, new_area_num_1based, lost_area_num_1based, dispersal_to, extirpation_from)
    tmprow
}

.bgb_m4_event_table <- function(events, row, attempts, duration, bottom_age, states, ranges, areas) {
    if (!length(events)) return(NA)
    history <- matrix(events, ncol = 3L, byrow = TRUE)
    rows <- lapply(seq_len(nrow(history)), function(i) {
        event <- history[i, ]
        .bgb_m4_event_row(row, attempts, duration, event[1L], event[2L], bottom_age - event[3L],
                          event[3L], states, ranges, areas)
    })
    table <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(table) <- c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based",
                      "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time",
                      "event_time", "event_type", "event_txt", "new_area_num_1based",
                      "lost_area_num_1based", "dispersal_to", "extirpation_from")
    row.names(table) <- NULL
    table
}

.bgb_m4_seed <- function(seed) {
    if (is.null(.bgb_m4_runtime$active)) set.seed(seed)
}

.bgb_m4_set_stratum <- function(stratum) {
    state <- .bgb_m4_runtime$active
    if (!is.null(state)) state$stratum <- stratum
}

.bgb_m4_sample <- function(x, size, replace, prob) {
    state <- .bgb_m4_runtime$active
    if (is.null(state)) return(sample(x, size, replace, prob))
    if (size != 1 || length(x) != length(prob) || any(!is.finite(prob)) || any(prob < 0))
        stop("invalid BSM endpoint probabilities", call. = FALSE)
    total <- sum(prob)
    if (total <= 0 || !is.finite(total)) stop("BSM endpoint has no positive probability", call. = FALSE)
    state$step <- state$step + 1
    u <- .bgb_kernel(25L, list(c(state$seed, state$map, state$stratum, state$step)))[[1L]]
    positive <- which(prob > 0)
    index <- match(TRUE, cumsum(prob[positive]) >= u * total, nomatch = length(positive))
    x[positive[index]]
}

.bgb_m4_record_branch <- function(row, trtable, q, states, ranges, areas, stratified, duration, from, to) {
    state <- .bgb_m4_runtime$active
    id <- match(TRUE, vapply(state$generators, identical, logical(1), q), nomatch = 0L)
    if (id == 0L) {
        id <- length(state$generators) + 1L
        state$generators[[id]] <- q
    }
    state$step <- state$step + 1
    i <- length(state$jobs) + 1L
    state$jobs[[i]] <- list(
        args = c(id, from, to, duration, state$map, state$stratum, state$step, trtable$node[row]),
        row = row, duration = duration, bottom = .bgb_m4_bottom_age(trtable, row, stratified, duration),
        states = states, ranges = ranges, areas = areas)
    structure(i, class = "bgb_bsm_pending")
}

.bgb_m4_validate_run <- function(seed, threads, tries, goal) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 || seed >= 2^53)
        stop("seedval must be a nonnegative number below 2^53", call. = FALSE)
    count_ok <- function(x, upper) is.numeric(x) && length(x) == 1L && is.finite(x) &&
        x >= 1 && x <= upper && x == trunc(x)
    if (!count_ok(threads, 256)) stop("num_threads must be an integer from 1 to 256", call. = FALSE)
    if (!count_ok(tries, .Machine$integer.max %/% 2L) || !count_ok(goal, .Machine$integer.max %/% 2L))
        stop("modern BSM map counts must be positive integers", call. = FALSE)
}

.bgb_m4_fill_events <- function(table, text) {
    column <- table$anagenetic_events_txt_below_node
    pending <- startsWith(column, ".bgb_m4_pending_")
    pending[is.na(pending)] <- FALSE
    if (any(pending)) {
        index <- as.integer(substring(column[pending], nchar(".bgb_m4_pending_") + 1L))
        table$anagenetic_events_txt_below_node[pending] <- text[index]
    }
    table
}

.bgb_m4_map_batch <- function(res, inputs, stratified, maps, seed, threads, maxtries, print_node) {
    previous <- .bgb_m4_runtime$active
    on.exit({ .bgb_m4_runtime$active <- previous })
    state <- new.env(parent = emptyenv())
    state$jobs <- state$generators <- list()
    state$seed <- floor(seed)
    .bgb_m4_runtime$active <- state
    endpoint_tables <- lapply(maps, function(map) {
        state$map <- map
        state$stratum <- state$step <- 0
        saved <- length(state$jobs)
        result <- try(if (stratified) {
            stochastic_mapping_on_stratified(res, inputs, maxtries, seed + map, print_node)
        } else {
            stochastic_map_given_inputs(inputs, maxtries = maxtries, seedval = seed + map,
                include_null_range = res$inputs$include_null_range,
                master_nodenum_toPrint = print_node, allow_null_tips = res$inputs$allow_null_tips)
        }, silent = TRUE)
        if (inherits(result, "try-error")) state$jobs <- head(state$jobs, saved)
        result
    })
    .bgb_m4_runtime$active <- previous
    request <- unlist(lapply(state$jobs, function(job) job$args), use.names = FALSE)
    if (length(request) != length(state$jobs) * 8) stop("malformed BSM branch request", call. = FALSE)
    histories <- tryCatch(.bgb_kernel(24L, c(list(c(state$seed, threads), request), state$generators)),
        error = function(error) {
            stop(structure(list(message = conditionMessage(error), call = NULL),
                           class = c("bgb_bsm_bridge_error", "error", "condition")))
        })
    text <- vapply(seq_along(histories), function(i) {
        job <- state$jobs[[i]]
        table <- .bgb_m4_event_table(histories[[i]], job$row, 1L, job$duration, job$bottom,
                                    job$states, job$ranges, job$areas)
        events_table_into_txt(table)
    }, character(1))
    lapply(endpoint_tables, function(result) {
        if (inherits(result, "try-error")) return(result)
        if (!stratified) return(.bgb_m4_fill_events(result, text))
        table <- .bgb_m4_fill_events(result$master_table_cladogenetic_events, text)
        result$master_table_cladogenetic_events <- table
        result$table_w_anagenetic_events <- events_txt_list_into_events_table(
            events_txt_list = table$anagenetic_events_txt_below_node, trtable = table,
            trtable_rownums_of_events_txt = seq_len(nrow(table)))
        result
    })
}

.bgb_m4_transition <- function(q, duration, cache) {
    hit <- match(duration, cache$lengths, nomatch = 0L)
    if (hit > 0L) return(cache$matrices[[hit]])
    p <- expokit_dgpadm_Qmat2(times = duration, Qmat = q, transpose_needed = TRUE)
    if (length(cache$lengths) == 32L) {
        cache$lengths <- cache$lengths[-1L]
        cache$matrices <- cache$matrices[-1L]
    }
    cache$lengths <- c(cache$lengths, duration)
    cache$matrices[[length(cache$lengths)]] <- p
    p
}

.bgb_m4_independent <- function(phy, q, cluster, cache) {
    if (identical(cluster, FALSE)) cluster <- NULL
    times <- reorder(phy, "pruningwise")$edge.length
    if (!is.null(cluster) || identical(.bgb_bsm_engine(), "R") ||
        !length(times) || !is.null(names(times)) || any(!is.finite(times)))
        return(calc_independent_likelihoods_on_each_branch(phy, q,
                   cluster_already_open = cluster, Qmat_is_sparse = FALSE))
    n <- nrow(q)
    matrices <- lapply(times, function(duration) .bgb_m4_transition(q, duration, cache))
    array(unlist(matrices, use.names = FALSE), c(n, n, length(times)))
}
