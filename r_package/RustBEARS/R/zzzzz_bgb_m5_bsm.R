# GPL >= 2. One run owns its deterministic, non-RNG BSM cache.
.bgb_m5_bsm_runtime <- new.env(parent = emptyenv())
.bgb_m5_bsm_runtime$scope <- NULL

.bgb_m5_plain_transition <- function(p) {
    is.double(p) && !is.object(p) && length(dim(p)) %in% c(2L, 3L) &&
        is.null(dimnames(p)) && all(is.finite(p)) && all(p >= 0)
}

.bgb_m5_mark_transitions <- function(inputs, stratified) {
    if (stratified) {
        lapply(inputs, function(layer) {
            if (is.list(layer)) layer$.bgb_m5_fast_pieces <- vapply(
                layer$independent_likelihoods_by_tree_piece_for_timeperiod_i,
                .bgb_m5_plain_transition, logical(1))
            layer
        })
    } else {
        inputs$.bgb_m5_fast_edges <- .bgb_m5_plain_transition(inputs$independent_likelihoods_on_each_branch)
        inputs
    }
}

.bgb_m5_one_hot <- function(v, n) {
    if (length(v) != n || anyNA(v)) return(integer())
    at <- which(v != 0)
    if (length(at) == 1L && v[at] == 1) at else integer()
}

.bgb_m5_row_product <- function(v, p, fast) {
    if (isTRUE(fast)) {
        at <- .bgb_m5_one_hot(v, nrow(p))
        if (length(at)) {
            row <- p[at, , drop = FALSE]
            row[row == 0] <- 0
            return(row)
        }
    }
    v %*% p
}

.bgb_m5_edge_product <- function(v, p, edge, fast) {
    if (isTRUE(fast)) {
        at <- .bgb_m5_one_hot(v, dim(p)[1L])
        if (length(at)) {
            row <- matrix(p[at, , edge], nrow = 1L)
            row[row == 0] <- 0
            return(row)
        }
    }
    v %*% p[, , edge]
}

.bgb_m5_root_transition <- function(q, length) {
    scope <- .bgb_m5_bsm_runtime$scope
    if (is.null(scope)) {
        p <- expokit_dgpadm_Qmat2(times = length, Qmat = q, transpose_needed = TRUE)
        return(list(p = p, fast = FALSE))
    }
    key <- list(q, length, .bgb_precision(), getOption("matprod"))
    hit <- which(vapply(scope$roots, function(x) identical(x$key, key, num.eq=FALSE), logical(1)))
    if (length(hit)) return(scope$roots[[hit[1L]]])
    p <- expokit_dgpadm_Qmat2(times = length, Qmat = q, transpose_needed = TRUE)
    value <- list(key = key, p = p, fast = .bgb_m5_plain_transition(p))
    scope$roots <- c(list(value), head(scope$roots, 31L))
    value
}

.bgb_m5_clado_ranges <- function(areas, states, maxareas, include_null_range) {
    scope <- .bgb_m5_bsm_runtime$scope
    key <- list(areas, states, maxareas, include_null_range, .bgb_precision())
    if (!is.null(scope)) {
        hit <- which(vapply(scope$ranges, function(x) identical(x$key, key), logical(1)))
        if (length(hit)) return(scope$ranges[[hit[1L]]])
    }
    if (is.null(states)) {
        states <- rcpp_areas_list_to_states_list(areas = areas, maxareas = maxareas,
                                                include_null_range = include_null_range)
        ranges <- areas_list_to_states_list_new(areas = areas, maxareas = maxareas,
                                                include_null_range = include_null_range, split_ABC = FALSE)
    } else {
        ranges <- lapply(unname(states), function(state) {
            if (length(state) == 1L && is.na(state)) "_" else paste(areas[1 + state], collapse = "")
        })
    }
    value <- list(key = key, states = states, ranges = ranges)
    if (!is.null(scope)) scope$ranges <- c(list(value), head(scope$ranges, 7L))
    value
}
