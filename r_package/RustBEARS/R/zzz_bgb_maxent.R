# SPDX-License-Identifier: GPL-3.0-or-later
# Clean-room R interface wrapper for discrete maximum-entropy prior solver
# (Della Pietra, Della Pietra & Lafferty 1997, IEEE TPAMI 19:380-393).
# Authored from scratch for RustBEARS without FD/rexpokit source code so the
# combined package remains cleanly licensed under GPL-3.0-or-later.

.rustbears_as_2d_numeric <- function(x, default_rownames = NULL, label = "input") {
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    } else if (!is.matrix(x)) {
        nm <- names(x)
        x <- matrix(x, nrow = 1L, ncol = length(x))
        dimnames(x) <- list(default_rownames, nm)
    }
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
        stop(sprintf("%s must contain finite numeric values without NAs", label))
    }
    storage.mode(x) <- "double"
    x
}

maxent <- function(constr, states, prior, tol = 1e-07, lambda = FALSE) {
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("tol must be a single positive finite number")
    }
    if (!is.logical(lambda) || length(lambda) != 1L || is.na(lambda)) {
        stop("lambda must be TRUE or FALSE")
    }

    target_mat <- .rustbears_as_2d_numeric(constr, default_rownames = "set1", label = "constr")
    feature_mat <- .rustbears_as_2d_numeric(states, default_rownames = "constraint", label = "states")
    if (ncol(feature_mat) == 1L && nrow(feature_mat) > 1L) {
        feature_mat <- t(feature_mat)
    }
    if (any(feature_mat <= 0)) {
        stop("states must contain strictly positive values")
    }

    num_features <- nrow(feature_mat)
    num_categories <- ncol(feature_mat)
    num_sets <- nrow(target_mat)
    if (ncol(target_mat) != num_features) {
        stop("ncol(constr) must match the number of constraint rows in states")
    }

    cat_labels <- colnames(feature_mat)
    feat_labels <- rownames(feature_mat)
    set_labels <- rownames(target_mat)

    if (missing(prior) || is.null(prior)) {
        prior_mat <- matrix(
            1.0 / num_categories,
            nrow = num_sets,
            ncol = num_categories,
            dimnames = list(set_labels, cat_labels)
        )
    } else {
        if (!is.matrix(prior) && !is.data.frame(prior)) {
            if (!is.numeric(prior) || anyNA(prior) || length(prior) != num_categories) {
                stop("vector prior must be finite numeric with length matching states")
            }
            prior_mat <- matrix(rep(as.double(prior), each = num_sets), nrow = num_sets, ncol = num_categories)
        } else {
            prior_mat <- .rustbears_as_2d_numeric(prior, label = "prior")
            if (ncol(prior_mat) == 1L && nrow(prior_mat) > 1L) {
                prior_mat <- t(prior_mat)
            }
            if (ncol(prior_mat) != num_categories) {
                stop("prior column count must match number of categories in states")
            }
            if (nrow(prior_mat) == 1L && num_sets > 1L) {
                prior_mat <- prior_mat[rep(1L, num_sets), , drop = FALSE]
            } else if (nrow(prior_mat) != num_sets) {
                stop("prior must have either 1 row or nrow(constr) rows")
            }
        }
        if (any(prior_mat < 0) || any(prior_mat > 1)) {
            stop("prior entries must lie within [0, 1]")
        }
        row_totals <- rowSums(prior_mat)
        if (any(row_totals <= 0)) {
            stop("each prior row must have positive total weight")
        }
        prior_mat <- prior_mat / row_totals
        dimnames(prior_mat) <- list(set_labels, cat_labels)
    }

    prob_out <- matrix(0.0, nrow = num_sets, ncol = num_categories, dimnames = list(set_labels, cat_labels))
    moment_out <- matrix(0.0, nrow = num_sets, ncol = num_features, dimnames = list(set_labels, feat_labels))
    entropy_out <- stats::setNames(numeric(num_sets), set_labels)
    iter_out <- stats::setNames(integer(num_sets), set_labels)
    lambda_out <- if (lambda) {
        matrix(
            0.0,
            nrow = num_sets,
            ncol = num_features + 2L,
            dimnames = list(set_labels, c("intercept", feat_labels, "prior"))
        )
    } else {
        NULL
    }

    for (row_idx in seq_len(num_sets)) {
        solved <- .bgb_itscale(
            feature_mat,
            num_categories,
            num_features,
            target_mat[row_idx, ],
            prior_mat[row_idx, ],
            tol
        )
        prob_out[row_idx, ] <- solved$prob
        moment_out[row_idx, ] <- solved$moments
        entropy_out[row_idx] <- solved$entropy
        iter_out[row_idx] <- solved$niter
        if (lambda) {
            design <- cbind(1.0, t(feature_mat), log(prior_mat[row_idx, ]))
            lambda_out[row_idx, ] <- stats::lm.fit(x = design, y = log(solved$prob))$coefficients
        }
    }

    if (num_sets == 1L) {
        out <- list(
            prob = prob_out[1L, ],
            moments = moment_out[1L, ],
            entropy = unname(entropy_out[1L]),
            iter = unname(iter_out[1L])
        )
        if (lambda) {
            out$lambda <- lambda_out[1L, ]
        }
        out$constr <- target_mat[1L, ]
        out$states <- feature_mat
        out$prior <- prior_mat[1L, ]
        return(out)
    }

    out <- list(
        prob = prob_out,
        moments = moment_out,
        entropy = entropy_out,
        iter = iter_out
    )
    if (lambda) {
        out$lambda <- lambda_out
    }
    out$constr <- target_mat
    out$states <- feature_mat
    out$prior <- prior_mat
    out
}
