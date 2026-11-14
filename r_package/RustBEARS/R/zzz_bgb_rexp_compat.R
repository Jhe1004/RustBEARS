# GPL >= 2. rexpokit 0.26.6.9 R compatibility shell, Nick Matzke.
# Native numerical calls replaced with Rust, helper names kept private.
coo2mat <- function (coomat, n = max(max(coomat[, 1]), max(coomat[, 2])), transpose_needed = FALSE) 
{
    defaults = "\n\t\n\t"
    outmat = matrix(double(length = n * n), nrow = n)
    ia = coomat[, 1]
    ja = coomat[, 2]
    a = coomat[, 3]
    if (transpose_needed == FALSE) {
        for (k in 1:length(ia)) {
            outmat[ia[k], ja[k]] = a[k]
        }
    }
    else {
        for (k in 1:length(ia)) {
            outmat[ja[k], ia[k]] = a[k]
        }
    }
    return(outmat)
}

expm <- function (x, t = 1) 
{
    if (nrow(x) != ncol(x)) 
        stop("Matrix exponentiation is only defined for square matrices.")
    ret <- expokit_dgpadm_Qmat(Qmat = x, t = t, transpose_needed = TRUE)
    return(ret)
}

expokit_dgexpv_Qmat <- function (Qmat = NULL, t = 2.1, inputprobs_for_fast = NULL, transpose_needed = TRUE, transform_to_coo_TF = TRUE, coo_n = NULL, anorm = NULL, check_for_0_rows = TRUE) 
{
    defaults = "\n\tQmat=NULL\n\tt = 2.1\n\tinputprobs_for_fast=NULL\n\ttranspose_needed=TRUE\n\ttransform_to_coo_TF=TRUE\n\tcoo_n=NULL\n\tanorm=NULL\n\tcheck_for_0_rows=FALSE\n\tcheck_for_0_rows=1e-15\n\t"
    matvec = Qmat
    if (is.null(matvec)) {
        warning("expokit_dgexpv_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead")
        matvec = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow = 4, byrow = TRUE)
    }
    if (is.null(t)) {
        stop("\nSTOP ERROR: expokit_dgexpv_Qmat() was provided a t (time or times list) with value NULL.  \n")
    }
    if (is.null(inputprobs_for_fast)) {
        if (.bgb_precision() != 1L && check_for_0_rows != FALSE) {
            rows_w_all_zeros_TF = .bgb_rexp_findrows_w_all_zeros(matvec)
            if (all(rows_w_all_zeros_TF == FALSE)) {
                pass = 1
            }
            else {
                rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]
                if (is.numeric(check_for_0_rows)) {
                  check_for_0_rows = check_for_0_rows
                }
                else {
                  check_for_0_rows = 1e-15
                }
                newrowvals = rep(check_for_0_rows, ncol(matvec))
                matvec[rows_allzero_indices, ] = newrowvals
                diagonal_val = -1 * sum(newrowvals[-1])
                matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
                warning(paste(sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldgexpv_', ...)\nand therefore expokit_wrapalldgexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep = ""))
            }
        }
    }
    if (transform_to_coo_TF == TRUE) {
        nz = sum(matvec != 0)
        ia = integer(length = nz)
        ja = integer(length = nz)
        a = double(length = nz)
        n = nrow(matvec)
    }
    else {
        n = coo_n
        nz = sum(matvec[, "a"] != 0)
    }
    ideg = as.integer(6)
    m = n - 1
    if (is.null(inputprobs_for_fast)) {
        v = double(n)
        v[1] = 1
    }
    else {
        v = double(n)
        v = inputprobs_for_fast
    }
    w = double(length = n)
    tol = as.numeric(0.01)
    lwsp = as.integer(n * (m + 2) + 5 * (m + 2)^2 + ideg + 1)
    wsp = double(length = lwsp)
    liwsp = max(m + 2, 7)
    iwsp = integer(length = liwsp)
    matvec = matvec
    if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE)) {
        tmatvec = t(matvec)
    }
    else {
        tmatvec = matvec
    }
    if ((exists("anorm") == FALSE) || is.null(anorm)) {
        if (transform_to_coo_TF == FALSE && transpose_needed == FALSE) {
            tmpQmat1 = coo2mat(matvec, n = coo_n)
            tmpQmat2 = t(tmpQmat1)
            anorm = as.numeric(norm(tmpQmat2, type = "O"))
        }
        else {
            anorm = as.numeric(norm(matvec, type = "O"))
        }
    }
    itrace = 0
    iflag = 0
    if (transform_to_coo_TF == TRUE) {
        tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
    }
    else {
        tmpmat_in_REXPOKIT_coo_fmt = matvec
    }
    ia = tmpmat_in_REXPOKIT_coo_fmt[, "ia"]
    ja = tmpmat_in_REXPOKIT_coo_fmt[, "ja"]
    a = tmpmat_in_REXPOKIT_coo_fmt[, "a"]
    if (is.null(inputprobs_for_fast)) {
        ret <- .bgb_rexp_call("R_dgexpv", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
        output_Pmat = matrix(ret$res, nrow = n, byrow = TRUE)
        return(output_Pmat)
    }
    else {
        v = inputprobs_for_fast
        ret <- .bgb_rexp_call("R_dgexpv", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
        w_output_probs = ret
        return(w_output_probs)
    }
}

expokit_dgpadm_Qmat <- function (Qmat = NULL, t = 2.1, transpose_needed = TRUE) 
{
    defaults = "\n\tQmat=NULL\n\tt = 2.1\n\ttranspose_needed=TRUE\n\t"
    if (is.null(Qmat)) {
        warning("expokit_dgpadm_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead")
        Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow = 4, byrow = TRUE)
    }
    if (is.null(t)) {
        stop("\nSTOP ERROR: expokit_dgpadm_Qmat() was provided a t (time or times list) with value NULL.  \n")
    }
    ideg = as.integer(6)
    m = as.integer(nrow(Qmat))
    res = double(length = m * m)
    matvec = Qmat
    if (transpose_needed == TRUE) {
        tmatvec = t(matvec)
        H = as.numeric(tmatvec)
    }
    else {
        H = as.numeric(matvec)
    }
    ldh = m
    if (!is.double(H)) 
        storage.mode(H) <- "double"
    res <- .bgb_rexp_call("R_dgpadm", as.integer(ideg), as.integer(m), as.double(t), H, as.integer(ldh))
    output = res$wsp
    ind <- res$ind
    output_Pmat_is = seq(ind, ind + m * m - 1, by = 1)
    output_Pmat = output[output_Pmat_is]
    output_Pmat = matrix(output_Pmat, nrow = m, byrow = TRUE)
    return(output_Pmat)
}

expokit_dmexpv_Qmat <- function (Qmat = NULL, t = 2.1, inputprobs_for_fast = NULL, transpose_needed = TRUE, transform_to_coo_TF = TRUE, coo_n = NULL, anorm = NULL, check_for_0_rows = TRUE) 
{
    defaults = "\n\tQmat=NULL\n\tt = 2.1\n\tinputprobs_for_fast=NULL\n\ttranspose_needed=TRUE\n\ttransform_to_coo_TF=TRUE\n\tcoo_n=NULL\n\tanorm=NULL\n\tcheck_for_0_rows=FALSE\n\tcheck_for_0_rows=1e-15\n\t"
    matvec = Qmat
    if (is.null(matvec)) {
        warning("expokit_dmexpv_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead")
        matvec = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow = 4, byrow = TRUE)
    }
    if (is.null(t)) {
        stop("\nSTOP ERROR: expokit_dmexpv_Qmat() was provided a t (time or times list) with value NULL.  \n")
    }
    if (is.null(inputprobs_for_fast)) {
        if (.bgb_precision() != 1L && check_for_0_rows != FALSE) {
            rows_w_all_zeros_TF = .bgb_rexp_findrows_w_all_zeros(matvec)
            if (all(rows_w_all_zeros_TF == FALSE)) {
                pass = 1
            }
            else {
                rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]
                if (is.numeric(check_for_0_rows)) {
                  check_for_0_rows = check_for_0_rows
                }
                else {
                  check_for_0_rows = 1e-15
                }
                newrowvals = rep(check_for_0_rows, ncol(matvec))
                matvec[rows_allzero_indices, ] = newrowvals
                diagonal_val = -1 * sum(newrowvals[-1])
                matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
                warning(paste(sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldmexpv_', ...)\nand therefore expokit_wrapalldmexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep = ""))
            }
        }
    }
    if (transform_to_coo_TF == TRUE) {
        nz = sum(matvec != 0)
        ia = integer(length = nz)
        ja = integer(length = nz)
        a = double(length = nz)
        n = nrow(matvec)
    }
    else {
        n = coo_n
        nz = sum(matvec[, "a"] != 0)
    }
    ideg = as.integer(6)
    m = n - 1
    if (is.null(inputprobs_for_fast)) {
        v = double(n)
        v[1] = 1
    }
    else {
        v = double(n)
        v = inputprobs_for_fast
    }
    tol = as.numeric(0.01)
    lwsp = as.double(n * (m + 2) + 5 * (m + 2)^2 + ideg + 1)
    wsp = double(length = lwsp)
    liwsp = max(m + 2, 7)
    iwsp = integer(length = liwsp)
    matvec = matvec
    if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE)) {
        tmatvec = t(matvec)
    }
    else {
        tmatvec = matvec
    }
    if ((exists("anorm") == FALSE) || is.null(anorm)) {
        if (transform_to_coo_TF == FALSE && transpose_needed == FALSE) {
            tmpQmat1 = coo2mat(matvec, n = coo_n)
            tmpQmat2 = t(tmpQmat1)
            anorm = as.numeric(norm(tmpQmat2, type = "O"))
        }
        else {
            anorm = as.numeric(norm(matvec, type = "O"))
        }
    }
    if (transform_to_coo_TF == TRUE) {
        tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
    }
    else {
        tmpmat_in_REXPOKIT_coo_fmt = matvec
    }
    ia = tmpmat_in_REXPOKIT_coo_fmt[, "ia"]
    ja = tmpmat_in_REXPOKIT_coo_fmt[, "ja"]
    a = tmpmat_in_REXPOKIT_coo_fmt[, "a"]
    if (is.null(inputprobs_for_fast)) {
        output_Pmat = .bgb_rexp_expokit_dmexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
        return(output_Pmat)
    }
    else {
        v = inputprobs_for_fast
        w_output_probs = .bgb_rexp_expokit_mydmexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
        return(w_output_probs)
    }
}

expokit_wrapalldgexpv_tvals <- function (Qmat = NULL, tvals = c(2.1), inputprobs_for_fast = NULL, transpose_needed = TRUE, transform_to_coo_TF = TRUE, coo_n = NULL, force_list_if_1_tval = FALSE, check_for_0_rows = TRUE) 
{
    defaults = "\n\tQmat=NULL\n\tt = 2.1\n\tinputprobs_for_fast=NULL\n\ttranspose_needed=TRUE\n\tCOO_needed=TRUE\n\ttransform_to_coo_TF=TRUE\n\tcoo_n=NULL\n\tforce_list_if_1_tval=FALSE\n\tcheck_for_0_rows=FALSE\n\tcheck_for_0_rows=1e-15\n\t"
    if (is.null(Qmat)) {
        warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
        Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow = 4, byrow = TRUE)
    }
    if (is.null(tvals)) {
        warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
        tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
    }
    matvec = Qmat
    if (is.null(inputprobs_for_fast)) {
        if (.bgb_precision() != 1L && check_for_0_rows != FALSE) {
            rows_w_all_zeros_TF = .bgb_rexp_findrows_w_all_zeros(matvec)
            if (all(rows_w_all_zeros_TF == FALSE)) {
                pass = 1
            }
            else {
                rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]
                if (is.numeric(check_for_0_rows)) {
                  check_for_0_rows = check_for_0_rows
                }
                else {
                  check_for_0_rows = 1e-15
                }
                newrowvals = rep(check_for_0_rows, ncol(matvec))
                matvec[rows_allzero_indices, ] = newrowvals
                diagonal_val = -1 * sum(newrowvals[-1])
                matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
                warning(paste(sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldgexpv_', ...)\nand therefore expokit_wrapalldgexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep = ""))
            }
        }
    }
    if (transform_to_coo_TF == TRUE) {
        nz = sum(matvec != 0)
        ia = integer(length = nz)
        ja = integer(length = nz)
        a = double(length = nz)
        n = nrow(matvec)
    }
    else {
        n = coo_n
        nz = sum(matvec[, "a"] != 0)
    }
    ideg = as.integer(6)
    m = n - 1
    if (is.null(inputprobs_for_fast)) {
        v = double(n)
        v[1] = 1
    }
    else {
        v = double(n)
        v = inputprobs_for_fast
    }
    w = double(length = n)
    tol = as.numeric(0.01)
    lwsp = as.integer(n * (m + 2) + 5 * (m + 2)^2 + ideg + 1)
    wsp = double(length = lwsp)
    liwsp = max(m + 2, 7)
    iwsp = integer(length = liwsp)
    if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE)) {
        tmatvec = t(matvec)
    }
    if ((exists("anorm") == FALSE) || is.null(anorm)) {
        if (transform_to_coo_TF == FALSE && transpose_needed == FALSE) {
            tmpQmat1 = coo2mat(matvec, n = coo_n)
            tmpQmat2 = t(tmpQmat1)
            anorm = as.numeric(norm(tmpQmat2, type = "O"))
        }
        else {
            anorm = as.numeric(norm(matvec, type = "O"))
        }
    }
    itrace = 0
    iflag = 0
    if (transform_to_coo_TF == TRUE) {
        tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
    }
    else {
        tmpmat_in_REXPOKIT_coo_fmt = matvec
    }
    ia = tmpmat_in_REXPOKIT_coo_fmt[, "ia"]
    ja = tmpmat_in_REXPOKIT_coo_fmt[, "ja"]
    a = tmpmat_in_REXPOKIT_coo_fmt[, "a"]
    num_tvals = length(tvals)
    if (is.null(inputprobs_for_fast)) {
        res = double(length = n * n)
        if ((num_tvals > 1) || (force_list_if_1_tval == TRUE)) {
            NA_matrix = matrix(NA, nrow = n, ncol = n)
            list_of_matrices_output = replicate(NA_matrix, n = num_tvals, simplify = FALSE)
            for (i in 1:num_tvals) {
                timeval = tvals[i]
                list_of_matrices_output[[i]] = .bgb_rexp_expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
            }
            return(list_of_matrices_output)
        }
        else {
            timeval = tvals
            output_Pmat = .bgb_rexp_expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
            return(output_Pmat)
        }
    }
    else {
        v = inputprobs_for_fast
        if ((num_tvals > 1) || (force_list_if_1_tval == TRUE)) {
            list_of_outprobs_output = matrix(NA, nrow = num_tvals, ncol = n)
            for (i in 1:num_tvals) {
                t = tvals[i]
                list_of_outprobs_output[i, ] = .bgb_rexp_expokit_mydgexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
            }
        }
        else {
            w_output_probs <- .bgb_rexp_expokit_mydgexpv_wrapper(n = n, m = m, t = tvals, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
            list_of_outprobs_output[1, ] = w_output_probs
            return(w_output_probs)
        }
        return(list_of_matrices_output)
    }
}

expokit_wrapalldmexpv_tvals <- function (Qmat = NULL, tvals = c(2.1), inputprobs_for_fast = NULL, transpose_needed = TRUE, transform_to_coo_TF = TRUE, coo_n = NULL, force_list_if_1_tval = FALSE, check_for_0_rows = TRUE) 
{
    defaults = "\n\tQmat=NULL\n\tt = 2.1\n\tinputprobs_for_fast=NULL\n\ttranspose_needed=TRUE\n\tCOO_needed=TRUE\n\ttransform_to_coo_TF=TRUE\n\tcoo_n=NULL\n\tforce_list_if_1_tval=FALSE\n\tcheck_for_0_rows=FALSE\n\tcheck_for_0_rows=1e-15\n\t"
    if (is.null(Qmat)) {
        warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
        Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow = 4, byrow = TRUE)
    }
    if (is.null(tvals)) {
        warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
        tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
    }
    matvec = Qmat
    if (is.null(inputprobs_for_fast)) {
        if (.bgb_precision() != 1L && check_for_0_rows != FALSE) {
            rows_w_all_zeros_TF = .bgb_rexp_findrows_w_all_zeros(matvec)
            if (all(rows_w_all_zeros_TF == FALSE)) {
                pass = 1
            }
            else {
                rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]
                if (is.numeric(check_for_0_rows)) {
                  check_for_0_rows = check_for_0_rows
                }
                else {
                  check_for_0_rows = 1e-15
                }
                newrowvals = rep(check_for_0_rows, ncol(matvec))
                matvec[rows_allzero_indices, ] = newrowvals
                diagonal_val = -1 * sum(newrowvals[-1])
                matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
                warning(paste(sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldmexpv_', ...)\nand therefore expokit_wrapalldmexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep = ""))
            }
        }
    }
    if (transform_to_coo_TF == TRUE) {
        nz = sum(matvec != 0)
        ia = integer(length = nz)
        ja = integer(length = nz)
        a = double(length = nz)
        n = nrow(matvec)
    }
    else {
        n = coo_n
        nz = sum(matvec[, "a"] != 0)
    }
    ideg = as.integer(6)
    m = n - 1
    if (is.null(inputprobs_for_fast)) {
        v = double(n)
        v[1] = 1
    }
    else {
        v = double(n)
        v = inputprobs_for_fast
    }
    tol = as.numeric(0.01)
    lwsp = as.integer(n * (m + 2) + 5 * (m + 2)^2 + ideg + 1)
    wsp = double(length = lwsp)
    liwsp = max(m + 2, 7)
    iwsp = integer(length = liwsp)
    if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE)) {
        tmatvec = t(matvec)
    }
    if ((exists("anorm") == FALSE) || is.null(anorm)) {
        if (transform_to_coo_TF == FALSE && transpose_needed == FALSE) {
            tmpQmat1 = coo2mat(matvec, n = coo_n)
            tmpQmat2 = t(tmpQmat1)
            anorm = as.numeric(norm(tmpQmat2, type = "O"))
        }
        else {
            anorm = as.numeric(norm(matvec, type = "O"))
        }
    }
    itrace = 0
    iflag = 0
    if (transform_to_coo_TF == TRUE) {
        tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
    }
    else {
        tmpmat_in_REXPOKIT_coo_fmt = matvec
    }
    ia = tmpmat_in_REXPOKIT_coo_fmt[, "ia"]
    ja = tmpmat_in_REXPOKIT_coo_fmt[, "ja"]
    a = tmpmat_in_REXPOKIT_coo_fmt[, "a"]
    num_tvals = length(tvals)
    if (is.null(inputprobs_for_fast)) {
        res = double(length = n * n)
        if ((num_tvals > 1) || (force_list_if_1_tval == TRUE)) {
            NA_matrix = matrix(NA, nrow = n, ncol = n)
            list_of_matrices_output = replicate(NA_matrix, n = num_tvals, simplify = FALSE)
            for (i in 1:num_tvals) {
                t = tvals[i]
                list_of_matrices_output[[i]] = .bgb_rexp_expokit_dmexpv_wrapper(n = n, m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
            }
            return(list_of_matrices_output)
        }
        else {
            t = tvals
            output_Pmat = .bgb_rexp_expokit_dmexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, ia = ia, ja = ja, a = a, nz = nz)
            return(output_Pmat)
        }
    }
    else {
        v = inputprobs_for_fast
        if ((num_tvals > 1) || (force_list_if_1_tval == TRUE)) {
            list_of_outprobs_output = matrix(NA, nrow = num_tvals, ncol = n)
            for (i in 1:num_tvals) {
                t = tvals[i]
                w_output_probs <- .bgb_rexp_expokit_mydmexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, itrace = itrace, iflag = iflag, ia = ia, ja = ja, a = a, nz = nz)
                list_of_outprobs_output[i, ] = w_output_probs
            }
        }
        else {
            t = tvals
            w_output_probs <- .bgb_rexp_expokit_mydmexpv_wrapper(n = n, m = m, t = t, v = v, tol = tol, anorm = anorm, wsp = wsp, lwsp = lwsp, iwsp = iwsp, liwsp = liwsp, itrace = itrace, iflag = iflag, ia = ia, ja = ja, a = a, nz = nz)
            list_of_outprobs_output[1, ] = w_output_probs
            return(w_output_probs)
        }
        return(list_of_matrices_output)
    }
}

mat2coo <- function (tmpmat) 
{
    defaults = "\n\ttmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)\n\t"
    numrows = nrow(tmpmat)
    numcols = ncol(tmpmat)
    if (numrows != numcols) {
        stop("ERROR! mat2coo(tmpmat) says that in tmpmat, nrows != ncols")
        return(NA)
    }
    numcells = numrows^2
    x = 1:numrows
    y = 1:numrows
    tmpmat_in_REXPOKIT_coo_fmt <- .bgb_rexp_rexpokit_as_coo(tmpmat)
    return(tmpmat_in_REXPOKIT_coo_fmt)
}

mat2coo_forloop <- function (tmpmat) 
{
    nz = sum(tmpmat != 0)
    ia = integer(length = nz)
    ja = integer(length = nz)
    a = integer(length = nz)
    count = 0
    for (i in 1:nrow(tmpmat)) {
        for (j in 1:ncol(tmpmat)) {
            if (tmpmat[i, j] != 0) {
                count = count + 1
                ia[count] = i
                ja[count] = j
                a[count] = tmpmat[i, j]
            }
        }
    }
    tmpmat_in_REXPOKIT_coo_fmt = cbind(ia, ja, a)
    return(tmpmat_in_REXPOKIT_coo_fmt)
}

SparseM_coo_to_REXPOKIT_coo <- function (tmpmat_in_SparseMcoo_fmt) 
{
    tmpcoo = tmpmat_in_SparseMcoo_fmt
    tmpmat_in_REXPOKIT_coo_fmt = cbind(tmpcoo@ia, tmpcoo@ja, tmpcoo@ra)
    colnames(tmpmat_in_REXPOKIT_coo_fmt) = c("ia", "ja", "a")
    return(tmpmat_in_REXPOKIT_coo_fmt)
}

.bgb_rexp_findrows_w_all_zeros <- function (matvec) 
{
    return(apply(X = matvec, MARGIN = 1, FUN = .bgb_rexp_row_allzero_TF))
}

.bgb_rexp_expokit_dmexpv_wrapper <- function (n, m, t, v, tol, anorm, wsp, lwsp, iwsp, liwsp, ia, ja, a, nz) 
{
    ret <- .bgb_rexp_call("R_dmexpv", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
    output_Pmat = matrix(ret$res, nrow = n, byrow = TRUE)
    return(output_Pmat)
}

.bgb_rexp_expokit_mydmexpv_wrapper <- function (n, m, t, v, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz) 
{
    ret <- .bgb_rexp_call("R_mydmexpv", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
    w_output_probs = matrix(ret, ncol = n, byrow = TRUE)
    return(w_output_probs)
}

.bgb_rexp_expokit_dgexpv_wrapper <- function (n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res) 
{
    ret <- .bgb_rexp_call('R_dgexpv',n,m,timeval,v,tol,anorm,wsp,lwsp,iwsp,liwsp,ia,ja,a,nz)
    tmpoutmat <- matrix(ret$res,nrow=n,byrow=TRUE)
    return(tmpoutmat)
}

.bgb_rexp_expokit_mydgexpv_wrapper <- function (n, m, t, v, tol, anorm, wsp, lwsp, iwsp, liwsp, ia, ja, a, nz) 
{
    res2 = NULL
    ret <- .bgb_rexp_call("R_mydgexpv", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
    w_output_probs = matrix(ret, ncol = n, byrow = TRUE)
    return(w_output_probs)
}

.bgb_rexp_rexpokit_as_coo <- function (x) 
{
    if (!is.double(x)) 
        storage.mode(x) <- "double"
    ret <- .bgb_rexp_call("R_rexpokit_as_coo", x)
    colnames(ret) <- c("ia", "ja", "a")
    return(ret)
}


.bgb_rexp_row_allzero_TF <- function(tmprow) all(tmprow == 0)
