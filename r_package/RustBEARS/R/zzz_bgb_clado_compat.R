# GPL >= 2. Original R compatibility wrappers by Nicholas J. Matzke.
# Extracted from cladoRcpp 0.15.1; native calls redirect to Rust.
areas_list_to_states_list_old <- function (areas = c("A", "B", "C"), maxareas = length(areas), include_null_range = TRUE, split_ABC = TRUE) 
{
    if (maxareas > length(areas)) {
        maxareas = length(areas)
    }
    nstates = numstates_from_numareas(numareas = length(areas), maxareas = maxareas, include_null_range = include_null_range)
    states_list = rep(NA, times = nstates)
    if (include_null_range == TRUE) {
        s = 1
    }
    else {
        s = 0
    }
    if (split_ABC == FALSE) {
        for (m in 1:maxareas) {
            states_matrix = combn(x = areas, m = m)
            tmp_states = apply(X = states_matrix, 2, paste, collapse = "")
            for (i in 1:length(tmp_states)) {
                states_list[[(s = s + 1)]] = tmp_states[[i]]
            }
        }
    }
    if (split_ABC == TRUE) {
        for (m in 1:maxareas) {
            states_matrix = combn(x = areas, m = m)
            tmp_states = apply(X = states_matrix, 2, paste, collapse = ",")
            for (i in 1:length(tmp_states)) {
                states_list[[(s = s + 1)]] = tmp_states[[i]]
            }
        }
        states_list = mapply(FUN = strsplit3, split = ",", states_list)
        names(states_list) = NULL
    }
    if (include_null_range == TRUE) {
        s = 1
        states_list[[s]] = c("_")
    }
    else {
        s = 1
    }
    return(states_list)
}

numstates_from_numareas <- function (numareas = 3, maxareas = numareas, include_null_range = FALSE) 
{
    nstates = 0
    for (i in 1:maxareas) {
        tmp_nstates = choose(n = numareas, k = i)
        nstates = nstates + tmp_nstates
    }
    if (include_null_range == TRUE) {
        nstates = nstates + 1
    }
    return(nstates)
}

rcpp_areas_list_to_states_list <- function (areas = c("A", "B", "C"), maxareas = length(areas), include_null_range = TRUE) 
{
    setup = "\n\tareas=c(\"A\",\"B\",\"C\")\n\trcpp_areas_list_to_states_list(areas=c(\"A\",\"B\",\"C\"), maxareas=length(areas), include_null_range=TRUE)\n\t\n\t"
    if (maxareas > length(areas)) {
        maxareas = length(areas)
    }
    nstates = numstates_from_numareas(numareas = length(areas), maxareas = maxareas, include_null_range = include_null_range)
    R_areas_indices = seq(0, length(areas) - 1, 1)
    R_states_list = as.list(rep(-1, times = nstates))
    R_maxareas = maxareas
    R_include_null_range = include_null_range
    R_states_list = .bgb_clado_call("cpp_areas_list_to_states_list", as.integer(R_areas_indices), as.integer(R_maxareas), as.logical(R_include_null_range))
    if (include_null_range == TRUE) {
        if (R_states_list[[1]] == -1) {
            R_states_list[[1]] = NA
        }
    }
    return(R_states_list)
}

rcpp_calc_anclikes_sp <- function (Rcpp_leftprobs, Rcpp_rightprobs, l, s = 1, v = 1, j = 0, y = 1, dmat = NULL, maxent01s = NULL, maxent01v = NULL, maxent01j = NULL, maxent01y = NULL, max_minsize_as_function_of_ancsize = NULL, Rsp_rowsums = rep(1, length(Rcpp_leftprobs)), printmat = FALSE) 
{
    if (all(c(length(Rcpp_leftprobs) == length(Rcpp_rightprobs), length(Rcpp_leftprobs) == length(l), length(l) == length(Rcpp_rightprobs))) == FALSE) {
        error_msg = paste("ERROR: These are not the same length (Rcpp_leftprobs, Rcpp_rightprobs, l): ", length(Rcpp_leftprobs), ", ", length(Rcpp_rightprobs), ", ", length(l), sep = "")
        return(stop(error_msg))
    }
    numareas = max(unlist(l), na.rm = TRUE) + 1
    if (is.null(max_minsize_as_function_of_ancsize)) {
        max_minsize_as_function_of_ancsize = rep(1, numareas)
    }
    if (is.null(dmat)) {
        dmat = matrix(data = 1, nrow = numareas, ncol = numareas)
    }
    if (is.null(maxent01s)) {
        maxent01s = matrix(0, nrow = numareas, ncol = numareas)
        maxent01s[, 1] = 1
    }
    if (is.null(maxent01v)) {
        numareas_smaller_vicariant_range = floor(numareas/2)
        maxent01v = matrix(0, nrow = numareas, ncol = numareas)
        maxent01v[, 1] = 1
    }
    if (is.null(maxent01j)) {
        maxent01j = matrix(0, nrow = numareas, ncol = numareas)
        maxent01j[, 1] = 1
    }
    if (is.null(maxent01y)) {
        maxent01y = matrix(0, nrow = numareas, ncol = numareas)
        maxent01y[, 1] = 1
    }
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    prob_ancestral_states = .bgb_clado_call("cpp_calc_anclikes_sp", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize), Rsp_rowsums = as.numeric(Rsp_rowsums))
    return(prob_ancestral_states)
}

rcpp_calc_anclikes_sp_COOprobs <- function (Rcpp_leftprobs, Rcpp_rightprobs, l, s = 1, v = 1, j = 0, y = 1, dmat = NULL, maxent01s = NULL, maxent01v = NULL, maxent01j = NULL, maxent01y = NULL, max_minsize_as_function_of_ancsize = NULL, printmat = TRUE) 
{
    numareas = max(unlist(l), na.rm = TRUE) + 1
    if (is.null(max_minsize_as_function_of_ancsize)) {
        max_minsize_as_function_of_ancsize = rep(1, numareas)
    }
    if (is.null(dmat)) {
        dmat = matrix(data = 1, nrow = numareas, ncol = numareas)
    }
    if (is.null(maxent01s)) {
        maxent01s = matrix(0, nrow = numareas, ncol = numareas)
        maxent01s[, 1] = 1
    }
    if (is.null(maxent01v)) {
        numareas_smaller_vicariant_range = floor(numareas/2)
        maxent01v = matrix(0, nrow = numareas, ncol = numareas)
        maxent01v[, 1] = 1
    }
    if (is.null(maxent01j)) {
        maxent01j = matrix(0, nrow = numareas, ncol = numareas)
        maxent01j[, 1] = 1
    }
    if (is.null(maxent01y)) {
        maxent01y = matrix(0, nrow = numareas, ncol = numareas)
        maxent01y[, 1] = 1
    }
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    list_weights_of_transitions = .bgb_clado_call("cpp_calc_anclikes_sp_COOprobs", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize))
    return(list_weights_of_transitions)
}

rcpp_calc_anclikes_sp_COOweights_faster <- function (Rcpp_leftprobs, Rcpp_rightprobs, l, s = 1, v = 1, j = 0, y = 1, dmat = NULL, maxent01s = NULL, maxent01v = NULL, maxent01j = NULL, maxent01y = NULL, max_minsize_as_function_of_ancsize = NULL, printmat = TRUE, m = NULL, m_null_range = TRUE, jts_matrix = NULL) 
{
    numareas = max(unlist(l), na.rm = TRUE) + 1
    if (is.null(max_minsize_as_function_of_ancsize)) {
        max_minsize_as_function_of_ancsize = rep(1, numareas)
    }
    if (is.null(dmat)) {
        dmat = matrix(data = 1, nrow = numareas, ncol = numareas)
    }
    if (is.null(maxent01s)) {
        maxent01s = matrix(0, nrow = numareas, ncol = numareas)
        maxent01s[, 1] = 1
    }
    if (is.null(maxent01v)) {
        numareas_smaller_vicariant_range = floor(numareas/2)
        maxent01v = matrix(0, nrow = numareas, ncol = numareas)
        maxent01v[, 1] = 1
    }
    if (is.null(maxent01j)) {
        maxent01j = matrix(0, nrow = numareas, ncol = numareas)
        maxent01j[, 1] = 1
    }
    if (is.null(maxent01y)) {
        maxent01y = matrix(0, nrow = numareas, ncol = numareas)
        maxent01y[, 1] = 1
    }
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    if (is.null(m) == TRUE) {
        COO_weights_columnar = .bgb_clado_call("cpp_calc_anclikes_sp_COOweights_faster", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize))
        return(COO_weights_columnar)
    }
    if (is.null(m) == FALSE) {
        COO_weights_columnar_1 = NULL
        COO_weights_columnar_2 = NULL
        COO_weights_columnar_3 = NULL
        COO_weights_columnar_4 = NULL
        if (m_null_range == TRUE) {
            numstates_ranges_times_ms = length(Rcpp_leftprobs) + 1
            numstates_ranges = numstates_ranges_times_ms/length(m)
            numstates_ranges_in_default_COOmat = numstates_ranges - 1
        }
        else {
            numstates_ranges_times_ms = length(Rcpp_leftprobs) + 0
            numstates_ranges = numstates_ranges_times_ms/length(m)
            numstates_ranges_in_default_COOmat = numstates_ranges - 0
        }
        dmat_orig = dmat
        start_index = 1
        end_index = numstates_ranges_in_default_COOmat
        j_events_TF = NULL
        trait_states = NULL
        for (i in 1:length(m)) {
            dmat_new = m[i] * dmat
            if (i > 1) {
                start_index = end_index + 1
                end_index = start_index + numstates_ranges - 1
            }
            tmpa = Rcpp_leftprobs[start_index:end_index]
            tmpb = Rcpp_rightprobs[start_index:end_index]
            tmp_COO_weights_columnar = .bgb_clado_call("cpp_calc_anclikes_sp_COOweights_faster", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat_new), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize))
            if (is.null(jts_matrix) == FALSE) {
                num_clado_events = length(tmp_COO_weights_columnar[[1]])
                tmp_j_events_TF = rep(FALSE, num_clado_events)
                tmp_trait_states = rep(i, num_clado_events)
                states_list_area_0based = l
                for (ci in 1:num_clado_events) {
                  ancstate = tmp_COO_weights_columnar[[1]][ci] + 1
                  Lstate = tmp_COO_weights_columnar[[2]][ci] + 1
                  Rstate = tmp_COO_weights_columnar[[3]][ci] + 1
                  anc_areas = c(l[[ancstate]])
                  L_areas = c(l[[Lstate]])
                  R_areas = c(l[[Rstate]])
                  desc_areas = sort(unique(c(L_areas, R_areas)))
                  desc_in_anc_TF = desc_areas %in% anc_areas
                  tmp_j_events_TF[ci] = all(desc_in_anc_TF) == FALSE
                }
                j_events_TF = c(j_events_TF, tmp_j_events_TF)
                trait_states = c(trait_states, tmp_trait_states)
            }
            if (m_null_range == FALSE) {
                tmp_COO_weights_columnar[[1]] = tmp_COO_weights_columnar[[1]] + (i - 1) * numstates_ranges
                tmp_COO_weights_columnar[[2]] = tmp_COO_weights_columnar[[2]] + (i - 1) * numstates_ranges
                tmp_COO_weights_columnar[[3]] = tmp_COO_weights_columnar[[3]] + (i - 1) * numstates_ranges
            }
            else {
                tmp_COO_weights_columnar[[1]] = tmp_COO_weights_columnar[[1]] + (i - 1) * numstates_ranges
                tmp_COO_weights_columnar[[2]] = tmp_COO_weights_columnar[[2]] + (i - 1) * numstates_ranges
                tmp_COO_weights_columnar[[3]] = tmp_COO_weights_columnar[[3]] + (i - 1) * numstates_ranges
            }
            COO_weights_columnar_1 = c(COO_weights_columnar_1, tmp_COO_weights_columnar[[1]])
            COO_weights_columnar_2 = c(COO_weights_columnar_2, tmp_COO_weights_columnar[[2]])
            COO_weights_columnar_3 = c(COO_weights_columnar_3, tmp_COO_weights_columnar[[3]])
            COO_weights_columnar_4 = c(COO_weights_columnar_4, tmp_COO_weights_columnar[[4]])
        }
        COO_weights_columnar = NULL
        COO_weights_columnar[[1]] = COO_weights_columnar_1
        COO_weights_columnar[[2]] = COO_weights_columnar_2
        COO_weights_columnar[[3]] = COO_weights_columnar_3
        COO_weights_columnar[[4]] = COO_weights_columnar_4
        new_COO_weights_columnar = COO_weights_columnar
        addl_COO_weights_columnar_1 = NULL
        addl_COO_weights_columnar_2 = NULL
        addl_COO_weights_columnar_3 = NULL
        addl_COO_weights_columnar_4 = NULL
        aind = 0
        if (is.null(jts_matrix) == FALSE) {
            for (jt_i in 1:length(m)) {
                for (jt_j in 1:length(m)) {
                  jts_fract = jts_matrix[jt_i, jt_j]
                  if (jt_i == jt_j) {
                    TF2 = trait_states == jt_i
                    TF = (j_events_TF + TF2) == 2
                    new_COO_weights_columnar[[4]][TF] = jts_fract * COO_weights_columnar[[4]][TF]
                    (next)()
                  }
                  if (jts_fract > 0) {
                    TF2 = trait_states == jt_i
                    starttrait_TF = (j_events_TF + TF2) == 2
                    TF3 = trait_states == jt_j
                    endtrait_TF = (j_events_TF + TF3) == 2
                    if (sum(starttrait_TF) == 0) {
                      (next)()
                    }
                    if (sum(endtrait_TF) == 0) {
                      (next)()
                    }
                    for (TF_in in 1:sum(starttrait_TF)) {
                      addl_COO_weights_columnar_1[(aind = aind + 1)] = COO_weights_columnar[[1]][starttrait_TF][TF_in]
                      addl_COO_weights_columnar_2[(aind)] = COO_weights_columnar[[2]][endtrait_TF][TF_in]
                      addl_COO_weights_columnar_3[aind] = COO_weights_columnar[[3]][endtrait_TF][TF_in]
                      addl_COO_weights_columnar_4[aind] = COO_weights_columnar[[4]][starttrait_TF][TF_in] * jts_fract
                    }
                  }
                }
            }
            COO_weights_columnar[[1]] = c(new_COO_weights_columnar[[1]], addl_COO_weights_columnar_1)
            COO_weights_columnar[[2]] = c(new_COO_weights_columnar[[2]], addl_COO_weights_columnar_2)
            COO_weights_columnar[[3]] = c(new_COO_weights_columnar[[3]], addl_COO_weights_columnar_3)
            COO_weights_columnar[[4]] = c(new_COO_weights_columnar[[4]], addl_COO_weights_columnar_4)
        }
        return(COO_weights_columnar)
    }
}

rcpp_calc_anclikes_sp_prebyte <- function (Rcpp_leftprobs, Rcpp_rightprobs, l, s = 1, v = 1, j = 0, y = 1, dmat = NULL, maxent01s = NULL, maxent01v = NULL, maxent01j = NULL, maxent01y = NULL, max_minsize_as_function_of_ancsize = NULL, Rsp_rowsums = rep(1, length(Rcpp_leftprobs)), printmat = FALSE) 
{
    if (all(c(length(Rcpp_leftprobs) == length(Rcpp_rightprobs), length(Rcpp_leftprobs) == length(l), length(l) == length(Rcpp_rightprobs))) == FALSE) {
        error_msg = paste("ERROR: These are not the same length (Rcpp_leftprobs, Rcpp_rightprobs, l): ", length(Rcpp_leftprobs), ", ", length(Rcpp_rightprobs), ", ", length(l), sep = "")
        return(stop(error_msg))
    }
    numareas = max(unlist(l), na.rm = TRUE) + 1
    if (is.null(max_minsize_as_function_of_ancsize)) {
        max_minsize_as_function_of_ancsize = rep(1, numareas)
    }
    if (is.null(dmat)) {
        dmat = matrix(data = 1, nrow = numareas, ncol = numareas)
    }
    if (is.null(maxent01s)) {
        maxent01s = matrix(0, nrow = numareas, ncol = numareas)
        maxent01s[, 1] = 1
    }
    if (is.null(maxent01v)) {
        numareas_smaller_vicariant_range = floor(numareas/2)
        maxent01v = matrix(0, nrow = numareas, ncol = numareas)
        maxent01v[, 1] = 1
    }
    if (is.null(maxent01j)) {
        maxent01j = matrix(0, nrow = numareas, ncol = numareas)
        maxent01j[, 1] = 1
    }
    if (is.null(maxent01y)) {
        maxent01y = matrix(0, nrow = numareas, ncol = numareas)
        maxent01y[, 1] = 1
    }
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    prob_ancestral_states = .bgb_clado_call("cpp_calc_anclikes_sp", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize), Rsp_rowsums = as.numeric(Rsp_rowsums))
    return(prob_ancestral_states)
}

rcpp_calc_anclikes_sp_rowsums <- function (Rcpp_leftprobs, Rcpp_rightprobs, l, s = 1, v = 1, j = 0, y = 1, dmat = NULL, maxent01s = NULL, maxent01v = NULL, maxent01j = NULL, maxent01y = NULL, max_minsize_as_function_of_ancsize = NULL, printmat = TRUE) 
{
    numareas = max(unlist(l), na.rm = TRUE) + 1
    if (is.null(max_minsize_as_function_of_ancsize)) {
        max_minsize_as_function_of_ancsize = rep(1, numareas)
    }
    if (is.null(dmat)) {
        dmat = matrix(data = 1, nrow = numareas, ncol = numareas)
    }
    if (is.null(maxent01s)) {
        maxent01s = matrix(0, nrow = numareas, ncol = numareas)
        maxent01s[, 1] = 1
    }
    if (is.null(maxent01v)) {
        numareas_smaller_vicariant_range = floor(numareas/2)
        maxent01v = matrix(0, nrow = numareas, ncol = numareas)
        maxent01v[, 1] = 1
    }
    if (is.null(maxent01j)) {
        maxent01j = matrix(0, nrow = numareas, ncol = numareas)
        maxent01j[, 1] = 1
    }
    if (is.null(maxent01y)) {
        maxent01y = matrix(0, nrow = numareas, ncol = numareas)
        maxent01y[, 1] = 1
    }
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    sp_rowsums = .bgb_clado_call("cpp_calc_anclikes_sp_rowsums", Rprintmat = as.integer(Rprintmat), cpp_leftprobs = as.numeric(tmpa), cpp_leftprobs = as.numeric(tmpb), l = l, s = as.numeric(s), v = as.numeric(v), j = as.numeric(j), y = as.numeric(y), dmat = as.matrix(dmat), maxent01s = as.matrix(maxent01s), maxent01v = as.matrix(maxent01v), maxent01j = as.matrix(maxent01j), maxent01y = as.matrix(maxent01y), max_minsize_as_function_of_ancsize = as.integer(max_minsize_as_function_of_ancsize))
    return(sp_rowsums)
}

rcpp_calc_anclikes_sp_using_COOprobs <- function (Rcpp_leftprobs, Rcpp_rightprobs, RCOO_left_i_list, RCOO_right_j_list, RCOO_probs_list, Rsp_rowsums, printmat = TRUE) 
{
    tmpa = Rcpp_leftprobs
    tmpb = Rcpp_rightprobs
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    R_anc_relprobs = .bgb_clado_call("cpp_calc_anclikes_sp_using_COOprobs", Rprintmat = as.integer(Rprintmat), leftprobs = as.numeric(tmpa), rightprobs = as.numeric(tmpb), RCOO_left_i_list = as.list(RCOO_left_i_list), RCOO_right_j_list = as.list(RCOO_right_j_list), RCOO_probs_list = as.list(RCOO_probs_list), Rsp_rowsums = as.numeric(Rsp_rowsums))
    return(R_anc_relprobs)
}

rcpp_calc_rowsums_for_COOweights_columnar <- function (COO_weights_columnar, numstates = 1 + max(sapply(X = COO_weights_columnar, FUN = max)[1:3]), printmat = TRUE) 
{
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
    RCOO_probs_list = COO_weights_columnar[[4]]
    rowsums = .bgb_clado_call("cpp_calc_rowsums_for_COOweights_columnar", RCOO_weights_columnar_anc_i_list = as.integer(RCOO_weights_columnar_anc_i_list), RCOO_probs_list = as.numeric(RCOO_probs_list), Rnumstates = as.integer(numstates))
    return(rowsums)
}

rcpp_calc_splitlikes_using_COOweights_columnar <- function (Rcpp_leftprobs, Rcpp_rightprobs, COO_weights_columnar, Rsp_rowsums, printmat = TRUE) 
{
    if (printmat == TRUE) {
        Rprintmat = 1
    }
    else {
        Rprintmat = 0
    }
    RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
    RCOO_left_i_list = COO_weights_columnar[[2]]
    RCOO_right_j_list = COO_weights_columnar[[3]]
    RCOO_probs_list = COO_weights_columnar[[4]]
    splitlikes = .bgb_clado_call("cpp_calc_splitlikes_using_COOweights_columnar", leftprobs = as.numeric(Rcpp_leftprobs), rightprobs = as.numeric(Rcpp_rightprobs), RCOO_weights_columnar_anc_i_list = as.integer(RCOO_weights_columnar_anc_i_list), RCOO_left_i_list = as.integer(RCOO_left_i_list), RCOO_right_j_list = as.integer(RCOO_right_j_list), RCOO_probs_list = as.numeric(RCOO_probs_list), Rsp_rowsums = as.numeric(Rsp_rowsums))
    return(splitlikes)
}

Rcpp_combn_zerostart <- function (n_to_choose_from, k_to_choose, maxlim = 1e+07) 
{
    defaults = "\n\tn_to_choose_from = 5\n\tk_to_choose = 6\n\tmaxlim=1e+07\n\t"
    n = n_to_choose_from
    m = k_to_choose
    if (n < 1) {
        txt = cat("\nERROR: n_to_choose_from must be > 0.  Your value is ", n_to_choose_from, "\n", sep = "")
        cat(txt)
        return(stop(txt))
    }
    if (m < 1) {
        txt = cat("\nERROR: k_to_choose must be > 0.  Your value is ", k_to_choose, "\n", sep = "")
        cat(txt)
        return(stop(txt))
    }
    if (n < m) {
        txt = cat("\nERROR: k_to_choose must be less than or equal to n_to_choose_from.\n", "Your values are:\n", "n_to_choose_from =\t\t", n_to_choose_from, "\n", "k_to_choose =\t\t", k_to_choose, "\n", sep = "")
        cat(txt)
        return(stop(txt))
    }
    predicted_number_of_cells_to_fill = choose(n, m)
    if ((predicted_number_of_cells_to_fill > maxlim) || (predicted_number_of_cells_to_fill < 1)) {
        txt = paste("ERROR: n=", n_to_choose_from, ", k=", k_to_choose, ", n choose k=", predicted_number_of_cells_to_fill, " > maxlim=", maxlim, "\nCalculating something this big may crash your computer!; (Note: (n choose k)<1 also raises this error.", sep = "")
        cat(txt)
        return(stop(txt))
    }
    outarray = .bgb_clado_call("cpp_combn_zerostart", as.integer(n), as.integer(m), as.integer(maxlim))
    return(outarray)
}

rcpp_convolve <- function (a, b) 
{
    convolve_result_vector = .bgb_clado_call("convolve3cpp", a = a, b = b)
    return(convolve_result_vector)
}

rcpp_mult2probvect <- function (a, b) 
{
    tcross_product_vector = .bgb_clado_call("mult2probvect", a = a, b = b)
    return(tcross_product_vector)
}

rcpp_states_list_to_DEmat <- function (areas_list, states_list, dmat, elist, amat = NULL, include_null_range = TRUE, normalize_TF = TRUE, makeCOO_TF = FALSE, min_precision = 1e-26) 
{
    if (is.null(amat)) {
        amat = dmat
        amat[is.numeric(amat)] = 0
        amat
    }
    R_areas_indices = areas_list
    R_states_list = states_list
    R_include_null_range = include_null_range
    R_include_null_range[is.na(R_include_null_range)] = -1
    if (normalize_TF == TRUE) {
        R_normalize_TF = 1
    }
    else {
        R_normalize_TF = 0
    }
    if (makeCOO_TF == TRUE) {
        DEmat_COO = .bgb_clado_call("cpp_states_list_to_DEmat_COO", as.list(R_areas_indices), as.list(R_states_list), as.matrix(dmat), as.numeric(elist), as.matrix(amat), as.integer(R_normalize_TF), as.numeric(min_precision))
        return(as.data.frame(DEmat_COO))
    }
    else {
        DEmat = .bgb_clado_call("cpp_states_list_to_DEmat", as.list(R_areas_indices), as.list(R_states_list), as.matrix(dmat), as.numeric(elist), as.matrix(amat), as.integer(R_normalize_TF))
        return(DEmat)
    }
    return()
}

strsplit3 <- function (x, ...) 
{
    out = strsplit(x, ...)[[1]]
    return(out)
}

