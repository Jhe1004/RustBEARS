# GPL >= 2. Legacy script bindings identified by R codetools.
if (getRversion() >= "2.15.1")
    utils::globalVariables(c(
        ".PlotPhyloEnv", "BSM_results", "COO_weights_columnar", "Pmat",
        "SSEsim_results_processed", "TF1", "TF2", "TF3",
        "TF4", "TF5", "TF6", "anc_statenums",
        "ancstate_accuracies_DECJ_subset", "ancstate_accuracies_DEC_subset", "areas", "branching.times",
        "edgelabels", "extract_lines_startstr_to_endstr", "fossils_older_than", "ia",
        "include_null_range", "inputs", "item", "ja",
        "list_of_matrices_w_variable_names", "newick_fns", "newvals", "newvals_st",
        "node_depth", "node_depth_edgelength", "node_height", "node_height_clado",
        "offending_entries", "on_NaN_error", "optim_results_mean", "optim_val",
        "printlevel", "res", "resDEC", "root",
        "siminf_stats_subset", "states_list_0based", "tr", "uniq_params_num",
        "unroot", "write_lines_good"
    ))
