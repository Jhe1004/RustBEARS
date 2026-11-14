"""RustBEARS: Native Rust computational engine and Python workflows for BioGeoBEARS models."""
from ._ffi import CompiledContext, context_stats
from .io import PhyloTree, TipRanges, read_geog, read_newick, read_times
from .model import BioGeoBEARS, BioGeoBEARSModelObject, EvaluationResult, FitResult, ModelConfig, Parameter, model_template

RustBEARS = BioGeoBEARS
from .operators import (
    Rcpp_combn_zerostart, coo2crs, expokit_dgexpv, expokit_dgexpv_Qmat, expokit_dgpadm_Qmat2, expokit_dmexpv_Qmat,
    maxent, rcpp_areas_list_to_states_list, rcpp_calc_anclikes_sp, rcpp_calc_anclikes_sp_COOprobs,
    rcpp_calc_anclikes_sp_COOweights_faster, rcpp_calc_anclikes_sp_rowsums, rcpp_calc_anclikes_sp_using_COOprobs,
    rcpp_calc_rowsums_for_COOweights_columnar, rcpp_calc_splitlikes_using_COOweights_columnar,
    rcpp_calc_splitlikes_using_cooWeights_columnar, rcpp_convolve, rcpp_mult2probvect, rcpp_states_list_to_DEmat,
)
from .stratified import section_the_tree
from .traits import TraitModel
from .uncertainty.contracts import RootPriorSpec, ContractMismatchError, UnsupportedModelContractError
from .uncertainty.simulate import SimulatedBatch, PreparedSimulation, simulate_datasets
from .uncertainty.parallel import ParallelConfig
from .uncertainty.bootstrap import BootstrapLRTResult, bootstrap_lrt
from .uncertainty.profile import ParameterInterval, ProfileCIResult, ProfileOptimizationError, profile_ci, profile_surface
from .uncertainty.checkpoint import CheckpointStore, CheckpointMismatchError
from .uncertainty.averaging import ModelAverageResult, model_average
from .uncertainty.diagnostics import BSMDiagnosticsResult, bsm_diagnostics, run_bsm_adaptive
from .uncertainty.multitree import MultiTreeFitResult, fit_over_trees, TreeFailureError

__version__ = "1.2.0.dev0"
from .bsm import BSMResult, count_ana_clado_events, get_inputs_for_stochastic_mapping, runBSM, simulate_source_areas_ana_clado
from .plotting import plot_bsm
from .stats import AkaikeWeights_on_summary_table, compare_models, lrttest


def define_BioGeoBEARS_run(tree, geog, model="DEC", **options):
    return BioGeoBEARS(tree, geog, model, **options)


def bears_optim_run(BioGeoBEARS_run_object, skip_optim=False, skip_optim_option="return_all", **options):
    model = BioGeoBEARS_run_object
    if skip_optim:
        out = model.evaluate(calc_ancprobs=skip_optim_option=="return_all", **options)
        return out if skip_optim_option=="return_all" else out.lnL
    return model.fit(**options)


def calc_loglike_sp(model, params=None, *, return_what="loglike", **options):
    out = model.evaluate(params,calc_ancprobs=return_what!="loglike",**options)
    return out.lnL if return_what=="loglike" else out


def calc_loglike_sp_stratified(model, params=None, **options):
    if model.sections is None:
        raise ValueError("expected a stratified BioGeoBEARS model")
    return calc_loglike_sp(model,params,**options)
