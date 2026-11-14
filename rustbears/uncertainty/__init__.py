"""Probability contracts and uncertainty utilities for RustBEARS."""
from .contracts import ContractMismatchError, RootPriorSpec, UnsupportedModelContractError
from .guards import PreparedFitComparison, prepare_comparison_fits, validate_compatible_fits
from .parallel import ParallelConfig
from .bootstrap import BootstrapLRTResult, bootstrap_lrt
from .profile import ParameterInterval, ProfileCIResult, ProfileOptimizationError, profile_ci, profile_surface
from .checkpoint import CheckpointStore, CheckpointMismatchError
from .averaging import ModelAverageResult, model_average
from .diagnostics import BSMDiagnosticsResult, bsm_diagnostics, run_bsm_adaptive


def __getattr__(name):
    if name in ("MultiTreeFitResult","fit_over_trees","TreeFailureError"):
        from . import multitree
        return getattr(multitree,name)
    raise AttributeError(name)

__all__ = ["RootPriorSpec", "ContractMismatchError", "UnsupportedModelContractError",
           "PreparedFitComparison", "prepare_comparison_fits", "validate_compatible_fits",
           "ParallelConfig", "BootstrapLRTResult", "bootstrap_lrt", "ParameterInterval", "ProfileCIResult",
           "ProfileOptimizationError", "profile_ci", "profile_surface", "CheckpointStore", "CheckpointMismatchError",
           "ModelAverageResult", "model_average", "BSMDiagnosticsResult", "bsm_diagnostics", "run_bsm_adaptive",
           "MultiTreeFitResult", "fit_over_trees", "TreeFailureError"]
