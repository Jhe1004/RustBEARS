\name{native_runBSM}
\alias{runBSM}
\title{Legacy or prepared native whole-tree stochastic mapping}
\usage{
runBSM(res, stochastic_mapping_inputs_list, maxnum_maps_to_try = 1,
  nummaps_goal = maxnum_maps_to_try, maxtries_per_branch = 40000,
  save_after_every_try = TRUE, savedir = getwd(), seedval = as.numeric(Sys.time()),
  wait_before_save = 0.1, master_nodenum_toPrint = 0, ..., bsm_mode = "legacy",
  num_threads = NULL, bsm_engine = c("legacy_r_bridge", "native_whole_tree"),
  target_rmcse = NULL, min_maps = 200L, batch_size = 200L, max_maps = 10000L,
  min_mean_count = 1)
}
\description{Preserve the existing R bridge or draw complete histories from an
owned prepared Rust sampler, optionally stopping at a mapping precision target.}
\arguments{
\item{res}{Fitted result.}
\item{stochastic_mapping_inputs_list}{Legacy mapping inputs; pass NULL for native mode.}
\item{maxnum_maps_to_try, nummaps_goal}{Legacy attempt budget and requested fixed map count.}
\item{maxtries_per_branch}{Legacy branch rejection limit.}
\item{save_after_every_try, savedir, wait_before_save, master_nodenum_toPrint}{Legacy
save and printing settings; native mode returns full histories in memory.}
\item{seedval}{Master seed. Native mode floors the clock-based default; explicit
seeds must satisfy the exact integer requirements of \code{draw_prepared_bsm}.}
\item{...}{Additional legacy bridge arguments.}
\item{bsm_mode}{Legacy bridge mode or \code{"rejection"} for native rejection bridges.}
\item{num_threads}{Native thread count, defaulting to the RustBEARS.threads option.}
\item{bsm_engine}{Explicit legacy-compatible R bridge or prepared native whole-tree engine.}
\item{target_rmcse}{NULL for fixed count; otherwise an approximate relative precision target.}
\item{min_maps, batch_size, max_maps}{Initial assessment count, incremental batch size and hard cap.}
\item{min_mean_count}{Statistics below this mean are marked not assessed.}
}
\details{Native batches retain global map indices. Changing batch size and
thread count preserves maps at the same global indices. Adaptive stopping
reports target_precision_met or max_maps_reached independently of route ranking.
Ordinary intervals evaluated at this random stopping time do not have a proven
95 percent coverage guarantee; inspect the empirical calibration record.}
\value{Legacy bridge output, or a BSMResult with full event tables, provenance,
and adaptive diagnostics when requested.}
\seealso{\code{bsm_diagnostics}, \code{prepare_tree_sampler}}
