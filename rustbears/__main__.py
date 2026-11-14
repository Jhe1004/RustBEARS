"""Run model fitting and optional BSM from a shell."""
import argparse
from pathlib import Path
import sys

from . import BioGeoBEARS, compare_models, plot_bsm


def _run_cli(argv=None):
    parser = argparse.ArgumentParser(prog="rustbears")
    parser.add_argument("--traceback", action="store_true", help="Show full Python traceback on error")
    cmd = parser.add_subparsers(dest="command", required=True)
    fit = cmd.add_parser("fit")
    fit.add_argument("--tree", required=True)
    fit.add_argument("--geog", required=True)
    fit.add_argument("--models", default="DEC")
    fit.add_argument("--max-range-size", type=int)
    fit.add_argument("--use-ambiguities", action="store_true", help="Expand '?' tip geography ambiguities")
    fit.add_argument("--precision", choices=["legacy_f32", "strict_f64"], default="strict_f64")
    fit.add_argument("--timeperiods")
    for name in ("distances", "envdistances", "areas-allowed", "areas-adjacency", "dispersal-multipliers"):
        fit.add_argument("--" + name)
    fit.add_argument("--out-dir", type=Path, required=True)
    fit.add_argument("--threads", type=int, default=1)
    fit.add_argument("--maps", type=int, default=0)
    fit.add_argument("--seed", type=int, default=12345)
    fit.add_argument("--plot", choices=["pie", "text"], default="pie")
    fit.add_argument(
        "--traceback",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Show full Python traceback on error",
    )
    args = parser.parse_args(argv)
    try:
        if args.maps < 0:
            raise ValueError("--maps must be nonnegative")
        model_names = [m.strip() for m in args.models.split(",") if m.strip()]
        if not model_names:
            raise ValueError("--models must name at least one model")
        args.out_dir.mkdir(parents=True, exist_ok=True)
        results = []
        opts = {
            k: getattr(args, k)
            for k in (
                "max_range_size",
                "use_ambiguities",
                "precision",
                "timeperiods",
                "distances",
                "envdistances",
                "areas_allowed",
                "areas_adjacency",
                "dispersal_multipliers",
            )
        }
        for name in model_names:
            model = BioGeoBEARS(args.tree, args.geog, name, **opts)
            result = model.fit(num_threads=args.threads)
            results.append(result)
            out = args.out_dir / model.model.name
            out.mkdir(exist_ok=True)
            result.evaluation.to_dataframe().to_csv(out / "ancestral.csv", index=False)
            result.plot(out / "ancestral.pdf", args.plot)
            if args.maps:
                bsm = result.run_bsm(args.maps, num_threads=args.threads, seed=args.seed)
                bsm.clado_events_df.to_csv(out / "clado.csv", index=False)
                bsm.ana_events_df.to_csv(out / "ana.csv", index=False)
                bsm.summary_counts["counts_list"].to_csv(out / "counts.csv")
                plot_bsm(bsm, out / "events.pdf")
        compare_models(results).to_csv(args.out_dir / "model-comparison.csv", index=False)
    except (OSError, ValueError, RuntimeError) as exc:
        if args.traceback:
            raise
        sys.stderr.write(f"rustbears: error: {exc}\n")
        raise SystemExit(1) from None


def main():
    _run_cli()


if __name__ == "__main__":
    main()
