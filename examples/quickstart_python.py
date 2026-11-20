"""Self-contained Python quickstart demonstrating survival conditioning, BSM, and uncertainty tools."""
import rustbears as rb


def main() -> None:
    tree = "((K_sp1:2.0,O_sp1:2.0):2.0,((O_sp2:1.0,M_sp1:1.0):1.5,(M_sp2:1.0,H_sp1:1.0):1.5):1.5);"
    geog = {
        "K_sp1": "1100", "O_sp1": "0100", "O_sp2": "0110",
        "M_sp1": "0010", "M_sp2": "0011", "H_sp1": "0001",
    }

    models = ["DEC", "DEC+J", "DIVALIKE", "BAYAREALIKE"]
    fits = [
        rb.BioGeoBEARS(
            tree=tree,
            geog=geog,
            model=m,
            max_range_size=4,
            likelihood_conditioning="tip_survival",
            root_prior=rb.RootPriorSpec("uniform_nonnull"),
        ).fit(num_threads=2)
        for m in models
    ]

    table = rb.compare_models(fits)
    print(table[["model", "lnL", "numparams", "d", "e", "j", "AIC"]])

    prof = fits[0].profile_ci(params=["d", "e"], n_grid=8, num_threads=2)
    print(prof.to_dataframe()[["param", "mle", "lower", "upper", "status_lower", "status_upper"]])

    avg = rb.model_average(fits, criterion="AIC")
    maps = avg.run_bsm(num_maps=100, num_threads=2, seed=42)
    diag = rb.bsm_diagnostics(maps)
    assert len(diag.event_summary_df) > 0
    print("PYTHON_QUICKSTART_PASS")


if __name__ == "__main__":
    main()
