# Maintenance and Supported Scope

`RustBEARS` is maintained by Jian He while preserving full attribution and scholarly co-citation of the foundational `BioGeoBEARS` framework by Nicholas J. Matzke (2013, 2014).

## Reporting Issues

Please open bug reports and installation questions on [GitHub Issues](https://github.com/Jhe1004/RustBEARS/issues) and include:
- `RustBEARS` version (or git commit SHA), operating system, and R / Python / `rustc` versions;
- Precision mode (`StrictF64` / `strict_f64` vs. `LegacyF32` / `legacy_f32`) and likelihood formulation (`tip_survival`, `unconditional_complete`, or `legacy_unconditioned`);
- A minimal reproducible script or tree/geography snippet.

## Supported Modelling Scope

- **Legacy compatibility mode (`legacy_unconditioned`)**: Supports the full historical `BioGeoBEARS` feature set, including fossil tips, non-monotone area masks, detection models, and trait-coupled workflows.
- **Survival-conditioned (`tip_survival`) & complete unconditional (`unconditional_complete`) workflows**: Support contemporary ultrametric trees across DEC, DIVALIKE, BAYAREALIKE, and `+J` model families, spatial/environmental/area modifiers (`+x`, `+n`, `+w`, `+u`, `areas_adjacency`), and multi-stratum time slicing (`timeperiods`) with chronologically monotone area emergence (`areas_allowed`). Non-monotone submergence masks, non-contemporaneous fossil tips, and trait-coupled models remain restricted to `legacy_unconditioned`.
