# Contributing to RustBEARS

Thank you for your interest in contributing to **`RustBEARS`**.

## Scope and Upstream Attribution

- `RustBEARS` maintains a memory-safe, unified Rust computational engine (`rust/`), a native Python package (`rustbears/`), and a drop-in function-compatible R bridge (`r_package/RustBEARS/`) for the [`BioGeoBEARS`](http://phylo.wikidot.com/biogeobears) historical biogeography framework created by Nicholas J. Matzke.
- Upstream copyright headers, license notices, and scholarly co-citation metadata must never be removed.

## Development and Verification Checklist

1. **Offline Rust Core**: All Rust dependencies are vendored under `rust/vendor/`. Do not introduce network-dependent build steps or unvendored crates without updating `THIRD_PARTY_LICENSES.md`.
2. **Minimum Toolchain**: Keep Rust code compatible with `rustc >= 1.82`, Python `>= 3.10`, and R `>= 4.0.0`.
3. **Full Verification**: Before opening a pull request, run:
   ```bash
   ./scripts/build_release_artifacts.sh
   ```
   and confirm that all Rust unit tests, R static-library builds (`--features r-blas`), and Python regression tests pass cleanly.
4. **Regression Tests**: Bug fixes affecting numerical kernels, FFI memory ownership, or file parsers must include a corresponding unit test under `rust/` or `tests/`.
