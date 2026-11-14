#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

SKIP_PY_TESTS=0
for arg in "$@"; do
  case "$arg" in
    --skip-python-tests) SKIP_PY_TESTS=1 ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

echo "[1/5] Synchronizing unified Rust workspace and license notices ..."
rm -rf r_package/RustBEARS/src/rust
rsync -a --exclude='target' --exclude='.DS_Store' rust/ r_package/RustBEARS/src/rust/
# R builds a static archive linked against R's own BLAS/LAPACK symbols
python3 - <<'PY'
from pathlib import Path
p = Path("r_package/RustBEARS/src/rust/crates/bgb-ffi/Cargo.toml")
p.write_text(p.read_text().replace('crate-type = ["staticlib", "cdylib"]', 'crate-type = ["staticlib"]'))
PY
cp -f NOTICE-Rust-port.md r_package/RustBEARS/inst/NOTICE-Rust-port.md
cp -f NOTICE-Rust-port.md rustbears/NOTICE-Rust-port.md
cp -f THIRD_PARTY_LICENSES.md r_package/RustBEARS/inst/THIRD_PARTY_LICENSES.md
cp -f THIRD_PARTY_LICENSES.md rustbears/THIRD_PARTY_LICENSES.md

echo "[2/5] Building & testing Rust core (offline vendored) ..."
(
  cd rust
  cargo --config 'source.crates-io.replace-with="vendored-sources"' \
        --config 'source.vendored-sources.directory="vendor"' \
        test --locked --offline
  cargo --config 'source.crates-io.replace-with="vendored-sources"' \
        --config 'source.vendored-sources.directory="vendor"' \
        build --release --locked --offline -p bgb-ffi
)

echo "[3/5] Verifying R-package staticlib build (--features r-blas) ..."
(
  cd r_package/RustBEARS/src/rust
  cargo --config 'source.crates-io.replace-with="vendored-sources"' \
        --config 'source.vendored-sources.directory="vendor"' \
        build --release --locked --offline -p bgb-ffi --features r-blas
  rm -rf target
)

if [ "$SKIP_PY_TESTS" -eq 1 ]; then
  echo "[4/5] Skipping Python regression test suite (--skip-python-tests; covered by Python CI job) ..."
else
  echo "[4/5] Running self-contained Python regression test suite ..."
  python3 -m unittest discover -s tests -p "test_*.py" -v
fi

echo "[5/5] Packaging reproducible anonymized R source tarball into dist/RustBEARS_1.1.4.tar.gz ..."
python3 - <<'PY'
import gzip
import os
from pathlib import Path
import stat
import tarfile

epoch = int(os.environ.get("SOURCE_DATE_EPOCH", "1700000000"))
dist_dir = Path("dist")
dist_dir.mkdir(exist_ok=True)
out_tar = dist_dir / "RustBEARS_1.1.4.tar.gz"
pkg_dir = Path("r_package") / "RustBEARS"

def should_skip(rel_parts):
    if any(part in {".DS_Store", "__pycache__", ".Rapp.history"} for part in rel_parts):
        return True
    rel_posix = "/".join(rel_parts)
    if rel_posix.startswith("src/rust/target") or rel_posix.endswith((".o", ".so", ".dylib", ".dll")):
        return True
    return False

entries = sorted([pkg_dir, *pkg_dir.rglob("*")], key=lambda p: p.relative_to(pkg_dir).as_posix())

with open(out_tar, "wb") as raw_fh:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as gz_fh:
        with tarfile.open(mode="w", fileobj=gz_fh, format=tarfile.PAX_FORMAT) as tf:
            for path in entries:
                rel = path.relative_to(pkg_dir)
                if rel.parts and should_skip(rel.parts):
                    continue
                arcname = "RustBEARS" if not rel.parts else f"RustBEARS/{rel.as_posix()}"
                ti = tf.gettarinfo(str(path), arcname=arcname)
                ti.uid = 0
                ti.gid = 0
                ti.uname = "root"
                ti.gname = "root"
                ti.mtime = epoch
                ti.pax_headers = {}
                if ti.isdir():
                    ti.mode = 0o755
                    tf.addfile(ti)
                elif ti.isfile():
                    is_exec = bool(path.stat().st_mode & stat.S_IXUSR)
                    ti.mode = 0o755 if is_exec else 0o644
                    with open(path, "rb") as f:
                        tf.addfile(ti, f)
PY

echo "Clean release build & verification succeeded: dist/RustBEARS_1.1.4.tar.gz"
