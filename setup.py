"""Build a Python-version-independent, platform-native Rust library."""
import gzip
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tarfile

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist
from wheel.bdist_wheel import bdist_wheel


def _native_lib_filename() -> str:
    if sys.platform == "win32":
        return "bgb_ffi.dll"
    if sys.platform == "darwin":
        return "libbgb_ffi.dylib"
    return "libbgb_ffi.so"


class NativeDistribution(Distribution):
    def has_ext_modules(self):
        return True


class BuildWheel(bdist_wheel):
    def get_tag(self):
        return ("py3", "none", super().get_tag()[2])


class ReproducibleSdist(sdist):
    """Write deterministic, root-anonymized source distributions under SOURCE_DATE_EPOCH."""

    def make_distribution(self):
        base_dir = self.distribution.get_fullname()
        self.make_release_tree(base_dir, self.filelist.files)
        try:
            dist_dir = Path(self.dist_dir)
            dist_dir.mkdir(parents=True, exist_ok=True)
            out_tar = dist_dir / f"{base_dir}.tar.gz"
            epoch = int(os.environ.get("SOURCE_DATE_EPOCH", "1700000000"))
            root_path = Path(base_dir)
            entries = sorted(
                [root_path, *root_path.rglob("*")],
                key=lambda p: p.relative_to(root_path).as_posix(),
            )
            with open(out_tar, "wb") as raw_fh:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as gz_fh:
                    with tarfile.open(mode="w", fileobj=gz_fh, format=tarfile.PAX_FORMAT) as tf:
                        for path in entries:
                            rel = path.relative_to(root_path)
                            if any(part in {".DS_Store", "__pycache__"} for part in rel.parts):
                                continue
                            arcname = base_dir if not rel.parts else f"{base_dir}/{rel.as_posix()}"
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
            self.archive_files = [str(out_tar)]
            self.distribution.dist_files.append(("sdist", "", str(out_tar)))
        finally:
            if not self.keep_temp:
                shutil.rmtree(base_dir, ignore_errors=True)


class BuildRust(build_py):
    def run(self):
        super().run()
        root = Path(__file__).resolve().parent
        rust = root / "rust"
        target = Path(os.environ.get("CARGO_TARGET_DIR", rust / "target")).resolve()
        env = dict(os.environ, CARGO_TARGET_DIR=str(target))
        subprocess.run(
            [
                "cargo",
                "--config",
                'source.crates-io.replace-with="vendored-sources"',
                "--config",
                'source.vendored-sources.directory="vendor"',
                "build",
                "--locked",
                "--offline",
                "--release",
                "-p",
                "bgb-ffi",
            ],
            cwd=rust,
            env=env,
            check=True,
        )
        name = _native_lib_filename()
        pkg_dir = Path(self.build_lib) / "rustbears"
        pkg_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target / "release" / name, pkg_dir / name)
        for doc in ("LICENSE", "NOTICE-Rust-port.md", "THIRD_PARTY_LICENSES.md"):
            src_doc = root / doc
            if src_doc.exists():
                shutil.copy2(src_doc, pkg_dir / doc)


setup(
    distclass=NativeDistribution,
    cmdclass={"build_py": BuildRust, "bdist_wheel": BuildWheel, "sdist": ReproducibleSdist},
)
