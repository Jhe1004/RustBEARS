#!/usr/bin/env python3
"""Synchronize the R-package Rust source mirror and verify every source byte."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'rust'
DESTINATION = ROOT / 'r_package/RustBEARS/src/rust'


def source_files(root):
    return {p.relative_to(root): p for p in root.rglob('*')
            if p.is_file() and not any(x in {'target', '.DS_Store', '__pycache__'} for x in p.relative_to(root).parts)}


def expected_bytes(relative, path):
    data = path.read_bytes()
    if relative.as_posix() == 'crates/bgb-ffi/Cargo.toml':
        data = data.replace(b'crate-type = ["staticlib", "cdylib"]', b'crate-type = ["staticlib"]')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--manifest', type=Path)
    args = parser.parse_args()
    source = source_files(SOURCE)
    if not args.check:
        DESTINATION.mkdir(parents=True, exist_ok=True)
        for relative in source_files(DESTINATION).keys() - source.keys():
            (DESTINATION / relative).unlink()
        for relative, path in source.items():
            target = DESTINATION / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            data = expected_bytes(relative, path)
            if not target.exists() or target.read_bytes() != data:
                target.write_bytes(data)
                shutil.copymode(path, target)
    target = source_files(DESTINATION)
    if target.keys() != source.keys():
        raise SystemExit('Rust mirror file inventory differs')
    manifest = {}
    for relative, path in source.items():
        expected = expected_bytes(relative, path)
        if target[relative].read_bytes() != expected:
            raise SystemExit(f'Rust mirror differs: {relative}')
        manifest[relative.as_posix()] = hashlib.sha256(expected).hexdigest()
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    print(f'Verified {len(manifest)} Rust mirror files (R staticlib manifest exception applied).')


if __name__ == '__main__':
    main()
