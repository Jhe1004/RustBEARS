#!/usr/bin/env python3
"""Run the full pytest suite against the installed wheel, outside the checkout."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
evidence = root / 'evidence'
evidence.mkdir(exist_ok=True)
with tempfile.TemporaryDirectory(prefix='rustbears-wheel-check-') as folder:
    work = Path(folder)
    shutil.copytree(root / 'tests', work / 'tests')
    shutil.copytree(root / 'scripts', work / 'scripts',
                    ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    subprocess.run([sys.executable, '-c',
        'import pathlib, importlib.metadata, rustbears; p=pathlib.Path(rustbears.__file__).resolve().parent; '
        'assert any(x in str(p) for x in ("site-packages","dist-packages")), p; '
        'assert (p/"LICENSE").is_file(); assert (p/"THIRD_PARTY_LICENSES.md").is_file(); '
        'assert rustbears.__version__ == importlib.metadata.version("rustbears")'],cwd=work,env=env,check=True)
    subprocess.run([sys.executable, '-m', 'pytest', '-q', 'tests',
        '--junitxml='+str(evidence / 'ci-python.xml')],cwd=work,env=env,check=True)
    subprocess.run([sys.executable, '-m', 'rustbears', '--help'],cwd=work,env=env,check=True)
