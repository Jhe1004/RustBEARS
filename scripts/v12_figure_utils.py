"""Small shared utilities for source-bound manuscript figures."""
import csv
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def export_rows(path, rows):
    if not rows:
        raise ValueError(f"empty source table: {path}")
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def copy_evidence(paths, output):
    """Retain exact source files, without depending on the original host path."""
    folder = Path(output) / 'source_files'
    folder.mkdir(exist_ok=True)
    inventory = []
    for index, path in enumerate(paths):
        path = Path(path)
        destination = folder / f'{index:03d}_{path.name}'
        shutil.copyfile(path, destination)
        inventory.append(dict(original=str(path.resolve()), archived=destination.relative_to(output).as_posix(), sha256=sha(path)))
    write_json(Path(output) / 'source_inventory.json', inventory)
    return inventory


def configure_plotting():
    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7})


def save_figure(fig, output, stem):
    output = Path(output)
    fig.canvas.draw()
    for extension in ('pdf', 'svg', 'png', 'tiff'):
        options={'dpi':600} if extension in ('png','tiff') else {}
        if extension=='tiff':options['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(output / f'{stem}.{extension}', **options)
    write_json(output / 'render_provenance.json', dict(
        outputs={f'{stem}.{extension}': sha(output / f'{stem}.{extension}') for extension in ('pdf', 'svg', 'png', 'tiff')},
        status='rendered_awaiting_text_collision_alignment_and_visual_QA'))


def audit_alignment(fig, output, stem, audit_scripts):
    import sys
    sys.path.insert(0, str(audit_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    require_matplotlib_panel_alignment(fig, json_out=str(Path(output) / f'{stem}.alignment.json'),
        overlay_svg=str(Path(output) / f'{stem}.alignment.svg'), tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, strict=True)
    write_json(Path(output) / 'alignment_toolchain.json', dict(
        helper=str(Path(audit_scripts).resolve() / 'audit_panel_alignment.py'),
        sha256=sha(Path(audit_scripts) / 'audit_panel_alignment.py'), tolerance_pt=1.5))
