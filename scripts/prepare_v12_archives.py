#!/usr/bin/env python3
"""Prepare local, explicitly unpublished code and research-data archives.

The code archive comes from one clean Git commit. Data files retain their
relative paths and SHA-256 provenance. This command never deposits an archive,
registers a DOI, changes repository visibility, or marks a release gate passed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zipfile


def sha(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):value.update(chunk)
    return value.hexdigest()


def run_git(repository,*args):
    return subprocess.check_output(['git','-C',str(repository),*args],text=True).strip()


def research_roots(study,extra_sources,output):
    """Explicit additional roots keep large server results in the data archive."""
    study=Path(study).resolve();output=Path(output).resolve()
    roots={name:study/name for name in ('evidence','case_studies','manuscript','formal_execution','performance_workflows') if (study/name).is_dir()}
    for item in extra_sources:
        if '=' not in item:raise ValueError('research source must be ARCHIVE_NAME=PATH')
        name,value=item.split('=',1)
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',name) or name in roots or name=='legacy_frozen_tables':
            raise ValueError('research source names must be safe, unique archive directories')
        path=Path(value).resolve()
        if not path.is_dir():raise FileNotFoundError(path)
        roots[name]=path
    for path in roots.values():
        if output.is_relative_to(path.resolve()):raise ValueError('archive output cannot be inside an included research source')
    return roots


def write_research_file(archive,source,archived):
    """Stream large files and hash the actual archived bytes; reject live changes."""
    source=Path(source);before=source.stat();digest=hashlib.sha256();size=0
    with source.open('rb') as incoming,archive.open(archived,'w',force_zip64=True) as outgoing:
        for chunk in iter(lambda:incoming.read(1024*1024),b''):
            outgoing.write(chunk);digest.update(chunk);size+=len(chunk)
    after=source.stat()
    if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino) or size!=before.st_size:
        raise RuntimeError('Research source changed while archiving: '+str(source))
    return dict(path=archived,bytes=size,sha256=digest.hexdigest())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--study-directory',type=Path,required=True)
    parser.add_argument('--output-directory',type=Path,required=True)
    parser.add_argument('--gate-status',type=Path,required=True)
    parser.add_argument('--research-source',action='append',default=[],metavar='ARCHIVE_NAME=PATH',
                        help='Additional explicit result directory, e.g. formal_suite=/server/project/results; repeat as needed')
    args=parser.parse_args()
    repo=args.repository.resolve();study=args.study_directory.resolve();out=args.output_directory.resolve()
    if run_git(repo,'status','--porcelain'):raise SystemExit('Commit the reviewed source changes before creating archives.')
    commit=run_git(repo,'rev-parse','HEAD')
    gates=json.loads(args.gate_status.read_text())
    roots=research_roots(study,args.research_source,out)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use an empty archive output directory.')
    out.mkdir(parents=True,exist_ok=True)
    code=out/f'RustBEARS-v1.2-development-{commit[:12]}-code.zip'
    subprocess.run(['git','-C',str(repo),'archive','--format=zip','--prefix=RustBEARS/',
                    '-o',str(code),commit],check=True)
    rows=[]
    data=out/'RustBEARS-v1.2-development-research-data.zip'
    with zipfile.ZipFile(data,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for directory,root in roots.items():
            for path in sorted(root.rglob('*')):
                if not path.is_file() or path.is_symlink() or '__pycache__' in path.parts or path.name=='.DS_Store':continue
                # Writing-strategy notes are not public research artifacts.
                if path.name=='project_context.md':continue
                relative=directory+'/'+path.relative_to(root).as_posix()
                rows.append(write_research_file(archive,path,relative))
        frozen_manifest=study/'evidence/frozen_17_tables.json'
        if frozen_manifest.exists():
            for item in json.loads(frozen_manifest.read_text())['files']:
                relative=Path(item['path'])
                source=(study.parent/relative).resolve()
                if (relative.is_absolute() or '..' in relative.parts
                        or not relative.parts[0].startswith('14_')
                        or not source.is_relative_to(study.parent)):
                    raise ValueError('Frozen table path is outside the declared protected stage.')
                digest=sha(source)
                if digest!=item['sha256']:raise ValueError(f'Protected table changed: {relative}')
                archived='legacy_frozen_tables/'+relative.as_posix()
                row=write_research_file(archive,source,archived)
                if row['sha256']!=digest:raise ValueError('Protected table changed during archiving')
                rows.append(dict(row,provenance='Unchanged read-only legacy supplementary table; see original table attribution.'))
        manifest=dict(status='development_snapshot_not_deposited',source_commit=commit,
            scientific_release_gates=gates,files=rows,research_source_roots={name:str(path) for name,path in roots.items()},
            licensing='Retain per-source license notices. Original BioGeoBEARS tutorial inputs retain their source attribution. Verify redistribution terms before a final Dryad deposit; this draft assigns no blanket license to third-party inputs.')
        archive.writestr('ARCHIVE_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        archive.writestr('README.md',
            '# RustBEARS development research snapshot\n\n'
            'This local archive is unpublished. No DOI has been assigned.\n'
            'Read ARCHIVE_MANIFEST.json and each experiment manifest before using results.\n'
            'Smoke tests, preliminary evaluations, failed attempts, and formal gate results\n'
            'retain their original labels. An included plan is not an executed experiment.\n'
            'Code, offline Rust dependencies, validation scripts and tutorials are in the companion code archive.\n')
    for path in (code,data):
        with zipfile.ZipFile(path) as archive:
            failed=archive.testzip()
            if failed:raise RuntimeError(f'Archive integrity failure: {path.name}: {failed}')
    index=dict(status='prepared_locally_unpublished',source_commit=commit,
        intended_destinations=dict(research_data='Dryad',code='Zenodo'),
        doi=dict(research_data=None,code=None),source_clean=True,
        gate_status=gates,archives=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in (code,data)],
        excluded=['installed runtimes','external downloaded software binaries','build caches','strategic writing notes'],
        research_files=len(rows))
    (out/'archive_index.json').write_text(json.dumps(index,ensure_ascii=False,indent=2)+'\n')
    (out/'README.md').write_text(
        '# Local archive preparation\n\nBoth ZIP files passed archive CRC verification.\n'
        'The code archive is an exact Git snapshot. The research archive records a SHA-256 for every input.\n'
        'This is a development snapshot. Formal calibration and release gates remain recorded in archive_index.json.\n'
        'Before deposit: complete all gates, replace provisional manuscript material, verify author metadata\n'
        'and redistribution licenses, then prepare a fresh immutable archive. No deposit or DOI exists yet.\n')
    print(json.dumps(index,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
