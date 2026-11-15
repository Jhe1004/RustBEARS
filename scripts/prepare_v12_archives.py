#!/usr/bin/env python3
"""Prepare local, explicitly unpublished code and research-data archives.

The code archive comes from one clean Git commit. Data files retain their
relative paths and SHA-256 provenance. This command never deposits an archive,
registers a DOI, changes repository visibility, or marks a release gate passed.
"""
import argparse
import hashlib
import json
import os
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


COMPILED_SUFFIXES={'.whl','.egg','.conda','.so','.dylib','.dll','.pyd',
                   '.exe','.o','.obj','.a','.lib','.pyc','.pyo'}


def research_exclusion_reason(source,root):
    """Exclude distribution artifacts, not scientific binary data/checkpoints.

    Directory names alone do not identify an R/conda/virtualenv runtime: use
    installation markers. The private figure-QA copies are a separate,
    explicitly unreviewed redistribution scope; their inventories stay included.
    """
    source=Path(source);root=Path(root);relative=source.relative_to(root)
    if source.is_symlink():return 'symlink_not_followed'
    if '__pycache__' in relative.parts:return 'python_bytecode_cache'
    if source.name=='.DS_Store':return 'platform_metadata'
    if source.name=='project_context.md':return 'private_writing_strategy'
    if 'qa_toolchain' in relative.parts:return 'private_qa_toolchain_license_not_reviewed'
    if source.name=='private_toolchain.tar.gz' and source.parent.name=='figure_qa':
        return 'private_qa_toolchain_license_not_reviewed'
    if source.suffix.lower() in COMPILED_SUFFIXES or re.search(r'\.so(?:\.\d+)+$',source.name,re.I):
        return 'compiled_binary_or_distribution_package'
    for directory in (source.parent,*source.parent.parents):
        if not directory.is_relative_to(root):break
        if (directory/'pyvenv.cfg').is_file() or (directory/'conda-meta/history').is_file():
            return 'installed_runtime_tree'
        if directory.name in ('site-packages','dist-packages'):
            return 'installed_python_packages'
        if (directory/'DESCRIPTION').is_file() and (directory/'Meta/package.rds').is_file():
            return 'installed_r_package'
    # Extensionless executables (e.g. RevBayes) otherwise evade suffix rules.
    # Do not exclude arbitrary .bin/.npz/.rds files merely for being binary.
    if source.suffix:return None
    with source.open('rb') as stream:
        magic=stream.read(4)
        if magic in (b'\x7fELF',b'\xfe\xed\xfa\xce',b'\xce\xfa\xed\xfe',
                     b'\xfe\xed\xfa\xcf',b'\xcf\xfa\xed\xfe',b'\xca\xfe\xba\xbe',
                     b'\xbe\xba\xfe\xca'):
            return 'native_executable_header'
        if magic[:2]==b'MZ':
            stream.seek(60);offset=stream.read(4)
            if len(offset)==4:
                stream.seek(int.from_bytes(offset,'little'))
                if stream.read(4)==b'PE\x00\x00':return 'native_executable_header'
    return None


def excluded_research_file(source,archived,reason):
    """Inventory skipped bytes without copying, modifying or following links."""
    source=Path(source);before=source.lstat();digest=hashlib.sha256();size=0
    if source.is_symlink():
        payload=os.fsencode(os.readlink(source));digest.update(payload);size=len(payload)
        kind='symlink_target_text'
    else:
        with source.open('rb') as incoming:
            for chunk in iter(lambda:incoming.read(1024*1024),b''):
                digest.update(chunk);size+=len(chunk)
        kind='file'
    after=source.lstat()
    if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino) or size!=before.st_size:
        raise RuntimeError('Excluded research source changed during inventory: '+str(source))
    return dict(path=archived,source_path=str(source.absolute()),bytes=size,
                sha256=digest.hexdigest(),reason=reason,content_kind=kind)


def archive_research_roots(archive,roots):
    """Return included and excluded inventories for exactly the selected roots."""
    rows=[];excluded=[]
    for directory,root in roots.items():
        for path in sorted(root.rglob('*')):
            if not path.is_file() and not path.is_symlink():continue
            relative=directory+'/'+path.relative_to(root).as_posix()
            reason=research_exclusion_reason(path,root)
            if reason:excluded.append(excluded_research_file(path,relative,reason))
            else:rows.append(write_research_file(archive,path,relative))
    return rows,excluded


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
    data=out/'RustBEARS-v1.2-development-research-data.zip'
    with zipfile.ZipFile(data,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        rows,excluded=archive_research_roots(archive,roots)
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
            excluded_files=excluded,excluded_file_count=len(excluded),excluded_bytes=sum(item['bytes'] for item in excluded),
            licensing='Retain per-source license notices. Original BioGeoBEARS tutorial inputs retain their source attribution. Verify redistribution terms before a final Dryad deposit; this draft assigns no blanket license to third-party inputs.')
        archive.writestr('ARCHIVE_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        archive.writestr('README.md',
            '# RustBEARS development research snapshot\n\n'
            'This local archive is unpublished. No DOI has been assigned.\n'
            'Read ARCHIVE_MANIFEST.json and each experiment manifest before using results.\n'
            'Smoke tests, preliminary evaluations, failed attempts, and formal gate results\n'
            'retain their original labels. An included plan is not an executed experiment.\n'
            'Excluded runtime/distribution/private-tool files remain at their original paths;\n'
            'ARCHIVE_MANIFEST.json records their paths, sizes, hashes and reasons.\n'
            'Code, offline Rust dependencies, validation scripts and tutorials are in the companion code archive.\n')
    for path in (code,data):
        with zipfile.ZipFile(path) as archive:
            failed=archive.testzip()
            if failed:raise RuntimeError(f'Archive integrity failure: {path.name}: {failed}')
    index=dict(status='prepared_locally_unpublished',source_commit=commit,
        intended_destinations=dict(research_data='Dryad',code='Zenodo'),
        doi=dict(research_data=None,code=None),source_clean=True,
        gate_status=gates,archives=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in (code,data)],
        excluded_files=excluded,excluded_file_count=len(excluded),excluded_bytes=sum(item['bytes'] for item in excluded),
        exclusion_scope='Research roots only: known compiled suffixes/extensionless native executable headers, structurally identified installed runtimes/packages, private qa_toolchain copies and figure_qa/private_toolchain.tar.gz, symlinks, Python caches and private writing notes. Other archives/data retain their own licenses; this is not a general license audit.',
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
