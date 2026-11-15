"""Archive-engineering fixtures; no scientific experiment or release is run."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

from scripts.prepare_v12_archives import (
    archive_research_roots, excluded_research_file, research_exclusion_reason,
)


class ArchiveExclusionsTests(unittest.TestCase):
    def write(self,root,name,payload):
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(payload);return p

    def test_cli_archive_retains_data_and_records_every_excluded_byte(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);repo=root/'repo';repo.mkdir()
            self.write(repo,'LICENSE',b'synthetic engineering fixture\n')
            for args in [('init','-q'),('add','LICENSE'),('-c','user.name=Archive Test','-c','user.email=archive-test@example.invalid','-c','commit.gpgsign=false','commit','-qm','fixture')]:
                subprocess.run(['git','-C',str(repo),*args],check=True,capture_output=True)
            study=root/'study';inputs=study/'evidence';extra=root/'extra_formal';extra.mkdir()
            included={
                'raw/replicates.npz':b'PK\x03\x04synthetic-array-bytes',
                'checkpoint/complete.rds':b'RDX3\nsynthetic-full-payload',
                'tables/all_failures.csv':b'rep,status\n1,failed\n',
                'logs/failed_attempt.log':b'error retained\n',
                'source_manifest.json':b'{"source":"fixture"}\n',
                'events.bin':bytes(range(64)),
                'raw/private_toolchain.tar.gz':b'ordinary scientific tar fixture',
                'qa_toolchain_inventory.json':b'[{"hash":"fixture"}]\n',
                'scripts/reproduce.py':b'print("fixture")\n',
            }
            excluded={
                'dist/package.whl':b'PK package fixture',
                'native/kernel.so':b'ELF fixture',
                'native/kernel.dll':b'PE fixture',
                'native/kernel.so.2':b'library fixture',
                'native/rb':b'\x7fELFexecutable fixture',
                'figures/qa_toolchain/private_validator.py':b'private-tool fixture',
                'tongfang/figure_qa/private_toolchain.tar.gz':b'private-tool archive fixture',
                'notes/project_context.md':b'private planning fixture',
            }
            for name,payload in included.items():self.write(inputs,name,payload)
            for name,payload in excluded.items():self.write(inputs,name,payload)
            self.write(extra,'checkpoint.json',b'{"failed":1}')
            self.write(extra,'bridge.dylib',b'dylib fixture')
            self.write(study,'gates.json',b'{"formal":false}')
            originals={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in study.rglob('*') if p.is_file()}
            originals.update({p:hashlib.sha256(p.read_bytes()).hexdigest() for p in extra.rglob('*') if p.is_file()})
            script=Path(__file__).resolve().parents[1]/'scripts/prepare_v12_archives.py'
            out=root/'archives'
            completed=subprocess.run([sys.executable,str(script),'--repository',str(repo),'--study-directory',str(study),'--output-directory',str(out),'--gate-status',str(study/'gates.json'),'--research-source','formal_suite='+str(extra)],capture_output=True,text=True,check=True)
            index=json.loads(completed.stdout)
            self.assertEqual(index['excluded_file_count'],len(excluded)+1)
            self.assertEqual(index['excluded_bytes'],sum(map(len,excluded.values()))+len(b'dylib fixture'))
            with zipfile.ZipFile(out/'RustBEARS-v1.2-development-research-data.zip') as archive:
                manifest=json.loads(archive.read('ARCHIVE_MANIFEST.json'))
                self.assertEqual(manifest['excluded_files'],index['excluded_files'])
                self.assertEqual(manifest['excluded_bytes'],index['excluded_bytes'])
                self.assertIsNone(archive.testzip())
                for name,payload in included.items():self.assertEqual(archive.read('evidence/'+name),payload)
                self.assertEqual(archive.read('formal_suite/checkpoint.json'),b'{"failed":1}')
                self.assertEqual({r['path'] for r in manifest['files']},set(archive.namelist())-{'README.md','ARCHIVE_MANIFEST.json'})
                for row in manifest['excluded_files']:
                    self.assertNotIn(row['path'],archive.namelist())
                    source=Path(row['source_path']);payload=source.read_bytes()
                    self.assertEqual(row['sha256'],hashlib.sha256(payload).hexdigest())
                    self.assertEqual(row['bytes'],len(payload));self.assertTrue(row['reason'])
            self.assertEqual(originals,{p:hashlib.sha256(p.read_bytes()).hexdigest() for p in originals})
            self.assertEqual(index['doi'],{'research_data':None,'code':None})

    def test_installed_package_markers_do_not_hide_source_or_scientific_checkpoints(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ['Rcheck/installed/DESCRIPTION','Rcheck/installed/Meta/package.rds','Rcheck/installed/R/code.rdb','Rcheck/source/DESCRIPTION','Rcheck/source/R/code.R','Rcheck/tests/comparison.rds','runtime/results/checkpoint.json','environment/pyvenv.cfg','environment/installed.txt','python/site-packages/module.py']:
                self.write(root,name,b'fixture')
            out=root.parent/(root.name+'-archive.zip')
            try:
                with zipfile.ZipFile(out,'w') as archive:rows,excluded=archive_research_roots(archive,{'study':root})
                self.assertEqual({r['path'] for r in rows},{'study/Rcheck/source/DESCRIPTION','study/Rcheck/source/R/code.R','study/Rcheck/tests/comparison.rds','study/runtime/results/checkpoint.json'})
                self.assertEqual(len(excluded),6)
                self.assertTrue(all(Path(r['source_path']).exists() for r in excluded))
            finally:out.unlink(missing_ok=True)

    def test_symlinks_are_not_followed_and_target_text_is_hashed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);target=self.write(root,'outside/data.bin',b'private target bytes')
            selected=root/'selected';selected.mkdir();link=selected/'external.bin'
            try:link.symlink_to(target)
            except OSError as error:self.skipTest('symlink creation not supported: '+str(error))
            reason=research_exclusion_reason(link,selected)
            self.assertEqual(reason,'symlink_not_followed')
            row=excluded_research_file(link,'selected/external.bin',reason)
            payload=os.fsencode(os.readlink(link))
            self.assertEqual(row['content_kind'],'symlink_target_text')
            self.assertEqual(row['sha256'],hashlib.sha256(payload).hexdigest())
            self.assertEqual(row['bytes'],len(payload));self.assertEqual(target.read_bytes(),b'private target bytes')

    def test_native_headers_preserve_unrelated_binary_and_source_archives(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            pe=bytearray(72);pe[:2]=b'MZ';pe[60:64]=(64).to_bytes(4,'little');pe[64:68]=b'PE\x00\x00'
            for name,payload in [('elf',b'\x7fELFtest'),('macho',b'\xcf\xfa\xed\xfetest'),('win',bytes(pe))]:
                self.assertEqual(research_exclusion_reason(self.write(root,name,payload),root),'native_executable_header')
            for name,payload in [('events.bin',b'\x7fELF arbitrary scientific binary'),('other.bin',b'MZnot-an-executable'),('posterior.npz',b'PK\x03\x04data'),('checkpoint.rds',b'RDX3\n'),('source.tar.gz',b'\x1f\x8bsource'),('run.log',b'ERROR retained')]:
                self.assertIsNone(research_exclusion_reason(self.write(root,name,payload),root))

    def test_mutated_excluded_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p=self.write(Path(directory),'native.so',b'binary')
            before=p.lstat();after=SimpleNamespace(st_size=before.st_size,st_mtime_ns=before.st_mtime_ns+1,st_ino=before.st_ino)
            # lstat is also used by is_symlink; keep that middle lookup genuine.
            with patch.object(Path,'is_symlink',return_value=False),patch.object(Path,'lstat',side_effect=[before,after]):
                with self.assertRaisesRegex(RuntimeError,'changed during inventory'):
                    excluded_research_file(p,'native.so','compiled_binary_or_distribution_package')


if __name__=='__main__':unittest.main()
