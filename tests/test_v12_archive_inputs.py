"""Additional server results must be included without path escapes or live writes."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from scripts.prepare_v12_archives import research_roots,write_research_file


class ArchiveInputsTests(unittest.TestCase):
    def test_default_protocols_and_explicit_remote_results_are_selected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);study=root/'study';study.mkdir()
            for name in ('evidence','formal_execution','performance_workflows'):(study/name).mkdir()
            server=root/'server_results';server.mkdir()
            selected=research_roots(study,[f'full_suite={server}'],root/'archives')
            self.assertEqual(set(selected),{'evidence','formal_execution','performance_workflows','full_suite'})
            self.assertEqual(selected['full_suite'],server.resolve())

    def test_unsafe_duplicate_and_recursive_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'evidence').mkdir()
            for alias in ('../escape','/absolute','evidence','legacy_frozen_tables'):
                with self.assertRaises(ValueError):research_roots(root,[f'{alias}={root}'],root/'archives')
            with self.assertRaisesRegex(ValueError,'inside'):
                research_roots(root,[f'full_suite={root}'],root/'archives')
            with self.assertRaises(FileNotFoundError):research_roots(root,[f'missing={root}/absent'],root/'archives')

    def test_streamed_archive_hash_matches_bytes_and_crc(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'events.bin';payload=bytes(range(256))*10000;source.write_bytes(payload)
            with zipfile.ZipFile(root/'data.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
                with patch.object(Path,'read_bytes',side_effect=AssertionError('whole-file read forbidden')):
                    row=write_research_file(archive,source,'full_suite/events.bin')
            self.assertEqual(row['sha256'],hashlib.sha256(payload).hexdigest())
            self.assertEqual(row['bytes'],len(payload))
            with zipfile.ZipFile(root/'data.zip') as archive:
                self.assertIsNone(archive.testzip());self.assertEqual(archive.read(row['path']),payload)


if __name__=='__main__':unittest.main()
