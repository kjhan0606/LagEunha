"""Synthetic fixtures for manual identity; no MPI/GPU/Slurm operations."""
import json
from pathlib import Path
import tempfile
import unittest

from audit_kepler_vertex_halo_manual_20261029 import launch_checks, FLAG_KEYS, MODE
from audit_kepler_vertex_halo_pair_20261029 import EXPECTED


class ContextTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = {'execution_mode': MODE, 'hostname': 'LagEunha',
            'GPU_UUID': 'GPU-00000000-0000-0000-0000-000000000000',
            'CUDA_VISIBLE_DEVICES': 'GPU-00000000-0000-0000-0000-000000000000',
            'manual_session_id': '00000000-0000-4000-8000-000000000000',
            'supervisor_pid': 123, 'SLURM_JOB_ID': ''}
        (self.root/'manual_execution.json').write_text(json.dumps(self.manifest))
        for name, expected in EXPECTED.items():
            (self.root/name).mkdir()
            (self.root/name/'launch_context.json').write_text(json.dumps({
                **self.manifest, **dict(zip(FLAG_KEYS, map(str, expected)))}))

    def change(self, name, key, value):
        path = self.root/name/'launch_context.json'
        row = json.loads(path.read_text())
        row[key] = value
        path.write_text(json.dumps(row))

    def test_complete_honest_context(self):
        self.assertTrue(all(launch_checks(self.root).values()))

    def test_missing_record_or_wrong_flags(self):
        self.change('vertex20', 'GFS_GPU_VERTEX_HALO_CERT', '0')
        self.assertFalse(all(launch_checks(self.root).values()))
        (self.root/'vertex40/launch_context.json').unlink()
        self.assertFalse(all(launch_checks(self.root).values()))

    def test_mixed_host_gpu_session(self):
        for key, value in [('hostname', 'other'), ('GPU_UUID', 'other'),
                           ('manual_session_id', 'other'), ('supervisor_pid', 124)]:
            with self.subTest(key=key):
                self.change('control40', key, value)
                self.assertFalse(all(launch_checks(self.root).values()))
                self.change('control40', key, self.manifest[key])

    def test_fabricated_slurm_or_missing_manifest_rejected(self):
        self.change('control40', 'SLURM_JOB_ID', '123456')
        self.assertFalse(all(launch_checks(self.root).values()))
        (self.root/'manual_execution.json').unlink()
        self.assertFalse(all(launch_checks(self.root).values()))

    def test_honestly_recorded_continuation(self):
        (self.root/'manual_continuation.json').write_text(json.dumps({
            **self.manifest, 'supervisor_pid': 124}))
        self.change('control40', 'supervisor_pid', 124)
        self.assertTrue(all(launch_checks(self.root).values()))
        continuation = {**self.manifest, 'supervisor_pid': 124, 'hostname': 'other'}
        (self.root/'manual_continuation.json').write_text(json.dumps(continuation))
        self.assertFalse(all(launch_checks(self.root).values()))


if __name__ == '__main__':
    unittest.main()
