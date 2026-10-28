"""Reject any non-pointer restart change, including metadata and truncation."""
from pathlib import Path
import struct
import tempfile
import unittest

from audit_kepler_restart_bytes_20261029 import compare


class RestartByteTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.a, self.b = [Path(self.tmp.name)/name for name in ['a.dat', 'b.dat']]
        self.raw = (struct.pack('<4i', -1, 2, 576, 1)+bytes(2*576)
                    +struct.pack('<ddfiiddff', 2.66, 1e-4, 1e-4, 384, 384, 6., 6., 0., 0.))
        self.a.write_bytes(self.raw)
        self.b.write_bytes(self.raw)

    def mutate(self, offset):
        data = bytearray(self.raw)
        data[offset] ^= 1
        self.b.write_bytes(data)

    def test_identical(self):
        self.assertTrue(compare(self.a, self.b)['passed'])

    def test_next_pointer_only_is_ignored(self):
        self.mutate(16)
        row = compare(self.a, self.b)
        self.assertTrue(row['passed'])
        self.assertEqual(row['records_with_different_next_pointer'], 1)

    def test_velocity_and_local_dt_are_not_ignored(self):
        for offset in [56, 128, 552]:
            with self.subTest(offset=offset):
                self.mutate(16+offset)
                row = compare(self.a, self.b)
                self.assertFalse(row['passed'])
                self.assertEqual(row['changed_records_excluding_pointer'], 1)

    def test_full_precision_footer_change_rejected(self):
        self.mutate(16+2*576+8)
        row = compare(self.a, self.b)
        self.assertFalse(row['passed'])
        self.assertFalse(row['checks']['full_precision_restart_footer_exact'])

    def test_truncation_and_layout_mismatch_rejected(self):
        self.b.write_bytes(self.raw[:-1])
        with self.assertRaises(ValueError):
            compare(self.a, self.b)
        self.mutate(8)
        with self.assertRaises(ValueError):
            compare(self.a, self.b)


if __name__ == '__main__':
    unittest.main()
