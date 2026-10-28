import struct
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from audit_kepler_geometry_force_dump_20261029 import MAGIC, compare_dumps, read_dump


def fixture(path, expanded=False, wrong_image=False, wrong_force=False):
    n, nf = 2, 2
    ids = [0, 1, 0, 1] if expanded else [0, 1, 1]
    x = [1., 2., 7., -4.] if expanded else [1., 2., -4.]
    y = [1.]*len(ids)
    index = 3 if expanded else 2
    if wrong_image: x[index] += 6.
    fields = [(ids, '<i8'), (x, '<f8'), (y, '<f8'), ([0, 1, 2], '<i4')]
    fields += [([i, i+.5], '<f8') for i in range(4)]
    fields += [([index, 0], '<i4') for _ in range(3)] + [([0, 0], '<i4')]
    fields += [([1., 2.], '<f8'), ([3., 5. if wrong_force else 4.], '<f8')]
    fields += [([0., 0.], '<f4'), ([0., 0.], '<f4'), ([.1, .1], '<f4')]
    with path.open('wb') as f:
        f.write(struct.pack('<8s6Q', MAGIC, n, len(ids), nf, 1, 0, 1))
        for values, dtype in fields: f.write(np.array(values, dtype=dtype).tobytes())


class DumpTests(unittest.TestCase):
    def test_padding_index_change_preserves_geometry_and_force(self):
        with TemporaryDirectory() as directory:
            a, b = Path(directory)/'a', Path(directory)/'b'
            fixture(a); fixture(b, expanded=True)
            self.assertTrue(compare_dumps(a, b)['passed'])

    def test_same_ID_with_wrong_periodic_image_rejected(self):
        with TemporaryDirectory() as directory:
            a, b = Path(directory)/'a', Path(directory)/'b'
            fixture(a); fixture(b, expanded=True, wrong_image=True)
            result = compare_dumps(a, b)
            self.assertFalse(result['passed'])
            self.assertFalse(result['checks']['neighbor_dx_exact'])

    def test_force_change_rejected(self):
        with TemporaryDirectory() as directory:
            a, b = Path(directory)/'a', Path(directory)/'b'
            fixture(a); fixture(b, wrong_force=True)
            self.assertFalse(compare_dumps(a, b)['checks']['ay_exact'])

    def test_truncated_or_wrong_format_rejected(self):
        with TemporaryDirectory() as directory:
            a = Path(directory)/'a'
            fixture(a); raw = a.read_bytes()
            for broken in [raw[:-1], b'BADMAGIC'+raw[8:]]:
                a.write_bytes(broken)
                with self.assertRaises(ValueError): read_dump(a)


if __name__ == '__main__': unittest.main()
