"""Strict full native snapshot reader preserving noncontiguous global IDs."""
from pathlib import Path
import struct
import numpy as np
from audit_kepler_native_radial_20261026 import read_snapshot

ID_MASK = np.uint64((1 << 56) - 1)


def ordered_global(path, expected_ids=None):
    path = Path(path)
    with path.open('rb') as stream:
        magic, n, stride, stress = struct.unpack('<4i', stream.read(16))
    if ((magic, stride, stress) != (-1, 576, 1) or n <= 0
            or path.stat().st_size != 16 + n * stride + 52):
        raise ValueError(f'Incomplete native stress-RK4 snapshot: {path}')
    p, meta = read_snapshot(path)
    ids = p['u4if'] & ID_MASK
    order = np.argsort(ids)
    sorted_ids = ids[order]
    if np.any(sorted_ids[1:] == sorted_ids[:-1]):
        raise ValueError(f'Duplicate global IDs: {path}')
    if expected_ids is not None and not np.array_equal(sorted_ids, expected_ids):
        raise ValueError(f'Changed global ID set: {path}')
    return p[order], meta
