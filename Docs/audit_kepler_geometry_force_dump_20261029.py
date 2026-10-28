"""Compare native read-only face CSR and per-cell force/work dumps."""
from pathlib import Path
import struct

import numpy as np

MAGIC = b'GFSFG01\0'


def read_dump(path):
    path = Path(path)
    with path.open('rb') as stream:
        raw = stream.read(56)
    if len(raw) != 56:
        raise ValueError('Incomplete face/force header')
    magic, nbp, total, faces, call, rank, version = struct.unpack('<8s6Q', raw)
    if (magic != MAGIC or version != 1 or not 0 < nbp <= total
            or not 0 < faces <= 32*nbp or not 1 <= call <= 4):
        raise ValueError('Unsupported face/force layout')
    fields = [('id', '<i8', total), ('x', '<f8', total), ('y', '<f8', total),
              ('offset', '<i4', nbp+1)]
    fields += [(name, '<f8', faces) for name in ['c1x', 'c1y', 'c2x', 'c2y']]
    fields += [(name, '<i4', faces) for name in ['neighbor', 'kp', 'km', 'is_ghost']]
    fields += [(name, '<f8', nbp) for name in ['ax', 'ay']]
    fields += [(name, '<f4', nbp) for name in ['die', 'dK', 'dt']]
    expected = 56 + sum(np.dtype(dtype).itemsize*count for _, dtype, count in fields)
    if path.stat().st_size != expected:
        raise ValueError('Incomplete or trailing face/force data')
    data = {'nbp': nbp, 'total': total, 'faces': faces, 'call': call, 'rank': rank}
    offset = 56
    for name, dtype, count in fields:
        data[name] = np.memmap(path, dtype=dtype, mode='r', offset=offset, shape=(count,))
        offset += np.dtype(dtype).itemsize*count
    o = data['offset']
    if o[0] != 0 or o[-1] != faces or np.any(np.diff(o) < 0):
        raise ValueError('Invalid CSR offsets')
    if not np.array_equal(np.sort(data['id'][:nbp]), np.arange(nbp)):
        raise ValueError('Incomplete owned cell IDs')
    for name in ['neighbor', 'kp', 'km']:
        if np.any((data[name] < -1) | (data[name] >= total)):
            raise ValueError('Invalid neighbor index')
    for name, dtype, _ in fields:
        if dtype.startswith('<f') and not np.isfinite(data[name]).all():
            raise ValueError(f'Nonfinite {name}')
    return data


def exact(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def compare_dumps(a_path, b_path):
    a, b = read_dump(a_path), read_dump(b_path)
    checks = {'same_owned_size_call_rank': all(a[k] == b[k] for k in ['nbp', 'faces', 'call', 'rank'])}
    if not checks['same_owned_size_call_rank']:
        return {'passed': False, 'checks': checks, 'paths': [str(a_path), str(b_path)]}
    n = a['nbp']
    checks['owned_ID_order_identical'] = exact(a['id'][:n], b['id'][:n])
    checks['CSR_offsets_identical'] = exact(a['offset'], b['offset'])
    for field in ['x', 'y']:
        checks['owned_'+field+'_exact'] = exact(a[field][:n], b[field][:n])
    for field in ['c1x', 'c1y', 'c2x', 'c2y', 'is_ghost', 'ax', 'ay', 'die', 'dK', 'dt']:
        checks[field+'_exact'] = exact(a[field], b[field])
    # A reduced halo legitimately changes array indices. Compare support IDs
    # and image displacements, retaining the per-owner face order and shape.
    for field in ['neighbor', 'kp', 'km']:
        normalized = []
        for data in [a, b]:
            owner = np.repeat(np.arange(n), np.diff(data['offset']))
            index = data[field]
            valid = index >= 0
            ids = np.full(len(index), -1, dtype=np.int64)
            dx, dy = np.zeros(len(index)), np.zeros(len(index))
            ids[valid] = data['id'][index[valid]]
            dx[valid] = data['x'][index[valid]] - data['x'][owner[valid]]
            dy[valid] = data['y'][index[valid]] - data['y'][owner[valid]]
            normalized.append((ids, dx, dy))
        for component, left, right in zip(['ID', 'dx', 'dy'], *normalized):
            checks[field+'_'+component+'_exact'] = exact(left, right)
    return {'passed': bool(all(checks.values())), 'checks': checks,
            'paths': [str(a_path), str(b_path)], 'force_call': a['call'],
            'owned_cells': n, 'faces': a['faces'], 'padding_counts': [a['total']-n, b['total']-n]}
