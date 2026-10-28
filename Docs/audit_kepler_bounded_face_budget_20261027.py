"""Bound active-face and neighborhood arithmetic on a native Kepler state."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Delaunay


N = 55351
STRIDE = 576


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('snapshot', type=Path)
    ap.add_argument('--span', type=float, default=6.5)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    if args.span <= 0 or args.output.exists():
        raise ValueError('Require positive span and fresh output path')
    path = args.snapshot.resolve()
    with path.open('rb') as stream:
        head = struct.unpack('<4i', stream.read(16))
        stream.seek(16 + N * STRIDE)
        t = struct.unpack('<d', stream.read(8))[0]
    if head != (-1, N, STRIDE, 1) or path.stat().st_size != 16 + N * STRIDE + 52:
        raise ValueError('Incomplete or wrong native snapshot')
    dtype = np.dtype({'names': ('id', 'x', 'y', 'hydro_dt'),
                      'formats': ('<u8', '<f8', '<f8', '<f4'),
                      'offsets': (8, 24, 32, 128), 'itemsize': STRIDE})
    p = np.memmap(path, dtype=dtype, mode='r', offset=16, shape=(N,))
    ids = p['id'].astype(np.uint64) & np.uint64((1 << 56) - 1)
    if not np.array_equal(ids, np.arange(N, dtype=np.uint64)):
        raise ValueError('Expected native particle order with complete IDs')
    xy = np.column_stack((p['x'].astype(float), p['y'].astype(float)))
    hydro = p['hydro_dt'].astype(float)
    if not np.all(np.isfinite(hydro) & (hydro > 0)):
        raise ValueError('Invalid per-cell hydro CFL')
    r = np.hypot(xy[:, 0] - 3, xy[:, 1] - 3)
    local = np.minimum(hydro, .06 * (r*r + .01175**2)**.75)
    base = float(local.min())
    level = np.floor(np.log2(local/base) + 1e-12).astype(int)
    triangles = Delaunay(xy).simplices
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    global_steps = int(np.ceil(args.span/base))
    edge_level = np.minimum(level[edges[:, 0]], level[edges[:, 1]])
    global_face_updates = len(edges) * global_steps
    block_face_updates = int(np.sum(np.ceil(args.span/(base*np.exp2(edge_level))), dtype=np.int64))
    rings = []
    current = level.copy()
    for depth in range(4):
        block_sites = int(np.sum(np.ceil(args.span/(base*np.exp2(current))), dtype=np.int64))
        rings.append({'neighbor_rings': depth,
                      'level0_sites': int(np.sum(current == 0)),
                      'global_site_updates': N * global_steps,
                      'block_site_updates': block_sites,
                      'site_arithmetic_ratio': N * global_steps/block_sites})
        previous = current.copy()
        np.minimum.at(current, edges[:, 0], previous[edges[:, 1]])
        np.minimum.at(current, edges[:, 1], previous[edges[:, 0]])
    result = {
        'snapshot': str(path), 'snapshot_sha256': digest(path), 'time': t,
        'span_for_frozen_budget': args.span, 'cells': N,
        'triangles': len(triangles), 'unique_delaunay_edges': len(edges),
        'directed_delaunay_edge_estimate': 2*len(edges),
        'global_dt': base, 'global_steps': global_steps,
        'global_face_updates': global_face_updates,
        'block_face_updates': block_face_updates,
        'face_arithmetic_ratio': global_face_updates/block_face_updates,
        'edge_level_histogram': [{'level': int(i), 'edges': int(v)}
                                 for i, v in enumerate(np.bincount(edge_level)) if v],
        'active_neighborhood_budget': rings,
        'scope': 'Delaunay-edge and 0-3-ring arithmetic only at one frozen native state. '
                 'Delaunay may include zero-length Voronoi faces; no asynchronous '
                 'flux register, moving-face prediction, GPU geometry or runtime measured.',
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
