"""Shared-face audit with explicit domain and force-call provenance."""
from pathlib import Path
import argparse
import json
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--call', type=int, default=1)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('cases', type=Path, nargs='+')
args = parser.parse_args()
results = {}
for case in args.cases:
    case = case.resolve()
    text = (case / 'params.dat').read_text()
    bounds = []
    for name in ('Xmin', 'Xmax', 'Ymin', 'Ymax'):
        match = re.search(r'^define KH ' + name + r'\s*=\s*([^\s#]+)', text, re.M)
        if match is None:
            raise ValueError(f'{case}: missing KH {name}; do not assume box24')
        bounds.append(float(match[1]))
    xmin, xmax, ymin, ymax = bounds
    domain = np.array([xmax-xmin, ymax-ymin])
    if not np.all(domain > 0):
        raise ValueError('invalid domain')
    files = sorted(case.glob(f'face_dump_r*_call{args.call}.tsv'))
    if not files:
        raise ValueError(f'{case}: missing requested face call {args.call}')
    f = np.vstack([np.loadtxt(p, skiprows=1, ndmin=2) for p in files])
    if np.any(f[:, 6] != 0):
        raise ValueError('wall faces require a separate contract')
    a, b = f[:, 4].astype(np.int64), f[:, 5].astype(np.int64)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    order = np.lexsort((hi, lo))
    starts = np.r_[0, np.flatnonzero((np.diff(lo[order]) != 0) | (np.diff(hi[order]) != 0))+1]
    counts = np.diff(np.r_[starts, len(f)])
    paired = starts[counts == 2]
    p, q = order[paired], order[paired+1]
    if not np.all((a[p] == b[q]) & (b[p] == a[q])):
        raise ValueError('duplicate face records do not identify opposite ends')
    s1, s2 = f[p, 10:12], f[q, 10:12]
    scale = np.maximum(np.linalg.norm(s1, axis=1), np.linalg.norm(s2, axis=1))
    absolute = np.linalg.norm(s1+s2, axis=1)
    relative = absolute/np.maximum(scale, 1e-300)
    c1 = .5*(f[p, 12:14]+f[p, 14:16])
    c2 = .5*(f[q, 12:14]+f[q, 14:16])
    dc = c1-c2
    dc -= domain*np.rint(dc/domain)
    orphan_count = int(np.sum(counts[counts != 2]))
    orphan_groups = [order[st:st+n] for st, n in zip(starts, counts) if n != 2]
    orphan = np.concatenate(orphan_groups) if orphan_groups else np.empty(0, dtype=int)
    orphan_lengths = np.linalg.norm(f[orphan,10:12], axis=1)
    force = -(f[p, 7, None]*s1 + f[q, 7, None]*s2)
    owners=f[:,4].astype(int)
    sx=np.bincount(owners,weights=f[:,10]);sy=np.bincount(owners,weights=f[:,11])
    surface_closure=np.hypot(sx,sy)
    row = {
        'case': str(case), 'force_call': args.call,
        'domain_source': str(case/'params.dat'), 'domain_lengths': domain.tolist(),
        'files': [str(p) for p in files], 'face_ends': len(f),
        'orphan_face_ends': orphan_count,
        'record_count_histogram': dict(zip(*[x.tolist() for x in np.unique(counts, return_counts=True)])),
        'geometry_abs_error_max': float(absolute.max()),
        'geometry_rel_error_max': float(relative.max()),
        'pairs_relative_gt_1e_minus8': int(np.sum(relative > 1e-8)),
        'orphan_length_max': float(orphan_lengths.max(initial=0)),
        'orphan_length_sum': float(orphan_lengths.sum()),
        'orphan_pressure_force_abs_sum': float(np.sum(abs(f[orphan,7])*orphan_lengths)),
        'orphan_work_abs_sum': float(abs(f[orphan,9]).sum()),
        'max_cell_surface_closure': float(surface_closure.max()),
        'tiny_face_caveat': 'Relative error alone can be large on negligible FP faces; inspect absolute geometry and force/work.',
        'centroid_periodic_distance_max': float(np.linalg.norm(dc, axis=1).max()),
        'pressure_difference_max': float(np.max(abs(f[p,7]-f[q,7]))),
        'sum_abs_pair_force_residual': float(np.linalg.norm(force,axis=1).sum()),
        'sum_abs_pair_work_residual': float(abs(f[p,9]+f[q,9]).sum()),
        'net_force_all_ends': (-np.sum(f[:,7,None]*f[:,10:12],axis=0)).tolist(),
        'sum_work_all_ends': float(f[:,9].sum()),
    }
    results[str(case)] = row
    print(case.name, 'orphans', orphan_count, 'max_absolute_geometry_error', row['geometry_abs_error_max'])
args.output.write_text(json.dumps(results, indent=2)+'\n')
