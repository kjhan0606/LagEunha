"""Keep matched large ends separate from extra floating-point contact ends.

A three-record ID pair is accepted only if exactly two non-tiny records
are opposite ends of the same face. Large split faces still fail.
"""
from pathlib import Path
import argparse
import json
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('case', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
f = np.loadtxt(args.case/'face_dump_r0_call1.tsv', skiprows=1)
a, b = f[:, 4].astype(np.int64), f[:, 5].astype(np.int64)
lo, hi = np.minimum(a, b), np.maximum(a, b)
order = np.lexsort((hi, lo))
starts = np.r_[0, np.flatnonzero((np.diff(lo[order]) != 0) | (np.diff(hi[order]) != 0))+1]
counts = np.diff(np.r_[starts, len(f)])
length = np.linalg.norm(f[:, 10:12], axis=1)
tolerance = 512*np.finfo(float).eps*24
bad, extras, large_pairs = [], [], []
for start, count in zip(starts[counts != 2], counts[counts != 2]):
    ix = order[start:start+count]
    large = ix[length[ix] > tolerance]
    small = ix[length[ix] <= tolerance]
    extras.extend(small.tolist())
    valid = len(large) == 0 or (len(large) == 2 and a[large[0]] == b[large[1]]
            and b[large[0]] == a[large[1]])
    if not valid:
        bad.append(dict(ids=[int(lo[ix[0]]), int(hi[ix[0]])],
                        lengths=length[ix].tolist(), records=f[ix].tolist()))
    elif len(large) == 2:
        large_pairs.append(large.tolist())
paired = starts[counts == 2]
pairs = np.column_stack((order[paired], order[paired+1]))
if large_pairs:
    pairs = np.vstack((pairs, np.array(large_pairs)))
p, q = pairs[:, 0], pairs[:, 1]
if not np.all((a[p] == b[q]) & (b[p] == a[q])):
    raise ValueError('matched ends have incompatible owner IDs')
geometry = np.linalg.norm(f[p, 10:12]+f[q, 10:12], axis=1)
c1 = .5*(f[p, 12:14]+f[p, 14:16])
c2 = .5*(f[q, 12:14]+f[q, 14:16])
delta = c1-c2
delta -= 24*np.rint(delta/24)
centroid_error = np.linalg.norm(delta, axis=1)
group_S = np.add.reduceat(f[order, 10:12], starts, axis=0)
group_force = np.add.reduceat(-f[order, 7, None]*f[order, 10:12], starts, axis=0)
group_work = np.add.reduceat(f[order, 9], starts)
owner = a
closure = np.hypot(np.bincount(owner, weights=f[:, 10]), np.bincount(owner, weights=f[:, 11]))
extra = np.array(extras, dtype=int)
row = dict(scope='initial native shared faces; large split groups rejected, extra FP ends explicitly retained',
    face_ends=len(f), groups_by_record_count={str(n): int(np.sum(counts == n)) for n in np.unique(counts)},
    unmatched_large_groups=bad, extra_tiny_ends=len(extra),
    extra_length_max=float(length[extra].max(initial=0)),
    extra_length_sum=float(length[extra].sum()),
    extra_absolute_pressure_force_sum=float(np.sum(abs(f[extra, 7])*length[extra])),
    extra_absolute_work_sum=float(np.sum(abs(f[extra, 9]))),
    paired_geometry_abs_max=float(geometry.max()), paired_centroid_periodic_max=float(centroid_error.max()),
    all_group_geometry_abs_max=float(np.linalg.norm(group_S, axis=1).max()),
    all_group_force_residual_abs_sum=float(np.linalg.norm(group_force, axis=1).sum()),
    all_group_work_residual_abs_sum=float(abs(group_work).sum()),
    surface_closure_max=float(closure.max()), geometry_tolerance=tolerance)
row['pass'] = bool(not bad and max(row[k] for k in
    ('extra_length_max', 'paired_geometry_abs_max', 'paired_centroid_periodic_max',
     'all_group_geometry_abs_max', 'surface_closure_max')) < tolerance)
args.output.write_text(json.dumps(row, indent=2)+'\n')
print(json.dumps({k: v for k, v in row.items() if k != 'unmatched_large_groups'}, indent=2))
raise SystemExit(0 if row['pass'] else 2)
