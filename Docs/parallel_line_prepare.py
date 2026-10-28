"""Batch the five audited straight-face velocity samples per cell.

The original per-face velocity calls and cached line constants remain the
reference. This helper only changes how independent sample evaluations are
batched before the original native cut/fill path consumes them.
"""
import os
import time
import numpy as np
from curved_disk_cell import ld, sub
from quadratic_velocity_reconstruction import velocity_samples
from kinetic_factor_rhs import ALPHA
import native_straight_faces as _native

_original_prepare = _native.prepare
_TIMES = (0.0, 0.5, 1.0, 0.25, 0.75)
_stats = {"status": "not_called"}


def prepare(cells, ids, positions, mass, factors, states, gen, priors,
            stage_cache, prior_cut):
    global _stats
    if stage_cache is None or hasattr(stage_cache, "native_straight_constants"):
        return _original_prepare(cells, ids, positions, mass, factors, states,
                                gen, priors, stage_cache, prior_cut)

    started = time.perf_counter()
    lookup = {int(g): i for i, g in enumerate(ids)}
    descriptors = []
    by_cell = [[] for _ in states]
    for i, cell in enumerate(cells):
        for line in cell.lines:
            if line.neighbor < 0:
                raise ValueError('unclipped square wall in physical mesh')
            j = lookup[line.neighbor]
            if ids[i] > ids[j]:
                continue
            key = (i, j)
            if key in stage_cache.lines:
                continue
            A = np.array([ld(v) for v in line.A])
            delta = np.array([ld(v) for v in sub(line.B, line.A)])
            _, normal, measure = line.quadrature(1)
            n = normal[0]
            length = measure[0]
            dx = positions[j] - positions[i]
            distance = float(np.linalg.norm(dx))
            if distance <= 0:
                raise ValueError('coincident fluid generators')
            points = np.asarray(
                [A + np.longdouble(t) * delta for t in _TIMES],
                dtype=np.longdouble)
            cuts = prior_cut(line, (i, j), priors)
            k = len(descriptors)
            descriptors.append(dict(i=i, j=j, key=key, A=A, delta=delta,
                normal=n, length=length, distance=distance, dx=dx,
                points=points, cuts=cuts, radial=[None, None]))
            by_cell[i].append((k, 0))
            by_cell[j].append((k, 1))

    cells_batched = 0
    for i, refs in enumerate(by_cell):
        if not refs:
            continue
        points = np.concatenate([descriptors[k]['points'] for k, _ in refs], axis=0)
        sampled = velocity_samples(states[i], factors[i], mass[i], points, ALPHA)[0]
        for order, (k, side) in enumerate(refs):
            lo = order * len(_TIMES)
            hi = lo + len(_TIMES)
            descriptors[k]['radial'][side] = sampled[lo:hi] @ descriptors[k]['normal']
        cells_batched += 1

    for item in descriptors:
        if item['radial'][0] is None or item['radial'][1] is None:
            raise RuntimeError('incomplete two-sided straight-face velocity batch')
        sample_normal = {
            t: (item['radial'][0][k], item['radial'][1][k])
            for k, t in enumerate(_TIMES)
        }
        stage_cache.lines[item['key']] = (
            item['A'], item['delta'], item['normal'], item['length'],
            item['distance'], item['dx'], sample_normal, item['cuts'])

    verify_count = max(0, int(os.environ.get('LAGFORCE_BATCH_FACE_VERIFY', '0')))
    exact = True
    max_abs = 0.0
    verify_faces = min(verify_count, len(descriptors))
    for item in descriptors[:verify_faces]:
        for k, t in enumerate(_TIMES):
            point = item['points'][k]
            for side, cell_i in enumerate((item['i'], item['j'])):
                reference = velocity_samples(states[cell_i], factors[cell_i],
                    mass[cell_i], point, ALPHA)[0][0] @ item['normal']
                actual = item['radial'][side][k]
                delta = float(np.max(np.abs(reference - actual)))
                max_abs = max(max_abs, delta)
                exact = exact and np.array_equal(reference, actual)

    prefetched = len(descriptors)
    before = stage_cache.line_hits
    result = _original_prepare(cells, ids, positions, mass, factors, states,
                               gen, priors, stage_cache, prior_cut)
    stage_cache.line_hits -= prefetched
    _stats = dict(status='passed_batched_line_sample_preparation',
        faces_batched=prefetched, cells_batched=cells_batched,
        point_evaluations=2 * len(_TIMES) * prefetched,
        scalar_velocity_calls_avoided=2 * len(_TIMES) * prefetched,
        batch_wall_seconds=time.perf_counter() - started,
        verification_faces=verify_faces,
        verification_bitwise_equal=exact if verify_faces else None,
        verification_max_abs_difference=max_abs if verify_faces else None,
        line_hit_count_adjustment=-prefetched,
        stage_cache_line_hits_before=before)
    return result


def provenance():
    return dict(_stats)
