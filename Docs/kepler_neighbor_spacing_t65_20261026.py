"""Measure radial/azimuthal neighbor proximity in accepted C010 snapshots.

The native IC is Cartesian, so "same shell" is an explicitly defined radial
band rather than an IC ring. This is a descriptive diagnostic, not a force or
volume causality test.
"""

from pathlib import Path
import json
import sys

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026'))
from audit_varpress_checkpoint_20261026 import load_complete

H = 24. / 256.
INITIAL = ROOT / 'core_refinement/fine_integrated/khout.000001.dat'
CURRENT = ROOT / 'rk4_diagnostics/varpress_w0_geom_cfl_long_20261026/C010/to10/khout.039445.dat'


def summary(a):
    return dict(min=float(np.min(a)), p01=float(np.quantile(a, .01)),
                p05=float(np.quantile(a, .05)), median=float(np.median(a)),
                p95=float(np.quantile(a, .95)), max=float(np.max(a)))


def neighbor_metrics(xy, radii, selected, threshold):
    distance, index = cKDTree(xy).query(xy[selected], k=513)
    distance, index = distance[:, 1:], index[:, 1:]
    same = np.abs(radii[index] - radii[selected, None]) < threshold
    same_distance = np.min(np.where(same, distance, np.inf), axis=1)
    cross_distance = np.min(np.where(~same, distance, np.inf), axis=1)
    if not (np.isfinite(same_distance).all() and np.isfinite(cross_distance).all()):
        raise ValueError('512 neighbors did not cover both radial classes')
    return dict(nearest=distance[:, 0], same=same_distance,
                cross=cross_distance)


def main():
    p0, m0 = load_complete(INITIAL)
    pt, mt = load_complete(CURRENT)
    if not np.array_equal(p0['u4if'], pt['u4if']):
        raise ValueError('Snapshot cell IDs do not match')
    xy0 = np.column_stack((p0['x'], p0['y'])).astype(float)
    xyt = np.column_stack((pt['x'], pt['y'])).astype(float)
    radius = np.hypot(xy0[:, 0] - 12., xy0[:, 1] - 12.)
    selected = np.flatnonzero((radius >= 2.) & (radius < 2.25))
    density_error = abs(pt['den'][selected].astype(float) /
                        p0['den'][selected].astype(float) - 1.)
    nearest0 = neighbor_metrics(xy0, radius, selected, H / 4.)['nearest']
    nearestt = neighbor_metrics(xyt, radius, selected, H / 4.)['nearest']
    bands = []
    for threshold in (.015, .020, H / 4., .030):
        start = neighbor_metrics(xy0, radius, selected, threshold)
        end = neighbor_metrics(xyt, radius, selected, threshold)
        bands.append(dict(radial_band_width=threshold,
                          initial_same_median=float(np.median(start['same'])),
                          initial_cross_median=float(np.median(start['cross'])),
                          initial_cross_closer_fraction=float(np.mean(start['cross'] < start['same'])),
                          current_same_median=float(np.median(end['same'])),
                          current_cross_median=float(np.median(end['cross'])),
                          current_cross_closer_fraction=float(np.mean(end['cross'] < end['same']))))
    decile = np.quantile(nearestt, .1)
    close = nearestt <= decile
    result = dict(initial=str(INITIAL), current=str(CURRENT),
                  time=float(mt['t']), initial_grid_spacing=H,
                  target='2 <= initial r < 2.25', cells=int(len(selected)),
                  initial_distinct_rounded_radii=int(len(np.unique(np.round(radius[selected], 6)))),
                  initial_nearest_distance=summary(nearest0),
                  current_nearest_distance=summary(nearestt),
                  current_fraction_nearest_below_half_initial_spacing=float(np.mean(nearestt < H / 2.)),
                  current_nearest_vs_absolute_density_error_spearman=float(spearmanr(nearestt, density_error).statistic),
                  closest_decile_cutoff=float(decile), closest_decile_cells=int(close.sum()),
                  closest_decile_mean_absolute_density_error=float(np.mean(density_error[close])),
                  remaining_mean_absolute_density_error=float(np.mean(density_error[~close])),
                  radial_band_tests=bands,
                  limitation='Cartesian IC has no physical shells; radial bands are descriptive. Correlation is not proof of causality.')
    output = Path(__file__).with_name('KEPLER_NEIGHBOR_SPACING_T65_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
