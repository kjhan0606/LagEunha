"""Per-ring Voronoi cell geometry over time from khout.NNNNNN.dat snapshots
(W = 0 cells restricted to the generator's disc of radius rho, as with
EUNHA_VORO_RHO). For each ring: number of real faces, shortest real face /
dx, generator-to-centroid offset / dx, and the cell elongation
sqrt(lambda_max/lambda_min) of the cell's second-moment tensor.

Usage: python cell_geometry_rings.py RUNDIR [cx cy rho_over_dx dx]
"""
import glob
import os
import struct
import sys

import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon

RUN = sys.argv[1]
CX = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
RHO = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
DX = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0625
RINGS = [(0.0, 0.5), (0.5, 0.55), (0.55, 0.6), (0.6, 0.75), (0.75, 1.0), (1.0, 1.5), (1.5, 1.9), (1.9, 2.1)]


def load(p):
    with open(p, 'rb') as f:
        _, n, ps, _ = struct.unpack('iiii', f.read(16))
        a = np.frombuffer(f.read(n * ps), dtype=np.uint8).reshape(n, ps)
        t = struct.unpack('d', f.read()[:8])[0]
    d = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
    return t, d(24) - CX, d(32) - CY


snaps = sorted(glob.glob(os.path.join(RUN, 'khout.*.dat')))
pick = [s for k, s in enumerate(snaps) if k in (0, 4, 9, 14, 16, 19) or s == snaps[-1]]
print("# t ring N  <faces>  med(min face)/dx  med|x_gen - x_cent|/dx  med elongation")
for p in sorted(set(pick)):
    t, x, y = load(p)
    P = np.stack([x, y], 1)
    vor = Voronoi(P)
    R = RHO * DX
    r = np.hypot(x, y)
    rows = {k: [] for k in range(len(RINGS))}
    ridges = {}
    for (i, j), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridges.setdefault(i, []).append((j, rv))
        ridges.setdefault(j, []).append((i, rv))
    for i in range(len(x)):
        k = next((m for m, (lo, hi) in enumerate(RINGS) if lo <= r[i] < hi), None)
        if k is None:
            continue
        reg = vor.regions[vor.point_region[i]]
        disc = Point(P[i]).buffer(R, 256)
        if -1 in reg or len(reg) == 0:
            cell = disc
            for j, _ in ridges.get(i, []):
                m = 0.5 * (P[i] + P[j]); nrm = (P[j] - P[i]) / np.linalg.norm(P[j] - P[i])
                t_ = np.array([-nrm[1], nrm[0]]); big = 10 * R
                hp = Polygon([m + big * t_, m + big * t_ - big * nrm, m - big * t_ - big * nrm, m - big * t_])
                cell = cell.intersection(hp)
        else:
            cell = Polygon(vor.vertices[reg]).intersection(disc)
        if cell.is_empty:
            continue
        nf, fmin = 0, np.inf
        for j, rv in ridges.get(i, []):
            m = 0.5 * (P[i] + P[j]); nrm = (P[j] - P[i]) / np.linalg.norm(P[j] - P[i])
            seg = cell.exterior.intersection(Polygon([m - 1e-9 * nrm + 10 * R * np.array([-nrm[1], nrm[0]]),
                                                       m + 1e-9 * nrm + 10 * R * np.array([-nrm[1], nrm[0]]),
                                                       m + 1e-9 * nrm - 10 * R * np.array([-nrm[1], nrm[0]]),
                                                       m - 1e-9 * nrm - 10 * R * np.array([-nrm[1], nrm[0]])]))
            L = seg.length
            if L > 1e-9 * R:
                nf += 1; fmin = min(fmin, L)
        c = np.array(cell.centroid.coords[0])
        xs = np.array(cell.exterior.coords) - c
        I = np.cov(xs.T) if len(xs) > 2 else np.eye(2)
        ev = np.linalg.eigvalsh(I)
        rows[k].append((nf, fmin / DX if np.isfinite(fmin) else 0, np.linalg.norm(P[i] - c) / DX,
                        np.sqrt(ev[-1] / max(ev[0], 1e-30))))
    for k, (lo, hi) in enumerate(RINGS):
        if rows[k]:
            a = np.array(rows[k])
            print(f"{t:6.2f} [{lo:4.2f},{hi:4.2f}) {len(a):4d}  {a[:, 0].mean():5.2f}  {np.median(a[:, 1]):8.3f}  "
                  f"{np.median(a[:, 2]):8.3f}  {np.median(a[:, 3]):7.3f}")
