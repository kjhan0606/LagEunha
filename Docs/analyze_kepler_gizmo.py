"""Analyse khout.NNNNNN.dat snapshots of the GIZMO particle Kepler disk (Hopkins 2015 §4.2.4).

Measures numerical-viscosity damage per orbit: surface density profile vs the initial one,
angular momentum by radius, and |v_phi - v_K|/v_K near r = 1.
Usage: python analyze_kepler_gizmo.py RUNDIR [cx cy eps]
"""
import glob
import os
import struct
import sys

import numpy as np

RUN = sys.argv[1]
CX = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
EPS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
EDGES = np.array([0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0])
NEFF = int(sys.argv[5]) if len(sys.argv) > 5 else 64


def initial_reference():
    """GIZMO particle IC: NEFF^2 lattice over -2..2, kept for 0.5<=r<=2, m = dx^2, circular orbits."""
    dx = 4.0 / NEFF
    g = -2.0 + (np.arange(NEFF) + 0.5) * dx
    X, Y = np.meshgrid(g, g)
    R = np.hypot(X, Y).ravel()
    R = R[(R >= 0.5) & (R <= 2.0)]
    m = dx * dx
    lz = m * R * R * (R * R + EPS * EPS) ** -0.75
    return np.histogram(R, bins=EDGES)[0] * m, lz.sum(), np.histogram(R, bins=EDGES, weights=lz)[0]


def load(path):
    with open(path, 'rb') as f:
        first = struct.unpack('i', f.read(4))[0]
        if first == -1:
            npart, psize, _stress = struct.unpack('iii', f.read(12))
        else:
            npart, psize = first, 360
        data = f.read(npart * psize)
        trailer = f.read()
    a = np.frombuffer(data, dtype=np.uint8).reshape(npart, psize)
    col = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
    x, y, m, vx, vy = col(24), col(32), col(48), col(56), col(64)
    t = struct.unpack('d', trailer[0:8])[0]
    return t, x, y, m, vx, vy


snaps = sorted(glob.glob(os.path.join(RUN, 'khout.*.dat')))
if not snaps:
    sys.exit(f"no snapshots in {RUN}")
m0, Lz0, Lb0 = initial_reference()
m0 = np.where(m0 > 0, m0, np.nan)
print("# t  orbits(r=1)  L1(|vphi-vK|/vK, 0.8<r<1.2)  dLz/Lz0  mass(r-bin)/initial for bins "
      + " ".join(f"[{a:g},{b:g})" for a, b in zip(EDGES[:-1], EDGES[1:])))
for p in snaps:
    t, x, y, m, vx, vy = load(p)
    rx, ry = x - CX, y - CY
    r = np.hypot(rx, ry)
    vphi = (rx * vy - ry * vx) / np.maximum(r, 1e-30)
    vk = r * (r * r + EPS * EPS) ** -0.75
    lz = m * (rx * vy - ry * vx)
    mb = np.histogram(r, bins=EDGES, weights=m)[0]
    sel = (r > 0.8) & (r < 1.2)
    L1 = np.mean(np.abs(vphi[sel] - vk[sel]) / vk[sel]) if sel.any() else np.nan
    ratio = mb / m0
    print(f"{t:8.2f} {t/(2*np.pi):6.2f} {L1:.3e} {(lz.sum()-Lz0)/Lz0:+.3e} "
          + " ".join('   -  ' if not np.isfinite(q) else f"{q:6.3f}" for q in ratio)
          + f"  (mass r<0.5: {mb[0]:.3e}, [2,3): {mb[-1]:.3e}, r>=3: {m[r >= EDGES[-1]].sum():.3e})")
    Lb = np.histogram(r, bins=EDGES, weights=lz)[0]
    print("         dLz(r-bin) = " + " ".join(f"{a:+.2e}" for a in (Lb - Lb0)))
