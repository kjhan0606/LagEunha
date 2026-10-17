"""Unlimited least-squares velocity gradient (Delaunay neighbours within the
cap reach, weight 1/d^2) vs the code's limited gradient vs the exact
Keplerian gradient, per radial ring.

If the unlimited LSF is close to exact where the code is not, the limiter
cut the gradient; if the unlimited LSF is also off, the one-sided stencil
itself is biased.
Usage: python lsf_unlimited_check.py SNAPSHOT [cx cy eps dmax]
"""
import struct
import sys

import numpy as np
from scipy.spatial import Delaunay

P = sys.argv[1]
CX = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
EPS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
DMAX = float(sys.argv[5]) if len(sys.argv) > 5 else 0.44

with open(P, 'rb') as f:
    _, n, ps, _ = struct.unpack('iiii', f.read(16))
    a = np.frombuffer(f.read(n * ps), dtype=np.uint8).reshape(n, ps)
    t = struct.unpack('d', f.read()[:8])[0]
d = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
x, y, vx, vy = d(24) - CX, d(32) - CY, d(56), d(64)
Gc = np.stack([d(352), d(360), d(368), d(376)], 1)
r = np.hypot(x, y)
Om = (r * r + EPS * EPS) ** -0.75
dOm = -1.5 * r * (r * r + EPS * EPS) ** -1.75
Ge = np.stack([-y * dOm * x / r, -Om - y * dOm * y / r, Om + x * dOm * x / r, x * dOm * y / r], 1)

tri = Delaunay(np.stack([x, y], 1))
indptr, nb = tri.vertex_neighbor_vertices
Gl = np.zeros((n, 4))
for i in range(n):
    js = nb[indptr[i]:indptr[i + 1]]
    dx, dy = x[js] - x[i], y[js] - y[i]
    dd = np.hypot(dx, dy)
    k = dd < DMAX
    js, dx, dy, dd = js[k], dx[k], dy[k], dd[k]
    w = 1.0 / dd ** 2
    M = np.array([[np.sum(w * dx * dx), np.sum(w * dx * dy)], [np.sum(w * dx * dy), np.sum(w * dy * dy)]])
    if len(js) < 2 or abs(np.linalg.det(M)) < 1e-30:
        Gl[i] = np.nan
        continue
    Mi = np.linalg.inv(M)
    for c, v in ((0, vx), (2, vy)):
        b = np.array([np.sum(w * dx * (v[js] - v[i])), np.sum(w * dy * (v[js] - v[i]))])
        Gl[i, c:c + 2] = Mi @ b

err = lambda G: np.linalg.norm(G - Ge, axis=1) / np.linalg.norm(Ge, axis=1)
er = np.stack([x / r, y / r], 1); ep = np.stack([-y / r, x / r], 1)
def shear(G):
    return (ep[:, 0] * (G[:, 0] * er[:, 0] + G[:, 1] * er[:, 1])
            + ep[:, 1] * (G[:, 2] * er[:, 0] + G[:, 3] * er[:, 1]))
se = shear(Ge)
print(f"# {P} t={t:.3f}")
print("#  ring           N   err code   err LSF-unlim   shear code/ex   shear LSF/ex   div code/sh   div LSF/sh")
edges = [0.49, 0.54, 0.6, 0.75, 1.0, 1.5, 1.9, 2.0]
for lo, hi in zip(edges[:-1], edges[1:]):
    s = (r >= lo) & (r < hi) & np.isfinite(Gl[:, 0])
    if s.sum() == 0:
        continue
    print(f"  [{lo:4.2f},{hi:4.2f}) {s.sum():5d}   {np.median(err(Gc)[s]):8.3e}   {np.median(err(Gl)[s]):8.3e}"
          f"      {np.median(shear(Gc)[s] / se[s]):7.3f}        {np.median(shear(Gl)[s] / se[s]):7.3f}"
          f"      {np.median(np.abs((Gc[:, 0] + Gc[:, 3]) / se)[s]):8.3e}  {np.median(np.abs((Gl[:, 0] + Gl[:, 3]) / se)[s]):8.3e}")
