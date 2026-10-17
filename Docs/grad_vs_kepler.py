"""Compare the code's limited velocity gradient (stress.gU**) with the exact
Keplerian gradient, per radial ring, in khout.NNNNNN.dat snapshots.

Exact field: v = Omega(r) (-y, x), Omega = (r^2 + eps^2)^(-3/4).
Reports, per ring, the median of |G_code - G_exact| / |G_exact| (Frobenius),
the median ratio of the shear component, and the median divergence (exact 0).
Usage: python grad_vs_kepler.py SNAPSHOT [cx cy eps]
"""
import struct
import sys

import numpy as np

P = sys.argv[1]
CX = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
EPS = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
R0 = 0.5173   # innermost initial ring of the 64^2 GIZMO lattice

with open(P, 'rb') as f:
    _, n, ps, _ = struct.unpack('iiii', f.read(16))
    a = np.frombuffer(f.read(n * ps), dtype=np.uint8).reshape(n, ps)
    t = struct.unpack('d', f.read()[:8])[0]
d = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
x, y = d(24) - CX, d(32) - CY
G = np.stack([d(352), d(360), d(368), d(376)], 1)          # gUxx gUxy gUyx gUyy
r = np.hypot(x, y)
Om = (r * r + EPS * EPS) ** -0.75
dOm = -1.5 * r * (r * r + EPS * EPS) ** -1.75                 # dOmega/dr
# v_x = -Om y, v_y = Om x ; d/dx = (x/r) d/dr
Ge = np.stack([-y * dOm * x / r, -Om - y * dOm * y / r,
               Om + x * dOm * x / r, x * dOm * y / r], 1)
err = np.linalg.norm(G - Ge, axis=1) / np.linalg.norm(Ge, axis=1)
# polar shear r dOmega/dr and divergence
er = np.stack([x / r, y / r], 1); ep = np.stack([-y / r, x / r], 1)
def comp(M, u, w):
    return (u[:, 0] * (M[:, 0] * w[:, 0] + M[:, 1] * w[:, 1])
            + u[:, 1] * (M[:, 2] * w[:, 0] + M[:, 3] * w[:, 1]))
shear_c = comp(G, ep, er); shear_e = comp(Ge, ep, er)       # d v_phi / d r
div_c = G[:, 0] + G[:, 3]
print(f"# {P}  t={t:.3f}")
print("#  ring            N   med|dG|/|G|  med(dvphi/dr code/exact)  med|div|/|shear|")
edges = [0, 0.49, 0.54, 0.6, 0.75, 1.0, 1.5, 1.9, 2.0, 2.3]
for lo, hi in zip(edges[:-1], edges[1:]):
    s = (r >= lo) & (r < hi)
    if s.sum() == 0:
        continue
    print(f"  [{lo:4.2f},{hi:4.2f})  {s.sum():5d}   {np.median(err[s]):9.3e}   "
          f"{np.median(shear_c[s] / shear_e[s]):9.3f}   {np.median(np.abs(div_c[s] / shear_e[s])):9.3e}")
