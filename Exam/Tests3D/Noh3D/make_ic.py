#!/usr/bin/env python3
"""3D (spherical) Noh implosion IC.

Periodic cube of side L (default 6) centred at (L/2, L/2, L/2): rho = 1,
P = 1e-6 (the value used by the 2D Noh IC in Exam/KH/util.c), gamma = 5/3,
v = -r_hat (|v| = 1). Exact solution: stagnation shock at R = t/3, post-shock
rho = 64, P = 64/3; pre-shock rho = (1 + t/r)^2.

The periodic box edge is not the Noh inflow boundary. Its disturbance moves
inward with the cold inflow at |v| = 1, so at time t only r < L/2 - t is
valid (r < 1 at the default L = 6, t_end = 2, where R = 2/3). The analysis
restricts itself to that sphere.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lattice  # noqa: E402
import ic_common  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=128, help="particles per side")
    ap.add_argument("--L", type=float, default=6.0)
    ap.add_argument("--p0", type=float, default=1e-6)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--tend", type=float, default=2.0)
    ap.add_argument("--lattice", default="cubic", choices=["cubic", "glass"])
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args(argv)
    n, g, L = a.n, a.gamma, a.L
    box = (0.0, L, 0.0, L, 0.0, L)
    x, y, z, info = lattice.positions(a.lattice, n, n, n, box)
    N = len(x)
    V = L ** 3 / N
    c = 0.5 * L
    dx = np.stack([x - c, y - c, z - c], axis=1)
    r = np.linalg.norm(dx, axis=1)
    vel = -dx / np.maximum(r, 1e-300)[:, None]
    vel[r < 1e-12] = 0.0
    rho = np.ones(N)
    d = dict(x=x, y=y, z=z, vx=vel[:, 0], vy=vel[:, 1], vz=vel[:, 2],
             mass=rho * V, u=np.full(N, a.p0 / ((g - 1) * 1.0)), rho=rho,
             vol=np.zeros(N), id=np.arange(N, dtype=np.int64))
    meta = dict(gamma=g, box=box, periodic=(1, 1, 1))
    rs = 0.5 * (g - 1) * a.tend
    rvalid = 0.5 * L - a.tend
    params = dict(test="Noh3D", n=n, L=L, p0=a.p0, rho0=1.0, v0=1.0,
                  center=[c, c, c], t_end=a.tend, R_shock_tend=rs,
                  r_valid_tend=rvalid, lattice=info,
                  cells_per_Rs=rs / (L / n))
    tot, bad = ic_common.report(d, meta)
    print("Noh: R_s(t_end)=%.4f = %.1f dx; valid region r < %.3f; post-shock rho=%.1f" %
          (rs, rs / (L / n), rvalid, ((g + 1) / (g - 1)) ** 3))
    if rvalid <= 1.2 * rs:
        print("WARNING: valid region barely contains the shock; increase L")
    out = a.out or "noh3d_n%d_L%g.l3d" % (n, L)
    ic_common.write(out, d, meta, params, tot)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
