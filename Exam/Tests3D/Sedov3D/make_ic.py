#!/usr/bin/env python3
"""3D Sedov-Taylor blast IC (periodic unit cube, gamma = 5/3).

rho = 1, P_amb = 1e-5 (>= 10x the 2D code's P = 1e-6 floor), v = 0.
E = 1 is deposited as a top hat over the N_hot particles nearest to the box
centre (default 64, the same count as the 2D Sedov IC in Exam/KH/util.c,
which spreads E over 64 central particles instead of one cell). The exact
number of hot particles (ties on a lattice) is reported and used.

Analytic shock radius R(t) = xi0 (E t^2/rho)^(1/5), xi0 = 1.15167 for
gamma = 5/3. Default t_end = 0.05 -> R = 0.3475, well inside the box.

--lattice cubic (default, cell-centred; the centre is a lattice vertex, so
the hot region is 8-fold symmetric) or glass (repulsive relaxation; see
common/lattice.py). The README/summary records which one was used.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lattice  # noqa: E402
import ic_common  # noqa: E402
import sedov_exact  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--lattice", default="cubic", choices=["cubic", "glass"])
    ap.add_argument("--E", type=float, default=1.0)
    ap.add_argument("--nhot", type=int, default=64)
    ap.add_argument("--pamb", type=float, default=1e-5)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--tend", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args(argv)
    n, g = a.n, a.gamma
    box = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    x, y, z, info = lattice.positions(a.lattice, n, n, n, box, seed=a.seed)
    N = len(x)
    V = 1.0 / N
    rho = np.ones(N)
    mass = rho * V
    r = np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2)
    order = np.argsort(r, kind="stable")
    rcut = r[order[a.nhot - 1]]
    hot = r <= rcut * (1 + 1e-9)           # include lattice ties
    nhot = int(hot.sum())
    u = np.full(N, a.pamb / ((g - 1) * 1.0))
    e_amb_hot = (mass[hot] * u[hot]).sum()
    u[hot] += a.E / mass[hot].sum()        # top hat in specific energy
    d = dict(x=x, y=y, z=z, vx=np.zeros(N), vy=np.zeros(N), vz=np.zeros(N),
             mass=mass, u=u, rho=rho, vol=np.zeros(N), id=np.arange(N, dtype=np.int64))
    meta = dict(gamma=g, box=box, periodic=(1, 1, 1))
    xi0 = sedov_exact.sedov_xi0(g)
    R = sedov_exact.shock_radius(a.tend, g, a.E, 1.0)
    params = dict(test="Sedov3D", n=n, lattice=info, E=a.E, nhot=nhot,
                  r_hot=float(rcut), pamb=a.pamb, rho0=1.0, center=[0.5, 0.5, 0.5],
                  t_end=a.tend, xi0=xi0, R_shock_tend=R,
                  E_thermal_ambient=float((mass * a.pamb / ((g - 1))).sum()))
    tot, bad = ic_common.report(d, meta)
    print("hot particles: %d within r<=%.4f (%.2f dx); P_hot/P_amb=%.3e" %
          (nhot, rcut, rcut * n, (g - 1) * u[hot][0] / a.pamb))
    print("Sedov: xi0=%.5f  R(t_end=%.3f)=%.4f  Eint-E=%.3e (ambient)" %
          (xi0, a.tend, R, tot["Eint"] - a.E))
    if R > 0.45:
        print("WARNING: shock reaches the periodic boundary before t_end")
    out = a.out or "sedov3d_n%d_%s.l3d" % (n, a.lattice)
    ic_common.write(out, d, meta, params, tot)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
