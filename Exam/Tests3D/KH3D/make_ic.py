#!/usr/bin/env python3
"""3D Kelvin-Helmholtz IC: McNally, Lyra & Passy (2012, ApJS 201, 18)
smoothed shear layer, extended uniformly in z.

Box [0,1] x [0,1] x [0,Lz], periodic. gamma = 5/3, P = 2.5,
rho1 = 1 (outer), rho2 = 2 (inner, 1/4 < y < 3/4), U1 = +0.5, U2 = -0.5,
exponential smoothing length L = 0.025 (McNally Eq. 1-2), seed
vy = 0.01 sin(4 pi x) (z-uniform, as in 2D). --znoise adds a random vz of
that rms (seeded) so genuinely 3D modes can develop; default 0.

Particles: cell-centred cubic lattice, n per unit length, mass = rho(y) dV
(the 2D KH IC in Exam/KH/util.c also uses a lattice with mass = rho*meanvol).
Default Lz = 1/4 (nz = n/4) keeps 128 per unit length at 1/4 of the cube's
cost; the linear mode is z-independent, so the growth rate is the same.
Linear growth rate of the seeded mode (k = 4 pi) from the compressible
eigen-solver common/kh_linear.py: sigma = 2.83 (the sharp incompressible
interface would give 5.92).
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lattice  # noqa: E402
import ic_common  # noqa: E402
import kh_linear  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=128, help="particles per unit length")
    ap.add_argument("--lz", type=float, default=0.25)
    ap.add_argument("--amp", type=float, default=0.01)
    ap.add_argument("--znoise", type=float, default=0.0)
    ap.add_argument("--L", type=float, default=0.025)
    ap.add_argument("--P0", type=float, default=2.5)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--tend", type=float, default=2.0)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--no-linear", action="store_true", help="skip the eigen-solve")
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args(argv)
    n, g = a.n, a.gamma
    nz = max(1, int(round(n * a.lz)))
    lz = nz / n
    box = (0.0, 1.0, 0.0, 1.0, 0.0, lz)
    x, y, z, info = lattice.cubic_lattice(n, n, nz, box) + (dict(kind="cubic"),)
    N = len(x)
    V = lz / N
    rho, U, _, _ = kh_linear.mcnally_profile(y, L=a.L)
    vy = a.amp * np.sin(4 * np.pi * x)
    vz = np.zeros(N)
    if a.znoise > 0:
        vz = a.znoise * np.random.default_rng(a.seed).standard_normal(N)
    u = a.P0 / ((g - 1) * rho)
    d = dict(x=x, y=y, z=z, vx=U, vy=vy, vz=vz, mass=rho * V, u=u, rho=rho,
             vol=np.zeros(N), id=np.arange(N, dtype=np.int64))
    meta = dict(gamma=g, box=box, periodic=(1, 1, 1))
    sig = None
    if not a.no_linear:
        sig, w, s0 = kh_linear.growth_rate(L=a.L, P0=a.P0, gamma=g, n=512)
    params = dict(test="KH3D", n=n, nz=nz, lz=lz, amp=a.amp, znoise=a.znoise,
                  L=a.L, P0=a.P0, rho1=1.0, rho2=2.0, U1=0.5, U2=-0.5, k=4 * np.pi,
                  t_end=a.tend, sigma_linear=sig, lattice=info,
                  cells_per_L=a.L * n)
    tot, bad = ic_common.report(d, meta)
    print("KH3D: %d x %d x %d, L/dx=%.2f, sigma_linear=%s" %
          (n, n, nz, a.L * n, "%.4f" % sig if sig else "skipped"))
    out = a.out or "kh3d_n%d_lz%g.l3d" % (n, lz)
    ic_common.write(out, d, meta, params, tot)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
