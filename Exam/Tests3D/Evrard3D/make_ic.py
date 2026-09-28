#!/usr/bin/env python3
"""Evrard (1988) adiabatic collapse IC.

Standard setup (Evrard 1988, MNRAS 235, 911; used in Springel 2005/2010 and
Hopkins 2015): gamma = 5/3, G = M = R = 1,
    rho(r) = M / (2 pi R^2 r)  for r < R  (so M(<r) = r^2),
    u = 0.05 G M / R,  v = 0.
E_pot(0) = -2/3 (continuum), E_th(0) = 0.05, E_tot(0) = -0.6167.

Particles: a cell-centred cubic lattice (or glass) in [-1,1]^3, cut at r<1,
then stretched radially r -> r^(3/2), which maps uniform density onto
rho ~ 1/r with equal particle masses m = M/N.

Gravity: G = 1 and a Plummer softening eps (default 0.01) are written into
the file header. The IC's E_pot is a direct softened sum when N <= 60000,
otherwise the spherical estimate -sum G M(<r_i) m / r_i.

Outer boundary: vacuum by default (particle-code convention). A Voronoi/
Laguerre code needs a bounded tessellation; --background uniform adds
low-density gas (rho_bg, u = 0.05) on a coarse lattice out to |x| < B.
Which one the 3D driver needs is an open decision (see README).
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lattice  # noqa: E402
import ic_common  # noqa: E402


def potential_energy(X, m, eps, G=1.0, direct_max=60000):
    N = len(m)
    if N <= direct_max:
        W = 0.0
        chunk = max(1, int(2e7 // N))
        for i0 in range(0, N, chunk):
            Xi = X[i0:i0 + chunk]
            d2 = ((Xi[:, None, :] - X[None, :, :]) ** 2).sum(-1) + eps * eps
            inv = 1.0 / np.sqrt(d2)
            for k in range(len(Xi)):
                inv[k, i0 + k] = 0.0
            W -= G * 0.5 * (m[i0:i0 + chunk, None] * m[None, :] * inv).sum()
        return float(W), "direct softened sum"
    r = np.linalg.norm(X, axis=1)
    o = np.argsort(r)
    Menc = np.cumsum(m[o]) - m[o]
    W = -G * np.sum(Menc * m[o] / np.maximum(r[o], 1e-12))
    return float(W), "spherical estimate (unsoftened)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=64,
                    help="lattice points per side of [-1,1]^3 before the cut (N ~ 0.52 n^3)")
    ap.add_argument("--lattice", default="cubic", choices=["cubic", "glass"])
    ap.add_argument("--eps", type=float, default=0.01, help="Plummer softening")
    ap.add_argument("--u0", type=float, default=0.05)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--background", default="none", choices=["none", "uniform"])
    ap.add_argument("--rho-bg", type=float, default=1e-4)
    ap.add_argument("--nbg", type=int, default=16, help="background lattice per side of [-B,B]^3")
    ap.add_argument("--B", type=float, default=2.0, help="half box size")
    ap.add_argument("--tend", type=float, default=3.0)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args(argv)
    g = a.gamma
    x, y, z, info = lattice.positions(a.lattice, a.n, a.n, a.n, (-1, 1, -1, 1, -1, 1))
    X = np.stack([x, y, z], axis=1)
    r = np.linalg.norm(X, axis=1)
    keep = r < 1.0
    X, r = X[keep], r[keep]
    X = X * np.sqrt(r)[:, None]            # |x| -> |x|^(3/2)
    N = len(X)
    m = np.full(N, 1.0 / N)
    u = np.full(N, a.u0)
    W, how = potential_energy(X, m, a.eps)
    nbg = 0
    if a.background == "uniform":
        xb, yb, zb, _ = lattice.positions("cubic", a.nbg, a.nbg, a.nbg, (-a.B, a.B) * 3)
        Xb = np.stack([xb, yb, zb], axis=1)
        Xb = Xb[np.linalg.norm(Xb, axis=1) > 1.05]
        mb = a.rho_bg * (2 * a.B / a.nbg) ** 3
        nbg = len(Xb)
        X = np.vstack([X, Xb])
        m = np.concatenate([m, np.full(nbg, mb)])
        u = np.concatenate([u, np.full(nbg, a.u0)])
        W, how = potential_energy(X, m, a.eps)
    Ntot = len(m)
    rr = np.linalg.norm(X, axis=1)
    rho = np.where(rr < 1, 1.0 / (2 * np.pi * np.maximum(rr, 1e-3)), a.rho_bg)
    B = a.B * 1.0000001
    d = dict(x=X[:, 0], y=X[:, 1], z=X[:, 2], vx=np.zeros(Ntot), vy=np.zeros(Ntot),
             vz=np.zeros(Ntot), mass=m, u=u, rho=rho, vol=np.zeros(Ntot),
             id=np.arange(Ntot, dtype=np.int64))
    meta = dict(gamma=g, box=(-B, B, -B, B, -B, B), periodic=(0, 0, 0), G=1.0,
                softening=a.eps)
    params = dict(test="Evrard3D", n=a.n, N_sphere=N, N_background=nbg,
                  lattice=info, eps=a.eps, u0=a.u0, G=1.0, M=1.0, R=1.0,
                  background=a.background, rho_bg=a.rho_bg if nbg else 0.0,
                  t_end=a.tend, Epot_method=how, Epot=W,
                  Epot_continuum=-2.0 / 3.0)
    tot, bad = ic_common.report(d, meta, pot_energy=W)
    print("Evrard: N_sphere=%d  N_bg=%d  eps=%.4g  Epot=%.5f (%s; continuum -0.66667)" %
          (N, nbg, a.eps, W, how))
    out = a.out or "evrard3d_n%d_%s%s.l3d" % (a.n, a.lattice, "_bg" if nbg else "")
    ic_common.write(out, d, meta, params, tot)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
