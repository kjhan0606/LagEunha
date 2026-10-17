#!/usr/bin/env python3
"""3D linear sound wave IC (periodic unit cube).

rho0 = 1, c0 = 1 (P0 = rho0 c0^2/gamma), gamma = 5/3, right-moving adiabatic
wave with amplitude A (default 1e-6):
    rho = rho0 (1 + A s),  v = v_boost + c0 A s k_hat,  P = P0 (rho/rho0)^gamma,
    s = sin(k . x),  k = 2 pi (1,0,0) [--dir x] or 2 pi (1,1,1) [--dir diag].
One period is T = lambda / c0 (1 for x, 1/sqrt(3) for diag); with a boost the
pattern is advected by v_boost as well, so the analytic solution at time t is
the IC profile shifted by (c0 k_hat + v_boost) t.

--mode mass     (default) cell-centred cubic lattice, equal volume, mass
                m_i = rho(x_i) dx^3. This is how the 2D KH/Gresho ICs are
                built (util.c: mass = rho*meanvol). Density is exact at t=0.
--mode displace equal-mass particles displaced by the linear Lagrangian
                displacement xi = +(A/k) cos(k.x0) k_hat.
--lattice glass only for A >= 1e-3: glass volume noise (~few %) swamps a
                1e-6 wave, because rho = m/V_voronoi in the code.
"""
import argparse
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lattice  # noqa: E402
import ic_common  # noqa: E402


def wave_setup(direction):
    if direction == "x":
        kv = 2 * np.pi * np.array([1.0, 0.0, 0.0])
    elif direction == "diag":
        kv = 2 * np.pi * np.array([1.0, 1.0, 1.0])
    else:
        raise ValueError(direction)
    return kv, kv / np.linalg.norm(kv), 2 * np.pi / np.linalg.norm(kv)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=32, help="particles per side")
    ap.add_argument("--amp", type=float, default=1e-6)
    ap.add_argument("--boost", type=float, default=0.0,
                    help="uniform bulk velocity along +x (Galilean test)")
    ap.add_argument("--dir", default="x", choices=["x", "diag"])
    ap.add_argument("--mode", default="mass", choices=["mass", "displace"])
    ap.add_argument("--lattice", default="cubic", choices=["cubic", "glass"])
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("-o", "--out", default=None)
    a = ap.parse_args(argv)
    if a.lattice == "glass" and a.amp < 1e-3:
        print("WARNING: glass + A<1e-3: glass volume noise will dominate the wave")
    n, g = a.n, a.gamma
    rho0, c0 = 1.0, 1.0
    P0 = rho0 * c0 ** 2 / g
    box = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    kv, khat, lam = wave_setup(a.dir)
    x, y, z, info = lattice.positions(a.lattice, n, n, n, box)
    V = 1.0 / n ** 3
    X = np.stack([x, y, z], axis=1)
    if a.mode == "mass":
        s = np.sin(X @ kv)
        rho = rho0 * (1 + a.amp * s)
        mass = rho * V
    else:
        X0 = X.copy()
        s0 = X0 @ kv
        X = X0 + (a.amp / np.linalg.norm(kv)) * np.cos(s0)[:, None] * khat[None, :]
        X = np.mod(X, 1.0)
        s = np.sin(s0)                       # Lagrangian phase
        rho = rho0 * (1 + a.amp * s)         # rho0/(1 + xi') = rho0 (1 + A sin) to 1st order
        mass = np.full(len(x), rho0 * V)
    P = P0 * (rho / rho0) ** g
    u = P / ((g - 1) * rho)
    vel = a.amp * c0 * s[:, None] * khat[None, :]
    vel[:, 0] += a.boost
    d = dict(x=X[:, 0], y=X[:, 1], z=X[:, 2], vx=vel[:, 0], vy=vel[:, 1],
             vz=vel[:, 2], mass=mass, u=u, rho=rho, vol=np.zeros(len(x)),
             id=np.arange(len(x), dtype=np.int64))
    meta = dict(gamma=g, box=box, periodic=(1, 1, 1))
    T = lam / c0
    params = dict(test="SoundWave3D", n=n, amp=a.amp, boost=a.boost, dir=a.dir,
                  mode=a.mode, lattice=info, rho0=rho0, c0=c0, P0=P0,
                  k=list(kv), period=T, t_end=T)
    tot, bad = ic_common.report(d, meta)
    print("wave: lambda=%.6f  period T=%.6f  (t_end=T)  kL/2pi=%s" % (lam, T, kv / (2 * np.pi)))
    out = a.out or "sw3d_n%d%s_%s.l3d" % (n, "_boost" if a.boost else "", a.dir)
    ic_common.write(out, d, meta, params, tot)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
