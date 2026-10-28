#!/usr/bin/env python3
"""Analyse a 3D Sedov-Taylor snapshot against the exact solution.

  analyze.py --snap run/final.l3d --ic sedov3d_n64_cubic.l3d [--snap2 run128/final.l3d] [-o sedov3d]
  analyze.py --selftest

Measures (at the snapshot time t, centre from the IC sidecar, default box centre):
  R_meas   swept-mass radius: particles with P/rho^g above the geometric mean of
           the ambient and the exact post-shock entropy are "shocked";
           R_meas = (3 M_shocked / 4 pi rho0)^(1/3)  (no bin-width bias)
  R_peak   radius of the binned density peak (reported only; biased ~dx/2 low)
  R_an     xi0 (E t^2/rho0)^(1/5), xi0 from common/sedov_exact.py (1.15167 for 5/3)
  rho_peak binned peak density (exact post-shock value 4 for gamma=5/3)
  L1_rho   volume-weighted mean |rho_i - rho_exact(r_i)| over r < 1.5 R_an
  dE/E0    (sum 1/2 m v^2 + m u)(t) vs the IC total (the IC file, or its sidecar)
  aniso    R_meas along the axes / along the body diagonals (cones of 10 deg)

PASS (at 64^3; tighter numbers apply automatically at >= 128^3):
  |R_meas/R_an - 1| <= 0.03 (0.02 at 128^3)
  |dE/E0| <= 1e-3
  rho_peak >= 2.5 (>= 3.0 at 128^3)
  |aniso - 1| <= 0.05 (lattice imprint)
  if --snap2 (higher resolution) is given: L1_rho(snap2) < L1_rho(snap)
Writes <o>_summary.json and <o>_profiles.png.
"""
import argparse
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lag3d_io  # noqa: E402
import ic_common  # noqa: E402
import radial  # noqa: E402
import sedov_exact  # noqa: E402


def measure(snap, E0=None, params=None):
    d, m = lag3d_io.read_l3d(snap)
    g, t = m["gamma"], m["time"]
    p = dict(E=1.0, rho0=1.0, center=[0.5, 0.5, 0.5], pamb=1e-5)
    if params:
        p.update({k: params[k] for k in p if k in params})
    r, vr, dx = radial.radii(d, p["center"], m["box"], m["periodic"])
    rho = lag3d_io.density(d)
    P = (g - 1) * rho * d["u"]
    w = d["vol"] if np.all(d["vol"] > 0) else d["mass"] / rho
    n = round(m["np"] ** (1 / 3))
    h = (m["box"][1] - m["box"][0]) / n
    edges = np.arange(0, 0.5 * (m["box"][1] - m["box"][0]), h)
    rc, rb, _ = radial.binned(r, rho, w, edges)
    _, vb, _ = radial.binned(r, vr, w, edges)
    _, pb, _ = radial.binned(r, P, w, edges)
    Rpk, rpk = radial.peak_radius(rc, rb)
    rho_an, v_an, p_an, Ran = sedov_exact.sedov_profile(r, t, g, p["E"], p["rho0"], p["pamb"])
    # Shock radius from the swept-up mass: a particle is shocked if its
    # entropy A = P/rho^g exceeds the geometric mean of the ambient value
    # and the exact post-shock value. The ambient density is uniform, so the
    # shocked mass fixes the radius: M_s = 4/3 pi R^3 rho0. Unlike the
    # density-peak radius this has no bin-width bias.
    A = P / rho ** g
    A_amb = p["pamb"] / p["rho0"] ** g
    rho2, _, p2, _ = sedov_exact.sedov_profile(np.array([Ran * (1 - 1e-9)]), t, g, p["E"], p["rho0"])
    A_thr = np.sqrt(A_amb * p2[0] / rho2[0] ** g)
    Ms = d["mass"][A > A_thr].sum()
    Rm = (3 * Ms / (4 * np.pi * p["rho0"])) ** (1 / 3)
    sel = r < 1.5 * Ran
    L1 = float(np.sum(w[sel] * np.abs(rho[sel] - rho_an[sel])) / np.sum(w[sel]))
    tot = lag3d_io.totals(d)
    E = tot["Ekin"] + tot["Eint"]
    dE = (E - E0) / E0 if E0 else float("nan")
    # anisotropy: axes vs body diagonals
    un = dx / np.maximum(r, 1e-300)[:, None]
    cmax = np.abs(un).max(axis=1)                       # cos to nearest axis
    cdiag = np.abs(un).sum(axis=1) / np.sqrt(3)         # cos to nearest diagonal
    cc = np.cos(np.radians(10))
    Rax = radial.peak_radius(*radial.binned(r[cmax > cc], rho[cmax > cc], w[cmax > cc], edges)[:2])[0]
    Rdg = radial.peak_radius(*radial.binned(r[cdiag > cc], rho[cdiag > cc], w[cdiag > cc], edges)[:2])[0]
    return dict(snap=snap, n=n, t=t, R_meas=float(Rm), R_peak=float(Rpk), R_an=float(Ran),
                R_ratio=float(Rm / Ran), rho_peak=float(rpk), L1_rho=L1,
                E=E, E0=E0, dE_over_E0=float(dE), aniso=float(Rax / Rdg),
                _prof=(rc, rb, vb, pb), _gamma=g, _p=p)


def criteria(res):
    hi = res["n"] >= 128
    c = dict(R=abs(res["R_ratio"] - 1) <= (0.02 if hi else 0.03),
             E=abs(res["dE_over_E0"]) <= 1e-3,
             peak=res["rho_peak"] >= (3.0 if hi else 2.5),
             aniso=abs(res["aniso"] - 1) <= 0.05)
    return c


def analyse(snap, ic=None, snap2=None, out="sedov3d", plot=True):
    params, E0 = None, None
    if ic:
        sc = ic_common.load_sidecar(ic)
        if sc:
            params = sc["params"]
        di, _ = lag3d_io.read_l3d(ic)
        ti = lag3d_io.totals(di)
        E0 = ti["Ekin"] + ti["Eint"]
    res = measure(snap, E0, params)
    crit = criteria(res)
    out_d = {k: v for k, v in res.items() if not k.startswith("_")}
    out_d["criteria"] = crit
    ok = all(crit.values())
    if snap2:
        r2 = measure(snap2, E0 * 1.0 if E0 else None, params)
        out_d["snap2"] = {k: v for k, v in r2.items() if not k.startswith("_")}
        out_d["criteria"]["L1_converges"] = bool(r2["L1_rho"] < res["L1_rho"])
        ok = ok and out_d["criteria"]["L1_converges"]
    out_d["PASS"] = bool(ok)
    with open(out + "_summary.json", "w") as f:
        json.dump(out_d, f, indent=1, default=float)
    if plot:
        plt = lag3d_io.get_plt()
        rc, rb, vb, pb = res["_prof"]
        rr = np.linspace(1e-4, rc.max(), 800)
        ra, va, pa, Ran = sedov_exact.sedov_profile(rr, res["t"], res["_gamma"],
                                                    res["_p"]["E"], res["_p"]["rho0"])
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        for a_, q, qa, lab in ((ax[0], rb, ra, "rho"), (ax[1], vb, va, "v_r"), (ax[2], pb, pa, "P")):
            a_.plot(rc, q, "o", ms=3, label="N=%d^3" % res["n"])
            a_.plot(rr, qa, "k-", label="exact")
            a_.set_xlabel("r"); a_.set_ylabel(lab)
        ax[0].legend()
        ax[0].set_title("t=%.4f R/R_an=%.3f dE/E0=%.1e" % (res["t"], res["R_ratio"], res["dE_over_E0"]))
        lag3d_io.savefig(fig, out + "_profiles.png")
    print(json.dumps(out_d, indent=1, default=float))
    print("PASS" if ok else "FAIL")
    return out_d


def selftest(tmp):
    """(1) exact-solution energy integral == E; (2) a synthetic snapshot built
    from the exact profile on a 64^3 lattice is recognised (R within 3%)."""
    from scipy.integrate import quad
    import make_ic
    os.makedirs(tmp, exist_ok=True)
    g, t = 5.0 / 3.0, 0.05
    xi0 = sedov_exact.sedov_xi0(g)
    R = sedov_exact.shock_radius(t, g)

    def e_dens(r):
        rho, v, p, _ = sedov_exact.sedov_profile(np.array([r]), t, g)
        return (0.5 * rho[0] * v[0] ** 2 + p[0] / (g - 1)) * 4 * np.pi * r * r
    Eint = quad(e_dens, 1e-6, R, limit=400, points=[0.9 * R, 0.99 * R])[0]
    print("selftest: xi0=%.5f (expect 1.1517), energy integral=%.6f (expect 1)" % (xi0, Eint))
    ok = abs(xi0 - 1.1517) < 5e-4 and abs(Eint - 1) < 2e-3
    ic = os.path.join(tmp, "sedov_ic_n64.l3d")
    make_ic.main(["--n", "64", "-o", ic])
    d, m = lag3d_io.read_l3d(ic)
    r = np.sqrt((d["x"] - .5) ** 2 + (d["y"] - .5) ** 2 + (d["z"] - .5) ** 2)
    rho, v, p, _ = sedov_exact.sedov_profile(r, t, g, p0=1e-5)
    d["rho"] = rho
    d["mass"] = rho / 64 ** 3          # mass consistent with volume
    d["u"] = p / ((g - 1) * rho)
    for k, c in (("vx", "x"), ("vy", "y"), ("vz", "z")):
        d[k] = v * (d[c] - .5) / np.maximum(r, 1e-300)
    snap = os.path.join(tmp, "sedov_synth_n64.l3d")
    lag3d_io.write_l3d(snap, d, time=t, gamma=g, box=m["box"])
    res = analyse(snap, ic, out=os.path.join(tmp, "sedov3d_selftest"))
    ok = ok and abs(res["R_ratio"] - 1) < 0.03 and res["rho_peak"] > 2.5 and abs(res["aniso"] - 1) < 0.05
    print("(synthetic snapshot energy is the lattice sampling of the exact solution; "
          "dE/E0=%.2e is reported, not tested)" % res["dE_over_E0"])
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snap")
    ap.add_argument("--snap2", default=None)
    ap.add_argument("--ic", default=None)
    ap.add_argument("-o", "--out", default="sedov3d")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tmp", default="/tmp/lag3d_selftest/sedov")
    a = ap.parse_args()
    if a.selftest:
        sys.path.insert(0, HERE)
        sys.exit(selftest(a.tmp))
    if not a.snap:
        ap.error("--snap required")
    r = analyse(a.snap, a.ic, a.snap2, a.out)
    sys.exit(0 if r["PASS"] else 2)
