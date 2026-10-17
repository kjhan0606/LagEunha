#!/usr/bin/env python3
"""Analyse a 3D Noh snapshot against the exact solution.

  analyze.py --snap run/final.l3d --ic noh3d_n128_L6.l3d [-o noh3d]
  analyze.py --selftest

Measured inside the valid sphere r < L/2 - v0 t:
  R_meas     from the shocked mass: particles with P/rho^g above the geometric
             mean of the inflow and post-shock entropies; they started inside
             R + v0 t, so R = (3 M_s / 4 pi rho0)^(1/3) - v0 t
  rho_plat   volume-weighted mean density for 0.3 R_an < r < 0.9 R_an
             (exact 64; the centre r < 0.3 R carries the known wall-heating dip)
  rho_centre mean density for r < 0.15 R_an (reported, wall heating)
  L1_pre     volume-weighted mean |rho - (1+t/r)^2| / (1+t/r)^2 for
             1.3 R_an < r < r_valid
  L1_rho     volume-weighted mean |rho - rho_exact| over r < r_valid
  dE/E0      whole periodic box, (sum 1/2 m v^2 + m u)(t) vs IC

PASS (at 128^3 with L = 6, R = 14 dx):
  |R_meas/R_an - 1| <= 0.05
  |rho_plat/64 - 1| <= 0.20
  L1_pre <= 0.05
  |dE/E0| <= 1e-3
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
import noh_exact  # noqa: E402


def analyse(snap, ic=None, out="noh3d", plot=True):
    d, m = lag3d_io.read_l3d(snap)
    g, t = m["gamma"], m["time"]
    L = m["box"][1] - m["box"][0]
    p = dict(rho0=1.0, v0=1.0, p0=1e-6, center=[L / 2] * 3)
    E0 = None
    if ic:
        sc = ic_common.load_sidecar(ic)
        if sc:
            p.update({k: sc["params"][k] for k in p if k in sc["params"]})
        di, _ = lag3d_io.read_l3d(ic)
        ti = lag3d_io.totals(di)
        E0 = ti["Ekin"] + ti["Eint"]
    r, vr, _ = radial.radii(d, p["center"], m["box"], m["periodic"])
    rho = lag3d_io.density(d)
    P = (g - 1) * rho * d["u"]
    w = d["vol"] if np.all(d["vol"] > 0) else d["mass"] / rho
    rho_an, v_an, p_an, Ran = noh_exact.noh_profile(r, t, g, p["rho0"], p["v0"])
    rvalid = 0.5 * L - p["v0"] * t
    A = P / rho ** g
    A_in = p["p0"] / p["rho0"] ** g
    comp = (g + 1) / (g - 1)
    A_post = (g - 1) * p["rho0"] * comp ** 3 * 0.5 * p["v0"] ** 2 / (p["rho0"] * comp ** 3) ** g
    A_thr = np.sqrt(A_in * A_post)
    sh = (A > A_thr) & (r < rvalid)
    Ms = d["mass"][sh].sum()
    Rm = (3 * Ms / (4 * np.pi * p["rho0"])) ** (1 / 3) - p["v0"] * t

    def wmean(sel, q):
        return float(np.sum(w[sel] * q[sel]) / np.sum(w[sel])) if sel.any() else float("nan")
    plat = (r > 0.3 * Ran) & (r < 0.9 * Ran)
    rho_plat = wmean(plat, rho)
    rho_cen = wmean(r < 0.15 * Ran, rho)
    pre = (r > 1.3 * Ran) & (r < rvalid)
    L1_pre = wmean(pre, np.abs(rho - rho_an) / rho_an)
    val = r < rvalid
    L1 = wmean(val, np.abs(rho - rho_an))
    tot = lag3d_io.totals(d)
    E = tot["Ekin"] + tot["Eint"]
    dE = (E - E0) / E0 if E0 else float("nan")
    res = dict(snap=snap, n=round(m["np"] ** (1 / 3)), t=t, L=L, r_valid=rvalid,
               R_meas=float(Rm), R_an=float(Ran), R_ratio=float(Rm / Ran),
               rho_plat=rho_plat, rho_centre=rho_cen, L1_pre=L1_pre, L1_rho=L1,
               E=E, E0=E0, dE_over_E0=float(dE))
    crit = dict(R=abs(res["R_ratio"] - 1) <= 0.05,
                plateau=abs(rho_plat / comp ** 3 - 1) <= 0.20,
                pre=L1_pre <= 0.05,
                E=abs(dE) <= 1e-3)
    res["criteria"] = crit
    res["PASS"] = bool(all(crit.values()))
    with open(out + "_summary.json", "w") as f:
        json.dump(res, f, indent=1, default=float)
    if plot:
        plt = lag3d_io.get_plt()
        h = L / res["n"]
        edges = np.arange(0, rvalid + h, 0.5 * h)
        rc, rb, _ = radial.binned(r, rho, w, edges)
        _, vb, _ = radial.binned(r, vr, w, edges)
        _, pb, _ = radial.binned(r, P, w, edges)
        rr = np.linspace(1e-3, rvalid, 1000)
        ra, va, pa, _ = noh_exact.noh_profile(rr, t, g, p["rho0"], p["v0"])
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        for a_, q, qa, lab in ((ax[0], rb, ra, "rho"), (ax[1], vb, va, "v_r"), (ax[2], pb, pa, "P")):
            a_.plot(rc, q, "o", ms=3, label="N=%d^3" % res["n"])
            a_.plot(rr, qa, "k-", label="exact")
            a_.set_xlabel("r"); a_.set_ylabel(lab)
        ax[0].legend()
        ax[0].set_title("t=%.3f plateau=%.1f R/R_an=%.3f" % (t, rho_plat, res["R_ratio"]))
        lag3d_io.savefig(fig, out + "_profiles.png")
    print(json.dumps(res, indent=1, default=float))
    print("PASS" if res["PASS"] else "FAIL")
    return res


def selftest(tmp):
    """Exact solution checks + a synthetic snapshot on a 64^3, L=6 lattice."""
    import make_ic
    os.makedirs(tmp, exist_ok=True)
    rho, v, pp, rs = noh_exact.noh_profile(np.array([0.1, 1.0]), 2.0)
    ok = abs(rho[0] - 64) < 1e-12 and abs(pp[0] - 64 / 3) < 1e-12 and abs(rs - 2 / 3) < 1e-12 \
        and abs(rho[1] - 9) < 1e-12
    print("exact: rho_post=%.3f (64) P_post=%.4f (21.333) R(2)=%.4f (0.6667) rho(1,2)=%.3f (9)"
          % (rho[0], pp[0], rs, rho[1]))
    ic = os.path.join(tmp, "noh_ic_n64.l3d")
    make_ic.main(["--n", "64", "-o", ic])
    d, m = lag3d_io.read_l3d(ic)
    g, t = 5.0 / 3.0, 2.0
    r = np.sqrt((d["x"] - 3) ** 2 + (d["y"] - 3) ** 2 + (d["z"] - 3) ** 2)
    ra, va, pa, _ = noh_exact.noh_profile(r, t, g, p0=1e-6)
    d["rho"] = ra
    d["mass"] = ra * (6.0 / 64) ** 3
    d["u"] = np.where(pa > 1e-6, pa, 1e-6 * ra ** g) / ((g - 1) * ra)
    for k, c in (("vx", "x"), ("vy", "y"), ("vz", "z")):
        d[k] = va * (d[c] - 3) / np.maximum(r, 1e-300)
    snap = os.path.join(tmp, "noh_synth_n64.l3d")
    lag3d_io.write_l3d(snap, d, time=t, gamma=g, box=m["box"])
    res = analyse(snap, ic, out=os.path.join(tmp, "noh3d_selftest"))
    ok = ok and res["criteria"]["R"] and res["criteria"]["plateau"] and res["criteria"]["pre"]
    print("(synthetic energy is not a hydro result; dE/E0 reported only)")
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snap")
    ap.add_argument("--ic", default=None)
    ap.add_argument("-o", "--out", default="noh3d")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tmp", default="/tmp/lag3d_selftest/noh")
    a = ap.parse_args()
    if a.selftest:
        sys.path.insert(0, HERE)
        sys.exit(selftest(a.tmp))
    if not a.snap:
        ap.error("--snap required")
    r = analyse(a.snap, a.ic, a.out)
    sys.exit(0 if r["PASS"] else 2)
