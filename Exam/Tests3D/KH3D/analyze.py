#!/usr/bin/env python3
"""Linear growth of the seeded KH mode in 3D vs the linear theory.

  analyze.py --snaps run/snap_*.l3d [--snaps2 run256/snap_*.l3d] [--ic kh3d_n128_lz0.25.l3d] [-o kh3d]
  analyze.py --selftest

Mode amplitude, McNally et al. (2012) Eq. 6-9 with cell volumes V_i
(= m_i/rho_i if the snapshot has no vol[]), summed over all z:
  s = sum vy V sin(4 pi x) E(y),  c = sum vy V cos(4 pi x) E(y),  dd = sum V E(y)
  E(y) = exp(-4 pi |y - 1/4|) for y < 1/2, exp(-4 pi |(1 - y) - 1/4|) otherwise
  M = 2 sqrt((s/dd)^2 + (c/dd)^2)          (M = 0.01 at t = 0 for the default seed)
sigma_fit is the least-squares slope of ln M(t) over the snapshots with
1.5 M(0) <= M <= 0.06 (after the initial transient, before saturation),
at least 4 points. sigma_lin comes from the IC sidecar or is recomputed by
common/kh_linear.py (compressible eigen-solve of the exact McNally profile;
2.83 for the default setup).

PASS:
  K1  |sigma_fit / sigma_lin - 1| <= 0.15 at n >= 128 per unit length
  K2  (with --snaps2 at higher resolution) the error |sigma_fit/sigma_lin - 1|
      does not grow with resolution
Writes <o>_summary.json and <o>_growth.png.
"""
import argparse
import glob
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lag3d_io  # noqa: E402
import ic_common  # noqa: E402
import kh_linear  # noqa: E402


def mode_amplitude(d):
    rho = lag3d_io.density(d)
    V = d["vol"] if np.all(d["vol"] > 0) else d["mass"] / rho
    y = np.mod(d["y"], 1.0)
    E = np.where(y < 0.5, np.exp(-4 * np.pi * np.abs(y - 0.25)),
                 np.exp(-4 * np.pi * np.abs((1 - y) - 0.25)))
    s = np.sum(d["vy"] * V * np.sin(4 * np.pi * d["x"]) * E)
    c = np.sum(d["vy"] * V * np.cos(4 * np.pi * d["x"]) * E)
    dd = np.sum(V * E)
    return 2 * np.hypot(s / dd, c / dd)


def series(paths):
    ts, Ms, n = [], [], None
    for p in sorted(paths):
        d, m = lag3d_io.read_l3d(p)
        ts.append(m["time"]); Ms.append(mode_amplitude(d))
        n = round((m["np"] / ((m["box"][5] - m["box"][4]))) ** (1 / 3))
    o = np.argsort(ts)
    return np.array(ts)[o], np.array(Ms)[o], n


def fit(t, M, lo=1.5, hi=0.06):
    sel = (M >= lo * M[0]) & (M <= hi)
    if sel.sum() < 4:
        return float("nan"), sel
    A = np.polyfit(t[sel], np.log(M[sel]), 1)
    return float(A[0]), sel


def analyse(snaps, snaps2=None, ic=None, out="kh3d", plot=True, sigma_lin=None):
    if sigma_lin is None and ic:
        sc = ic_common.load_sidecar(ic)
        if sc and sc["params"].get("sigma_linear"):
            sigma_lin = sc["params"]["sigma_linear"]
    if sigma_lin is None:
        sigma_lin = kh_linear.growth_rate(n=512)[0]
    t, M, n = series(snaps)
    sf, sel = fit(t, M)
    res = dict(n=n, sigma_lin=sigma_lin, sigma_fit=sf, ratio=sf / sigma_lin,
               M0=float(M[0]), M_max=float(M.max()), t=t.tolist(), M=M.tolist(),
               n_fit_points=int(sel.sum()))
    err = abs(sf / sigma_lin - 1)
    crit = dict(K1=bool(n >= 128 and err <= 0.15))
    if snaps2:
        t2, M2, n2 = series(snaps2)
        sf2, sel2 = fit(t2, M2)
        res["run2"] = dict(n=n2, sigma_fit=sf2, ratio=sf2 / sigma_lin, t=t2.tolist(), M=M2.tolist())
        crit["K2"] = bool(abs(sf2 / sigma_lin - 1) <= err + 0.01)
    res["criteria"] = crit
    res["PASS"] = bool(all(crit.values()))
    with open(out + "_summary.json", "w") as f:
        json.dump(res, f, indent=1)
    if plot:
        plt = lag3d_io.get_plt()
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.semilogy(t, M, "o-", ms=3, label="n=%d" % n)
        if snaps2:
            ax.semilogy(res["run2"]["t"], res["run2"]["M"], "s-", ms=3, label="n=%d" % res["run2"]["n"])
        if np.isfinite(sf):
            t0 = t[sel][0]
            ax.semilogy(t, M[sel][0] * np.exp(sigma_lin * (t - t0)), "k--",
                        label="linear theory sigma=%.2f" % sigma_lin)
        ax.set_xlabel("t"); ax.set_ylabel("mode amplitude M")
        ax.set_title("sigma_fit=%.3f  ratio=%.3f" % (sf, sf / sigma_lin))
        ax.legend(fontsize=8)
        lag3d_io.savefig(fig, out + "_growth.png")
    print(json.dumps({k: v for k, v in res.items() if k not in ("t", "M", "run2")}, indent=1))
    print("PASS" if res["PASS"] else "FAIL")
    return res


def selftest(tmp):
    import make_ic
    os.makedirs(tmp, exist_ok=True)
    ic = os.path.join(tmp, "kh_ic_n128.l3d")
    make_ic.main(["--n", "128", "--lz", str(4.0 / 128), "-o", ic])
    d, m = lag3d_io.read_l3d(ic)
    sc = ic_common.load_sidecar(ic)
    sl = sc["params"]["sigma_linear"]
    M0 = mode_amplitude(d)
    print("IC mode amplitude M0 = %.5f (expect 0.01)" % M0)
    ok = abs(M0 - 0.01) < 2e-4 and abs(sl - 2.83) < 0.03

    def make_series(tag, sigma):
        paths = []
        for i, tt in enumerate(np.arange(0, 1.51, 0.05)):
            dd = dict(d)
            amp = min(0.01 * np.exp(sigma * tt), 0.15)  # exponential, then saturated
            dd["vy"] = d["vy"] * amp / 0.01
            p = os.path.join(tmp, "%s_%03d.l3d" % (tag, i))
            lag3d_io.write_l3d(p, dd, time=tt, gamma=m["gamma"], box=m["box"])
            paths.append(p)
        return paths
    good = analyse(make_series("good", 0.95 * sl), ic=ic, out=os.path.join(tmp, "good"))
    ok &= good["PASS"] and abs(good["ratio"] - 0.95) < 0.03
    bad = analyse(make_series("bad", 0.7 * sl), ic=ic, out=os.path.join(tmp, "bad"), plot=False)
    ok &= not bad["PASS"]
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snaps", nargs="*", default=[])
    ap.add_argument("--snaps2", nargs="*", default=None)
    ap.add_argument("--ic", default=None)
    ap.add_argument("-o", "--out", default="kh3d")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tmp", default="/tmp/lag3d_selftest/kh")
    a = ap.parse_args()
    if a.selftest:
        sys.path.insert(0, HERE)
        sys.exit(selftest(a.tmp))
    snaps = [p for s in a.snaps for p in (glob.glob(s) or [s])]
    if not snaps:
        ap.error("--snaps required")
    snaps2 = [p for s in a.snaps2 for p in (glob.glob(s) or [s])] if a.snaps2 else None
    r = analyse(snaps, snaps2, a.ic, a.out)
    sys.exit(0 if r["PASS"] else 2)
