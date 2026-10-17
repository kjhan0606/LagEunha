#!/usr/bin/env python3
"""Analyse 3D sound-wave runs: L1 error vs analytic, convergence order,
Galilean invariance.

  analyze.py --snaps run32/final.l3d run64/final.l3d run128/final.l3d \
             [--boosted run64b/final.l3d] [--ic ic_n32_x.json] [-o sw3d]
  analyze.py --selftest

Each snapshot is compared with the IC profile advected by (c0 k_hat +
v_boost) t, evaluated at the snapshot particle positions (Lagrangian
positions, Eulerian analytic field). The wave parameters come from the
snapshot's JSON sidecar if present, else from --ic, else defaults
(A=1e-6, x direction, no boost).

Outputs <o>_summary.json, <o>_convergence.png, <o>_profile.png.

PASS (all must hold):
  P1  convergence order of L1(rho) between the two finest runs >= 1.7
  P2  L1(rho)/(A rho0) <= 1e-2 at 64^3 (or the finest run if 64^3 is absent)
  P3  boosted run: L1_boost / L1_rest (same N) <= 1.25
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

DEF = dict(amp=1e-6, dir="x", boost=0.0, rho0=1.0, c0=1.0)


def wave_params(path, icjson=None):
    p = dict(DEF)
    sc = ic_common.load_sidecar(path)
    if sc is None and icjson:
        with open(icjson) as f:
            sc = json.load(f)
    if sc:
        p.update({k: v for k, v in sc["params"].items() if k in DEF})
    return p


def errors(d, t, p, gamma):
    kv = 2 * np.pi * (np.array([1.0, 0, 0]) if p["dir"] == "x" else np.ones(3))
    khat = kv / np.linalg.norm(kv)
    X = np.stack([d["x"], d["y"], d["z"]], axis=1)
    shift = p["c0"] * khat * t + np.array([p["boost"], 0, 0]) * t
    s = np.sin((X - shift) @ kv)
    rho_an = p["rho0"] * (1 + p["amp"] * s)
    vpar_an = p["c0"] * p["amp"] * s
    rho = lag3d_io.density(d)
    vpar = (np.stack([d["vx"] - p["boost"], d["vy"], d["vz"]], axis=1) @ khat)
    w = d["vol"] if np.all(d["vol"] > 0) else d["mass"] / rho
    w = w / w.sum()
    L1r = float(np.sum(w * np.abs(rho - rho_an)))
    L1v = float(np.sum(w * np.abs(vpar - vpar_an)))
    # fitted amplitude and phase of the wave in the snapshot
    c = np.cos((X - shift) @ kv)
    A_s = 2 * np.sum(w * (rho / p["rho0"] - 1) * s)
    A_c = 2 * np.sum(w * (rho / p["rho0"] - 1) * c)
    return dict(L1_rho=L1r, L1_v=L1v, L1_rho_over_A=L1r / (p["amp"] * p["rho0"]),
                amp_ratio=float(np.hypot(A_s, A_c) / p["amp"]),
                phase_err=float(np.arctan2(A_c, A_s)), rho=rho, rho_an=rho_an,
                phase=np.mod((X - shift) @ kv, 2 * np.pi))


def analyse(snaps, boosted=None, icjson=None, out="sw3d", plot=True):
    rows = []
    for path in snaps:
        d, m = lag3d_io.read_l3d(path)
        p = wave_params(path, icjson)
        e = errors(d, m["time"], p, m["gamma"])
        n = round(m["np"] ** (1 / 3))
        rows.append(dict(path=path, n=n, t=m["time"], **{k: v for k, v in e.items()
                                                       if not isinstance(v, np.ndarray)}))
        rows[-1]["_e"] = e
    rows.sort(key=lambda r: r["n"])
    ns = np.array([r["n"] for r in rows], float)
    L1 = np.array([r["L1_rho"] for r in rows])
    order_pairs = [float(-np.log(L1[i + 1] / L1[i]) / np.log(ns[i + 1] / ns[i]))
                   for i in range(len(rows) - 1)]
    order_fit = float(-np.polyfit(np.log(ns), np.log(L1), 1)[0]) if len(rows) > 1 else float("nan")
    ref = [r for r in rows if r["n"] == 64] or rows[-1:]
    res = dict(runs=[{k: v for k, v in r.items() if k != "_e"} for r in rows],
               order_pairs=order_pairs, order_fit=order_fit)
    P1 = bool(order_pairs and order_pairs[-1] >= 1.7)
    P2 = bool(ref[0]["L1_rho_over_A"] <= 1e-2)
    res["P1_order"] = P1
    res["P2_L1_64"] = P2
    passed = P1 and P2
    if boosted:
        db, mb = lag3d_io.read_l3d(boosted)
        pb = wave_params(boosted, icjson)
        if pb["boost"] == 0:
            print("WARNING: boosted snapshot has boost=0 in its sidecar")
        eb = errors(db, mb["time"], pb, mb["gamma"])
        nb = round(mb["np"] ** (1 / 3))
        rest = [r for r in rows if r["n"] == nb]
        ratio = eb["L1_rho"] / rest[0]["L1_rho"] if rest else float("nan")
        res["boosted"] = dict(path=boosted, n=nb, boost=pb["boost"], L1_rho=eb["L1_rho"],
                              L1_v=eb["L1_v"], ratio_to_rest=ratio)
        res["P3_galilean"] = bool(ratio <= 1.25)
        passed = passed and res["P3_galilean"]
    res["PASS"] = bool(passed)
    with open(out + "_summary.json", "w") as f:
        json.dump(res, f, indent=1)
    if plot:
        plt = lag3d_io.get_plt()
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.loglog(ns, L1, "o-", label="L1(rho)")
        ax.loglog(ns, L1[0] * (ns / ns[0]) ** -2, "k--", label="N^-2")
        if boosted:
            ax.loglog([res["boosted"]["n"]], [res["boosted"]["L1_rho"]], "rs", label="boosted")
        ax.set_xlabel("N per side"); ax.set_ylabel("L1 density error"); ax.legend()
        ax.set_title("3D sound wave, order=%.2f" % (order_pairs[-1] if order_pairs else np.nan))
        lag3d_io.savefig(fig, out + "_convergence.png")
        fig, ax = plt.subplots(figsize=(6, 4))
        for r in rows:
            e = r["_e"]
            idx = np.argsort(e["phase"])[:: max(1, len(e["phase"]) // 4000)]
            ax.plot(e["phase"][idx], e["rho"][idx] - 1, ".", ms=1, label="N=%d" % r["n"])
        ph = np.linspace(0, 2 * np.pi, 200)
        A = wave_params(rows[-1]["path"], icjson)["amp"]
        ax.plot(ph, A * np.sin(ph), "k-", label="analytic")
        ax.set_xlabel("phase k.(x - c t)"); ax.set_ylabel("rho - rho0"); ax.legend(markerscale=6)
        lag3d_io.savefig(fig, out + "_profile.png")
    print(json.dumps({k: v for k, v in res.items() if k != "runs"}, indent=1))
    for r in res["runs"]:
        print("  N=%4d t=%.4f L1(rho)=%.3e L1/A=%.3e amp_ratio=%.4f" %
              (r["n"], r["t"], r["L1_rho"], r["L1_rho_over_A"], r["amp_ratio"]))
    print("PASS" if passed else "FAIL")
    return res


def selftest(tmp):
    """Synthetic snapshots: analytic wave at t=T plus an error that scales as
    N^-2 (and a boosted copy); the script must recover order 2 and pass."""
    import make_ic
    os.makedirs(tmp, exist_ok=True)
    snaps = []
    for n in (16, 32, 64):
        f = os.path.join(tmp, "sw_n%d.l3d" % n)
        make_ic.main(["--n", str(n), "-o", f])
        d, m = lag3d_io.read_l3d(f)
        X = np.stack([d["x"], d["y"], d["z"]], axis=1)
        # at t = T the exact profile equals the IC; add a 2nd-order phase error
        s = np.sin(2 * np.pi * X[:, 0] - 0.5 * (16.0 / n) ** 2 * 0.1)
        d["rho"] = 1 + 1e-6 * s
        lag3d_io.write_l3d(f, d, time=1.0, gamma=m["gamma"], box=m["box"])
        snaps.append(f)
    fb = os.path.join(tmp, "sw_n32_boost.l3d")
    make_ic.main(["--n", "32", "--boost", "1.0", "-o", fb])
    d, m = lag3d_io.read_l3d(fb)
    X = np.stack([d["x"], d["y"], d["z"]], axis=1)
    # boosted: after t=T with v_b=1 the pattern moved by 2 box lengths
    d["x"] = np.mod(d["x"] + 1.0 * 1.0, 1.0)
    s = np.sin(2 * np.pi * (d["x"] - 2.0) - 0.5 * (16.0 / 32) ** 2 * 0.1)
    d["rho"] = 1 + 1e-6 * s
    lag3d_io.write_l3d(fb, d, time=1.0, gamma=m["gamma"], box=m["box"])
    res = analyse(snaps, fb, out=os.path.join(tmp, "sw3d_selftest"))
    ok = abs(res["order_pairs"][-1] - 2) < 0.1 and res["PASS"]
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snaps", nargs="*", default=[])
    ap.add_argument("--boosted", default=None)
    ap.add_argument("--ic", default=None, help="IC JSON sidecar with the wave params")
    ap.add_argument("-o", "--out", default="sw3d")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tmp", default="/tmp/lag3d_selftest/sw")
    a = ap.parse_args()
    if a.selftest:
        sys.path.insert(0, HERE)
        sys.exit(selftest(a.tmp))
    if not a.snaps:
        ap.error("--snaps required")
    r = analyse(a.snaps, a.boosted, a.ic, a.out)
    sys.exit(0 if r["PASS"] else 2)
