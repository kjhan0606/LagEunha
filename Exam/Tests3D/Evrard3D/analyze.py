#!/usr/bin/env python3
"""Analyse an Evrard collapse run: energy curves and the t = 0.8 profile.

  analyze.py --log run/log [--snap run/snap_t0.800.l3d] [-o evrard3d]
  analyze.py --table energies.txt ...     (columns: t Ekin Eth Epot Etot)
  analyze.py --selftest

References (nothing here is digitised from a paper):
  * ref/evrard1d_N2000_energy.txt, ref/evrard1d_N2000_prof_t0.800.txt:
    the 1D spherical solution computed by evrard1d_ref.c in this directory
    (Lagrangian, compatible energy, exact self-gravity; dE/|E0| = 1.5e-5 at
    t = 3). Its t = 0.8 profile matches HydroCode1D to 0.3% in log rho
    (selftest). Key numbers: Ekin max 0.450 at t = 0.88, Eth max 1.758 at
    t = 1.07, Epot min -2.543 at t = 1.05, Etot = -0.61667.
  * ref/evrardCollapse3D_exact.txt (getReference.sh): HydroCode1D t = 0.8
    profile, the SWIFT reference file.
A 3D run at finite resolution and softening will not reach the 1D central
density; the energy-curve criteria below are therefore loose.

PASS:
  E1  max_{t<=t_end} |Etot(t) - Etot(0)| / |Etot(0)| <= 0.01      (primary)
  E2  |t(Ekin max) - 0.88| <= 0.10  and  |t(Eth max) - 1.07| <= 0.15
  E3  mean_t |Eth(t) - Eth_ref(t)| / max(Eth_ref) <= 0.15 over the run
  E4  (only with --snap near t=0.8) mean |log10(rho/rho_ref)| over
      0.05 < r < 0.8 <= 0.10
E3 and E4 are this suite's choices, not community thresholds.
"""
import argparse
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import lag3d_io  # noqa: E402
import radial  # noqa: E402

REF_E = os.path.join(HERE, "ref", "evrard1d_N2000_energy.txt")
REF_P = os.path.join(HERE, "ref", "evrard1d_N2000_prof_t0.800.txt")
REF_H = os.path.join(HERE, "ref", "evrardCollapse3D_exact.txt")


def load_energies(log=None, table=None):
    if table:
        a = np.loadtxt(table)
        return dict(t=a[:, 0], Ekin=a[:, 1], Eint=a[:, 2], Epot=a[:, 3], Etot=a[:, 4])
    e = lag3d_io.parse_energy_log(log)
    need = ("t", "Ekin", "Eint", "Epot")
    if not all(k in e for k in need):
        raise SystemExit("log has no [E3D] lines with %s" % (need,))
    if "Etot" not in e:
        e["Etot"] = e["Ekin"] + e["Eint"] + e["Epot"]
    return e


def energy_metrics(e, ref):
    t = e["t"]
    E0 = e["Etot"][0]
    dE = np.abs(e["Etot"] - E0) / abs(E0)
    tk = float(t[np.argmax(e["Ekin"])])
    tt = float(t[np.argmax(e["Eint"])])
    rt, rk, rth, rp = ref[:, 0], ref[:, 1], ref[:, 2], ref[:, 3]
    sel = t <= rt.max()
    Eth_ref = np.interp(t[sel], rt, rth)
    L1th = float(np.mean(np.abs(e["Eint"][sel] - Eth_ref)) / rth.max())
    tk_ref = float(rt[np.argmax(rk)])
    tt_ref = float(rt[np.argmax(rth)])
    m = dict(E0=float(E0), max_dE=float(dE.max()), t_Ekin_max=tk, t_Eth_max=tt,
             t_Ekin_max_ref=tk_ref, t_Eth_max_ref=tt_ref,
             Ekin_max=float(e["Ekin"].max()), Eth_max=float(e["Eint"].max()),
             Epot_min=float(e["Epot"].min()), Ekin_max_ref=float(rk.max()),
             Eth_max_ref=float(rth.max()), Epot_min_ref=float(rp.min()),
             L1_Eth=L1th, t_end=float(t.max()))
    c = dict(E1=m["max_dE"] <= 0.01,
             E2=abs(tk - 0.88) <= 0.10 and abs(tt - 1.07) <= 0.15,
             E3=L1th <= 0.15)
    return m, c


def profile_metrics(snap):
    d, meta = lag3d_io.read_l3d(snap)
    g = meta["gamma"]
    tot = lag3d_io.totals(d)
    com = np.array([np.sum(d["mass"] * d[k]) for k in ("x", "y", "z")]) / tot["M"]
    r, vr, _ = radial.radii(d, com)
    rho = lag3d_io.density(d)
    P = (g - 1) * rho * d["u"]
    w = d["vol"] if np.all(d["vol"] > 0) else d["mass"] / rho
    edges = np.geomspace(0.01, 1.5, 50)
    rc, rb, cnt = radial.binned(r, rho, w, edges)
    _, vb, _ = radial.binned(r, vr, w, edges)
    _, pb, _ = radial.binned(r, P, w, edges)
    ref = np.loadtxt(REF_P)
    ok = (rc > 0.05) & (rc < 0.8) & np.isfinite(rb) & (cnt > 5)
    err = float(np.mean(np.abs(np.log10(rb[ok] / np.interp(rc[ok], ref[:, 0], ref[:, 1])))))
    return dict(t=meta["time"], L1_logrho=err), (rc, rb, vb, pb)


def compare_1d_to_hydrocode1d():
    """Selftest helper: our 1D reference against the published HydroCode1D."""
    if not os.path.exists(REF_H):
        return None
    h = np.loadtxt(REF_H)
    mine = np.loadtxt(REF_P)
    sel = (h[:, 0] > 0.02) & (h[:, 0] < 0.9)
    out = {}
    for k, name in ((1, "rho"), (3, "P")):
        out["L1_log" + name] = float(np.mean(np.abs(np.log10(
            np.interp(h[sel, 0], mine[:, 0], mine[:, k]) / h[sel, k]))))
    v = np.interp(h[sel, 0], mine[:, 0], mine[:, 2])
    out["L1_v_rel"] = float(np.mean(np.abs(v - h[sel, 2])) / np.mean(np.abs(h[sel, 2])))
    return out


def analyse(log=None, table=None, snap=None, out="evrard3d", plot=True):
    ref = np.loadtxt(REF_E)
    e = load_energies(log, table)
    m, c = energy_metrics(e, ref)
    res = dict(energy=m)
    prof = None
    if snap:
        pm, prof = profile_metrics(snap)
        res["profile"] = pm
        c["E4"] = bool(pm["L1_logrho"] <= 0.10 and abs(pm["t"] - 0.8) < 0.05)
    res["criteria"] = {k: bool(v) for k, v in c.items()}
    res["PASS"] = bool(all(c.values()))
    with open(out + "_summary.json", "w") as f:
        json.dump(res, f, indent=1)
    if plot:
        plt = lag3d_io.get_plt()
        fig, ax = plt.subplots(1, 2 if prof is None else 3, figsize=(13 if prof else 10, 4))
        for k, lab, col in (("Ekin", "E_kin", "C0"), ("Eint", "E_th", "C1"),
                            ("Epot", "E_pot", "C2"), ("Etot", "E_tot", "k")):
            ax[0].plot(e["t"], e[k], "-", color=col, label=lab + " (run)")
        for j, col in ((1, "C0"), (2, "C1"), (3, "C2"), (4, "k")):
            ax[0].plot(ref[:, 0], ref[:, j], "--", color=col, lw=1)
        ax[0].set_xlabel("t"); ax[0].set_ylabel("energy"); ax[0].legend(fontsize=7)
        ax[0].set_title("dashed: 1D reference (evrard1d_ref N=2000)")
        ax[1].semilogy(e["t"], np.abs(e["Etot"] - e["Etot"][0]) / abs(e["Etot"][0]) + 1e-16)
        ax[1].axhline(0.01, color="r", ls=":")
        ax[1].set_xlabel("t"); ax[1].set_ylabel("|dEtot/Etot0|")
        if prof is not None:
            rc, rb, vb, pb = prof
            r1 = np.loadtxt(REF_P)
            ax[2].loglog(rc, rb, "o", ms=3, label="run")
            ax[2].loglog(r1[:, 0], r1[:, 1], "k-", lw=1, label="1D ref")
            if os.path.exists(REF_H):
                h = np.loadtxt(REF_H)
                ax[2].loglog(h[:, 0], h[:, 1], "r:", lw=1, label="HydroCode1D")
            ax[2].set_xlim(0.01, 1.5); ax[2].set_ylim(1e-2, 1e4)
            ax[2].set_xlabel("r"); ax[2].set_ylabel("rho"); ax[2].legend(fontsize=7)
        lag3d_io.savefig(fig, out + "_energy.png")
    print(json.dumps(res, indent=1))
    print("PASS" if res["PASS"] else "FAIL")
    return res


def selftest(tmp):
    os.makedirs(tmp, exist_ok=True)
    ok = True
    ref = np.loadtxt(REF_E)
    drift = abs(ref[-1, 4] - ref[0, 4]) / abs(ref[0, 4])
    print("1D reference: Etot drift to t=3 = %.2e" % drift)
    ok &= drift < 1e-4
    cmp = compare_1d_to_hydrocode1d()
    if cmp is None:
        print("HydroCode1D file absent (run getReference.sh); comparison skipped")
    else:
        print("1D reference vs HydroCode1D at t=0.8:", cmp)
        ok &= cmp["L1_logrho"] < 0.02 and cmp["L1_logP"] < 0.03 and cmp["L1_v_rel"] < 0.03
    # synthetic "good" run: the reference with noise, written as [E3D] log lines
    rng = np.random.default_rng(0)
    logp = os.path.join(tmp, "good.log")
    with open(logp, "w") as f:
        for i, row in enumerate(ref[::3]):
            n = 1 + 0.01 * rng.standard_normal(3)
            ek, eth, ep = row[1] * n[0], row[2] * n[1], row[3]
            et = row[4] * (1 + 2e-3 * row[0] / 3)
            f.write("Time= %g step= %d\n[E3D] step= %d t= %.6e dt= 1e-4 Ekin= %.8e Eint= %.8e Epot= %.8e Etot= %.8e\n"
                    % (row[0], i, i, row[0], ek, eth, ep, et))
    good = analyse(log=logp, out=os.path.join(tmp, "good"))
    ok &= good["PASS"]
    # profile code path: an IC lattice with the reference t=0.8 density
    sys.path.insert(0, HERE)
    import make_ic
    icf = os.path.join(tmp, "ev_ic_n24.l3d")
    make_ic.main(["--n", "24", "-o", icf])
    d, mt = lag3d_io.read_l3d(icf)
    r = np.sqrt(d["x"] ** 2 + d["y"] ** 2 + d["z"] ** 2)
    pr = np.loadtxt(REF_P)
    d["rho"] = np.interp(r, pr[:, 0], pr[:, 1])
    sn = os.path.join(tmp, "ev_synth_t0.8.l3d")
    lag3d_io.write_l3d(sn, d, time=0.8, gamma=mt["gamma"], box=mt["box"], periodic=(0, 0, 0))
    pm, _ = profile_metrics(sn)
    print("synthetic t=0.8 snapshot: mean|log10 rho/ref| = %.4f" % pm["L1_logrho"])
    ok &= pm["L1_logrho"] < 0.02
    # synthetic "bad" run: 3% energy drift and a late collapse
    badp = os.path.join(tmp, "bad.txt")
    b = ref.copy()
    b[:, 0] *= 1.25
    b[:, 4] = b[:, 4] * (1 + 0.03 * b[:, 0] / 3)
    np.savetxt(badp, b)
    bad = analyse(table=badp, out=os.path.join(tmp, "bad"), plot=False)
    ok &= (not bad["PASS"]) and (not bad["criteria"]["E1"]) and (not bad["criteria"]["E2"])
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log")
    ap.add_argument("--table")
    ap.add_argument("--snap")
    ap.add_argument("-o", "--out", default="evrard3d")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--tmp", default="/tmp/lag3d_selftest/evrard")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest(a.tmp))
    if not (a.log or a.table):
        ap.error("--log or --table required")
    r = analyse(a.log, a.table, a.snap, a.out)
    sys.exit(0 if r["PASS"] else 2)
