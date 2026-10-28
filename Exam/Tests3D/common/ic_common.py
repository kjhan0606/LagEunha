"""Shared helpers for the IC generators: argument defaults, totals report,
and a JSON sidecar with the parameters so analysis scripts do not have to
guess them."""
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lag3d_io  # noqa: E402


def report(d, meta, extra=None, pot_energy=None):
    t = lag3d_io.totals(d)
    n = len(d["x"])
    Etot = t["Ekin"] + t["Eint"] + (pot_energy or 0.0)
    print("np=%d  M=%.12g  P=(%.3e, %.3e, %.3e)" % (n, t["M"], t["Px"], t["Py"], t["Pz"]))
    print("Ekin=%.12g  Eint=%.12g%s  Etot=%.12g" % (
        t["Ekin"], t["Eint"],
        "" if pot_energy is None else "  Epot=%.12g" % pot_energy, Etot))
    bad = []
    for k in ("x", "y", "z", "vx", "vy", "vz", "mass", "u"):
        if not np.all(np.isfinite(d[k])):
            bad.append("non-finite " + k)
    if np.any(d["mass"] <= 0):
        bad.append("mass<=0")
    if np.any(d["u"] < 0):
        bad.append("u<0")
    b = meta["box"]
    for ax, lo, hi in (("x", b[0], b[1]), ("y", b[2], b[3]), ("z", b[4], b[5])):
        if np.any(d[ax] < lo) or np.any(d[ax] >= hi):
            bad.append("%s outside box" % ax)
    ids = d.get("id", np.arange(n))
    if len(np.unique(ids)) != n:
        bad.append("duplicate ids")
    out = dict(t, np=n, Etot=Etot)
    if pot_energy is not None:
        out["Epot"] = pot_energy
    if extra:
        out.update(extra)
    if bad:
        print("IC CHECK FAILED: " + ", ".join(bad))
    else:
        print("IC checks: finite, m>0, u>=0, inside box, unique ids -> OK")
    return out, bad


def write(path, d, meta, params, totals):
    lag3d_io.write_l3d(path, d, time=0.0, gamma=meta["gamma"], box=meta["box"],
                       periodic=meta["periodic"], G=meta.get("G", 0.0),
                       softening=meta.get("softening", 0.0))
    side = os.path.splitext(path)[0] + ".json"
    with open(side, "w") as f:
        json.dump(dict(params=params, meta=meta, totals=totals), f, indent=1,
                  default=float)
    print("wrote %s and %s" % (path, side))


def load_sidecar(path):
    side = os.path.splitext(path)[0] + ".json"
    if os.path.exists(side):
        with open(side) as f:
            return json.load(f)
    return None
