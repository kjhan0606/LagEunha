"""Read-only radial and boundary audit for native Kepler snapshots of any N."""

import argparse
import json
from pathlib import Path
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


def ordered(path):
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, stride, stress) != (-1, 576, 1) or path.stat().st_size != 16+n*stride+52:
        raise ValueError(f"Incomplete or incompatible snapshot {path}")
    particles, meta = read_snapshot(path)
    ids = (particles["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if not np.array_equal(np.sort(ids), np.arange(n)):
        raise ValueError(f"Invalid cell IDs {path}")
    return particles[np.argsort(ids)], meta


EPS = .01175
GAMMA_MINUS_ONE = float(np.float32(1.666667))-1.


def mean(values, weights):
    return float(np.sum(values*weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def taper_speed(r, r0=2.35, r1=2.75):
    z = -1. / np.sqrt(r*r+EPS*EPS)
    z0, z1 = -1. / np.sqrt(r0*r0+EPS*EPS), -1. / np.sqrt(r1*r1+EPS*EPS)
    u = np.clip((z-z0)/(z1-z0), 0., 1.)
    smooth = 1.-u**3*(10.+u*(-15.+6.*u))
    return np.where(r<=r0, 1., np.where(r>=r1, 0., smooth))


def potential(r, mode):
    z = -1. / np.sqrt(r*r+EPS*EPS)
    if mode == "baseline":
        return z
    r0, r1 = 2.35, 2.75
    z0, z1 = -1. / np.sqrt(r0*r0+EPS*EPS), -1. / np.sqrt(r1*r1+EPS*EPS)
    u = np.clip((z-z0)/(z1-z0), 0., 1.)
    u2, u4 = u*u, u*u*u*u
    h = u + u4*(-5.+u*(6.+u*(-2.+u*(100./7.+u*(-75./2.+
        u*(115./3.+u*(-18.+u*36./11.)))))))
    return np.where(r<=r0, z,
                    np.where(r>=r1, z0+(z1-z0)*181./462., z0+(z1-z0)*h))


def energy(p, mode):
    x = p["x"].astype(float)-3.
    y = p["y"].astype(float)-3.
    r = np.hypot(x,y)
    m = p["mass"].astype(float)
    v2 = p["vx"].astype(float)**2+p["vy"].astype(float)**2
    return float(np.sum(.5*m*v2+p["ie"].astype(float)+m*potential(r,mode),
                        dtype=np.longdouble))


def region(mask, m, density_change, pressure_change, entropy_log,
           vr, l0, lt, orbit_error, radius_current):
    if not np.any(mask):
        return None
    w = m[mask]
    initial_l = np.sum(l0[mask]*w, dtype=np.longdouble)
    final_l = np.sum(lt[mask]*w, dtype=np.longdouble)
    return {
        "cells": int(np.sum(mask)),
        "initial_mass": float(np.sum(w, dtype=np.longdouble)),
        "density_change_abs_mass_weighted": mean(abs(density_change[mask]),w),
        "pressure_change_abs_mass_weighted": mean(abs(pressure_change[mask]),w),
        "entropy_log_signed_mass_weighted": mean(entropy_log[mask],w),
        "entropy_log_abs_mass_weighted": mean(abs(entropy_log[mask]),w),
        "radial_velocity_abs_mass_weighted": mean(abs(vr[mask]),w),
        "radial_velocity_abs_max": float(np.max(abs(vr[mask]))),
        "angular_momentum_relative": float(final_l/initial_l-1.) if abs(initial_l)>0 else None,
        "orbit_position_error_abs_mass_weighted": mean(orbit_error[mask],w),
        "current_radius_minmax": [float(np.min(radius_current[mask])),
                                  float(np.max(radius_current[mask]))],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--mode", choices=("baseline","taper"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    p0, m0 = ordered(args.initial)
    p, mt = ordered(args.endpoint)
    t = float(mt["t"])
    if len(p0)!=len(p) or float(m0["t"])!=0. or not t>0:
        raise ValueError("Expected matching exact t=0 and evolved snapshots")
    finite = all(np.isfinite(q[field]).all() for q in (p0,p) for field in
                 ("x","y","vx","vy","mass","volume","ie","den","pressure","w2"))
    mass_exact = bool(np.array_equal(p0["mass"],p["mass"]))
    w0 = bool(np.all(p0["w2"]==0) and np.all(p["w2"]==0))
    positive = bool(np.all(p0["ie"]>0) and np.all(p["ie"]>0) and
                    np.all(p0["volume"]>0) and np.all(p["volume"]>0))
    if not (finite and mass_exact and w0 and positive):
        raise ValueError("Native state failed finite/fixed-mass/W0/positive gate")
    x0, y0 = p0["x"].astype(float)-3., p0["y"].astype(float)-3.
    x, y = p["x"].astype(float)-3., p["y"].astype(float)-3.
    r0, rt = np.hypot(x0,y0), np.hypot(x,y)
    m = p0["mass"].astype(float)
    omega = (r0*r0+EPS*EPS)**(-.75)
    if args.mode == "taper":
        omega *= taper_speed(r0)
    theta = t*omega
    xb = x0*np.cos(theta)-y0*np.sin(theta)
    yb = x0*np.sin(theta)+y0*np.cos(theta)
    dx = (x-xb+3.)%6.-3.
    dy = (y-yb+3.)%6.-3.
    orbit_error = np.hypot(dx,dy)
    vr = np.divide(x*p["vx"].astype(float)+y*p["vy"].astype(float),rt,
                   out=np.zeros_like(rt),where=rt>0)
    l0 = x0*p0["vy"].astype(float)-y0*p0["vx"].astype(float)
    lt = x*p["vy"].astype(float)-y*p["vx"].astype(float)
    density_change = p["den"].astype(float)/p0["den"].astype(float)-1.
    pressure_change = p["pressure"].astype(float)/p0["pressure"].astype(float)-1.
    entropy = np.log(p["ie"].astype(float)/p0["ie"].astype(float)) + \
              GAMMA_MINUS_ONE*np.log(p["volume"].astype(float)/p0["volume"].astype(float))
    edge = np.minimum.reduce((p["x"],p["y"],6.-p["x"],6.-p["y"]))
    bins = []
    for lo,hi in ((0.,.1),(.1,.2),(.2,.5),(.5,.75),(.75,1.),(1.,2.),
                  (2.,2.35),(2.35,2.75),(2.75,5.)):
        row = region((r0>=lo)&(r0<hi),m,density_change,pressure_change,
                     entropy,vr,l0,lt,orbit_error,rt)
        if row is not None:
            row.update(initial_radius_lo=lo,initial_radius_hi=hi)
            bins.append(row)
    e0, et = energy(p0,args.mode), energy(p,args.mode)
    # Cell-centre crossings measure Lagrangian transport across the inner
    # disk radius. They are not a reconstructed face flux through that circle.
    inner0, innert = r0 < .5, rt < .5
    into_inner = ~inner0 & innert
    out_of_inner = inner0 & ~innert
    inner_crossings = {
        "inward_count": int(np.sum(into_inner)),
        "inward_mass": float(np.sum(m[into_inner], dtype=np.longdouble)),
        "outward_count": int(np.sum(out_of_inner)),
        "outward_mass": float(np.sum(m[out_of_inner], dtype=np.longdouble)),
        "net_inward_mass": float(np.sum(m[into_inner], dtype=np.longdouble) -
                                 np.sum(m[out_of_inner], dtype=np.longdouble)),
        "initial_inner_mass_by_cell_center": float(np.sum(m[inner0], dtype=np.longdouble)),
        "current_inner_mass_by_cell_center": float(np.sum(m[innert], dtype=np.longdouble)),
        "definition": "Fixed masses of cells whose generators cross r=0.5; not a face-integrated Eulerian flux",
    }
    result = {
        "initial": str(args.initial), "endpoint": str(args.endpoint),
        "mode": args.mode, "time": t, "cells": len(p),
        "energy_initial": e0, "energy_final": et,
        "energy_change_relative": (et-e0)/abs(e0),
        "bins_by_initial_radius": bins,
        "current_edge_within_0p2": region(edge<.2,m,density_change,pressure_change,
                                           entropy,vr,l0,lt,orbit_error,rt),
        "initial_disk_to_inner_count": int(np.sum((r0>=.5)&(r0<2.)&(rt<.5))),
        "initial_disk_to_outer_count": int(np.sum((r0>=.5)&(r0<2.)&(rt>=2.))),
        "inner_disk_edge_cell_center_crossings": inner_crossings,
        "checks": {
            "finite": finite,"fixed_mass_exact": mass_exact,"W0_exact": w0,
            "positive_volume_internal_energy": positive,
            "rho_mass_over_volume_max_relative_below_1e_minus_6":
                bool(np.max(abs(p["den"].astype(float)/(m/p["volume"].astype(float))-1.))<1e-6),
        },
        "scope": "ID-matched native radial/edge diagnostics with analytic baseline or opt-in tapered external potential. Initial-radius bins do not establish causal propagation; total energy is a finite-volume ledger check, not a global benchmark pass.",
    }
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")
    print(json.dumps({"time":t,"mode":args.mode,"cells":len(p),
                      "energy_change_relative":result["energy_change_relative"],
                      "bins":[{"lo":b["initial_radius_lo"],"hi":b["initial_radius_hi"],
                               "cells":b["cells"],
                               "density":b["density_change_abs_mass_weighted"],
                               "entropy":b["entropy_log_abs_mass_weighted"],
                               "vr":b["radial_velocity_abs_mass_weighted"]} for b in bins]},indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("Native radial audit failed")


if __name__ == "__main__":
    main()
