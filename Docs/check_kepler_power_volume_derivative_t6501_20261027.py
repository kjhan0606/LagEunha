"""Check the analytic power-face area rate against endpoint areas."""

import json
from pathlib import Path

import numpy as np

from audit_kepler_fresh_weighted_flux_proxy_t4952_20261027 import touching_finite_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RECEIPT = ROOT/"gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1/audit_t6501.json"


def scatter(i, j, value, n):
    out = np.zeros(n)
    np.add.at(out, i, value)
    np.add.at(out, j, -value)
    return out


def compare(predict, target):
    slope = np.dot(predict, target)/np.dot(predict, predict)
    cosine = np.dot(predict, target)/(
        np.linalg.norm(predict)*np.linalg.norm(target))
    return (float(np.sqrt(np.mean((predict-target)**2))),
            float(slope), float(cosine))


def main():
    audit = json.loads(RECEIPT.read_text())
    n = audit["n"]
    p0, _ = ordered(Path(audit["initial"]), n)
    a, _ = ordered(Path(audit["parent"]), n)
    b, _ = ordered(Path(audit["endpoint"]), n)
    active_mask = np.fromfile(RECEIPT.parent/"disk_active_mask.bin", dtype=np.uint8)
    active = np.flatnonzero(active_mask)
    lookup = np.full(n, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    xy0 = np.column_stack((a["x"].astype(float)-3., a["y"].astype(float)-3.))
    xy1 = np.column_stack((b["x"].astype(float)-3., b["y"].astype(float)-3.))
    w0 = a["w2"].astype(float)
    w1 = b["w2"].astype(float)
    dt = audit["time"]-audit["parent_time"]
    area0, _, tri, centers = power_areas(xy0, w0, active, lookup, True)
    area1, _ = power_areas(xy1, w1, active, lookup)
    i, j, c, length = touching_finite_faces(tri, centers, active_mask != 0)
    d = np.linalg.norm(xy0[j]-xy0[i], axis=1)
    v = (xy1-xy0)/dt
    wd = (w1-w0)/dt
    fixed = (np.einsum("ij,ij->i", xy0[j], v[j])
             - np.einsum("ij,ij->i", xy0[i], v[i])
             - np.einsum("ij,ij->i", v[j]-v[i], c))/d
    fixed_rate = scatter(i, j, length*fixed, n)
    weight_rate = scatter(i, j, length*(wd[i]-wd[j])/(2*d), n)
    truth = area1/area0-1.
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    first_ids = np.flatnonzero((r0 >= .5) & (r0 < .52))
    linear_first = dt*(fixed_rate[first_ids]+weight_rate[first_ids])/area0[lookup[first_ids]]
    for fraction in (.01, .1, .5, 1.):
        fraction_area, _ = power_areas(xy0+fraction*(xy1-xy0),
                                      w0+fraction*(w1-w0), active, lookup)
        finite_first = (fraction_area[lookup[first_ids]]/area0[lookup[first_ids]]-1.)/fraction
        print("first_ring_fraction", fraction, "linear_vs_finite", compare(linear_first, finite_first))
    print("area_native_start", np.max(np.abs(area0/a["volume"][active]-1.)))
    print("area_native_end", np.max(np.abs(area1/b["volume"][active]-1.)))
    for lo, hi in ((.5, .52), (.52, .55), (.55, .65), (.5, 2.)):
        take = (r0 >= lo) & (r0 < hi)
        area = np.zeros(n)
        area[active] = area0
        target = truth[lookup[np.flatnonzero(take)]]
        f = dt*fixed_rate[take]/area[take]
        w = dt*weight_rate[take]/area[take]
        print(lo, hi, "actual_native_vs_power", compare(
            b["volume"][take].astype(float)/a["volume"][take].astype(float)-1., target))
        for label, prediction in (("f+w", f+w), ("f-w", f-w),
                                  ("-f+w", -f+w), ("-f-w", -f-w)):
            print(" ", label, compare(prediction, target))


if __name__ == "__main__":
    main()
