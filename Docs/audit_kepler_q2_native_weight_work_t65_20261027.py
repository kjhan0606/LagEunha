"""Pressure work of the proposed Wdot on an accepted W=0 Kepler state."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from export_kepler_q2_weight_rate_csr_t65_20261027 import finite_power_edges
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = 0.00193983


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("reference_rates", type=Path)
    parser.add_argument("case", choices=("aligned", "jitter"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if receipt["time"] != 6.5 or not all(receipt["checks"].values()):
        raise ValueError("Expected an accepted t=6.5 native checkpoint")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Native endpoint changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3, p0["y"].astype(float)-3))
    xy = np.column_stack((p["x"].astype(float)-3, p["y"].astype(float)-3))
    r0 = np.linalg.norm(xy0, axis=1)
    active = (r0 >= 0.5) & (r0 < 2.0)
    with np.load(args.reference_rates) as saved:
        rate = saved[args.case].copy()
        if not np.array_equal(saved["active"].astype(bool), active):
            raise ValueError("Rate active IDs differ from initial disk")
    if not np.all(np.isfinite(rate)) or np.any(rate[~active] != 0):
        raise ValueError("Invalid or nonzero inactive rate")
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(active.sum())
    area, meta, triangles, centers = power_areas(
        xy, np.zeros(N), np.flatnonzero(active), lookup, return_geometry=True)
    if meta["hidden_active"] or np.any(area <= 0):
        raise ValueError("Invalid W0 power diagram")
    i, j, c1, c2 = finite_power_edges(triangles, centers, active)
    length = np.linalg.norm(c2-c1, axis=1)
    distance = np.linalg.norm(xy[j]-xy[i], axis=1)
    g = length/(2.0*distance)
    face_volume_rate = g*(rate[i]-rate[j])
    weighted_volume_rate = np.zeros(N)
    np.add.at(weighted_volume_rate, i, face_volume_rate)
    np.add.at(weighted_volume_rate, j, -face_volume_rate)
    pressure = p["pressure"].astype(float)
    face_pressure_work = (pressure[i]-pressure[j])*face_volume_rate
    sum_by_face = float(np.sum(face_pressure_work, dtype=np.longdouble))
    sum_by_cell = float(np.sum(pressure*weighted_volume_rate, dtype=np.longdouble))
    initial_thermal_energy = float(np.sum(p["ie"].astype(float), dtype=np.longdouble))
    result = {
        "case": args.case,
        "accepted_endpoint_sha256": receipt["endpoint_sha256"],
        "reference_rates_sha256": digest(args.reference_rates),
        "active_cells": int(active.sum()),
        "undirected_faces_touching_active": int(len(i)),
        "rate_abs_max": float(np.max(np.abs(rate))),
        "total_weight_volume_rate": float(np.sum(weighted_volume_rate, dtype=np.longdouble)),
        "disk_weight_volume_rate_sum": float(np.sum(weighted_volume_rate[active], dtype=np.longdouble)),
        "pressure_weight_power_by_faces": sum_by_face,
        "pressure_weight_power_by_cells": sum_by_cell,
        "pressure_weight_power_abs_face_sum": float(
            np.sum(np.abs(face_pressure_work), dtype=np.longdouble)),
        "pressure_weight_power_face_cell_difference": abs(sum_by_face-sum_by_cell),
        "one_stage_dt": DT,
        "one_stage_unbalanced_energy_over_current_thermal":
            float(abs(DT*sum_by_face)/initial_thermal_energy),
        "one_stage_unbalanced_energy_over_abs_total_energy":
            float(abs(DT*sum_by_face)/5.352708),
        "scope": "First-order pressure work at an accepted native W0 t=6.5 state using independent CPU power geometry and the validated Jacobi16 rate. This is the energy term that a fixed-mass weighted evolution must account for; no weighted hydro has run.",
    }
    if (len(i) != 148557 or abs(result["total_weight_volume_rate"]) > 1e-10
            or result["pressure_weight_power_face_cell_difference"] > 1e-12):
        raise ValueError("Face/cell pressure-weight work audit failed")
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
