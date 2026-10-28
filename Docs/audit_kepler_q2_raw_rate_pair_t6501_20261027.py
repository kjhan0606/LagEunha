"""Compare native projected and raw-rate replay from the same q2 checkpoint."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_bounded_native_evolution_20261027 import energy
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


BINS = ((.5, .52), (.52, .55), (.55, .65), (.5, 2.))
STAGE = re.compile(r"\[GPU_VOLW2\] stage=(\d+) active=(\d+) faces=(\d+) "
                   r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=([\deE+.-]+) s")


def weighted_abs(values, masses):
    return float(np.average(np.abs(values), weights=masses))


def stages(run):
    rows = []
    for line in (run / "log").read_text(errors="replace").splitlines():
        match = STAGE.fullmatch(line)
        if match:
            rows.append({"stage": int(match.group(1)), "active": int(match.group(2)),
                         "faces": int(match.group(3)), "alpha": float(match.group(4)),
                         "work_rate": float(match.group(5)),
                         "seconds": float(match.group(6))})
    if len(rows) != 4 or [row["stage"] for row in rows] != list(range(4)):
        raise ValueError("Expected four native RK4 weight stages")
    return rows


def metrics(initial, start, end, r0):
    x0 = initial["x"].astype(float)-3.
    y0 = initial["y"].astype(float)-3.
    ell0 = x0*initial["vy"].astype(float)-y0*initial["vx"].astype(float)
    xp = start["x"].astype(float)-3.
    yp = start["y"].astype(float)-3.
    ep = xp*start["vy"].astype(float)-yp*start["vx"].astype(float)
    xe = end["x"].astype(float)-3.
    ye = end["y"].astype(float)-3.
    ee = xe*end["vy"].astype(float)-ye*end["vx"].astype(float)
    mass = initial["mass"].astype(float)
    bins = {}
    for lo, hi in BINS:
        select = (r0 >= lo) & (r0 < hi)
        m = mass[select]
        if not np.all(ell0[select] > 0):
            raise ValueError("Expected positive initial disk orbital angular momentum")
        bins[f"{lo:g}:{hi:g}"] = {
            "sites": int(select.sum()),
            "volume_step_rms": float(np.sqrt(np.average(
                (end["volume"][select].astype(float)
                 /start["volume"][select].astype(float)-1.)**2, weights=m))),
            "volume_step_mean_abs": weighted_abs(
                end["volume"][select].astype(float)
                /start["volume"][select].astype(float)-1., m),
            "density_error_mean_abs": weighted_abs(
                end["den"][select].astype(float)
                /initial["den"][select].astype(float)-1., m),
            "pressure_error_mean_abs": weighted_abs(
                end["pressure"][select].astype(float)
                /initial["pressure"][select].astype(float)-1., m),
            "ell_step_mean_abs_relative_initial": weighted_abs(
                (ee[select]-ep[select])/ell0[select], m),
            "ell_error_mean_abs_relative_initial": weighted_abs(
                ee[select]/ell0[select]-1., m),
        }
    return bins


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline_receipt", type=Path)
    parser.add_argument("raw_receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    base = json.loads(args.baseline_receipt.read_text())
    raw = json.loads(args.raw_receipt.read_text())
    if (base["time"] != 6.501 or raw["time"] != 6.501
            or base["parent_time"] != 6.5 or raw["parent_time"] != 6.5
            or base["initial_sha256"] != raw["initial_sha256"]
            or base["parent_sha256"] != raw["parent_sha256"]
            or base["binary_sha256"] == raw["binary_sha256"]
            or not all(base["checks"].values())):
        raise ValueError("Invalid projected/raw native comparison lineage")
    for receipt in (base, raw):
        for key, sha_key in (("initial", "initial_sha256"),
                            ("parent", "parent_sha256"),
                            ("endpoint", "endpoint_sha256")):
            if digest(Path(receipt[key])) != receipt[sha_key]:
                raise ValueError(f"Changed {key} snapshot")
    n = base["n"]
    initial, _ = ordered(Path(base["initial"]), n)
    start, _ = ordered(Path(base["parent"]), n)
    projected, _ = ordered(Path(base["endpoint"]), n)
    unprojected, _ = ordered(Path(raw["endpoint"]), n)
    r0 = np.hypot(initial["x"].astype(float)-3.,
                  initial["y"].astype(float)-3.)
    base_stages = stages(Path(base["endpoint"]).parent)
    raw_stages = stages(Path(raw["endpoint"]).parent)
    if any(row["alpha"] != 0. for row in raw_stages):
        raise ValueError("Raw-rate native stage applied nonzero projection")
    dt = raw["time"]-raw["parent_time"]
    weights = np.array([1., 2., 2., 1.])/6.
    work_integral = float(dt*np.dot(weights,
                                   [row["work_rate"] for row in raw_stages]))
    parent_energy = energy(start)
    projected_energy = energy(projected)
    raw_energy = energy(unprojected)
    results = {
        "baseline_receipt": str(args.baseline_receipt),
        "raw_receipt": str(args.raw_receipt),
        "initial_sha256": base["initial_sha256"],
        "shared_parent_sha256": base["parent_sha256"],
        "baseline_endpoint_sha256": base["endpoint_sha256"],
        "raw_endpoint_sha256": raw["endpoint_sha256"],
        "time_start": 6.5, "time_end": 6.501,
        "baseline_checks": base["checks"], "raw_checks": raw["checks"],
        "baseline_total_energy_relative_t0": base["total_energy_relative"],
        "raw_total_energy_relative_t0": raw["total_energy_relative"],
        "parent_energy": parent_energy,
        "baseline_energy_step": projected_energy-parent_energy,
        "raw_energy_step": raw_energy-parent_energy,
        "raw_minus_baseline_energy_step": raw_energy-projected_energy,
        "raw_minus_baseline_energy_step_over_abs_parent":
            (raw_energy-projected_energy)/abs(parent_energy),
        "native_weight_stages_baseline": base_stages,
        "native_weight_stages_raw": raw_stages,
        "raw_RK4_pressure_weight_work_integral": work_integral,
        "baseline": metrics(initial, start, projected, r0),
        "raw": metrics(initial, start, unprojected, r0),
        "scope": "Same accepted t=6.5 parent. Work integral is signed RK4 "
                 "pressure-weight rate diagnostic, not the full discrete "
                 "total-energy error. The existing projected_weight_work "
                 "gate uses a broad 1e-10 absolute threshold and must not "
                 "be interpreted as zero work in raw mode.",
    }
    with args.output.open("x") as out:
        json.dump(results, out, indent=2, allow_nan=False)
        out.write("\n")
    for key in ("0.5:0.52", "0.52:0.55", "0.5:2"):
        a, b = results["baseline"][key], results["raw"][key]
        print(key, "volume_rms projected/raw", a["volume_step_rms"],
              b["volume_step_rms"], "density projected/raw",
              a["density_error_mean_abs"], b["density_error_mean_abs"],
              "ell_increment projected/raw", a["ell_step_mean_abs_relative_initial"],
              b["ell_step_mean_abs_relative_initial"])
    print("raw_work_integral", work_integral,
          "total_energy_relative_t0", base["total_energy_relative"],
          raw["total_energy_relative"])


if __name__ == "__main__":
    main()
