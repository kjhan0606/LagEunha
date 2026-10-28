"""Show radial pressure, density and orbital error hidden by diskwide averages."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


BINS = ((.5, .52), (.52, .55), (.55, .65), (.65, .75), (.75, 1.), (1., 2.), (.5, 2.))


def weighted_abs(values, mass):
    return float(np.average(np.abs(values), weights=mass))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipts", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if len(args.receipts) < 2:
        raise ValueError("Need at least two native snapshots")
    receipts = [json.loads(path.read_text()) for path in args.receipts]
    initial_sha = receipts[0]["initial_sha256"]
    n = receipts[0]["n"]
    initial_path = Path(receipts[0]["initial"])
    if digest(initial_path) != initial_sha:
        raise ValueError("Changed initial snapshot")
    p0, _ = ordered(initial_path, n)
    x0 = p0["x"].astype(float)-3.
    y0 = p0["y"].astype(float)-3.
    r0 = np.hypot(x0, y0)
    l0 = x0*p0["vy"].astype(float)-y0*p0["vx"].astype(float)
    m0 = p0["mass"].astype(float)
    result = {"initial_sha256": initial_sha, "n": n, "rows": []}
    for path, receipt in zip(args.receipts, receipts):
        if (receipt["initial_sha256"] != initial_sha or receipt["n"] != n
                or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
            raise ValueError(f"Changed or mismatched endpoint: {path}")
        end, _ = ordered(Path(receipt["endpoint"]), n)
        x = end["x"].astype(float)-3.
        y = end["y"].astype(float)-3.
        l = x*end["vy"].astype(float)-y*end["vx"].astype(float)
        bins = {}
        for lo, hi in BINS:
            selected = (r0 >= lo) & (r0 < hi)
            mass = m0[selected]
            if not np.all(l0[selected] > 0):
                raise ValueError("Expected positive initial disk orbital motion")
            density_error = end["den"][selected].astype(float)/p0["den"][selected].astype(float)-1.
            pressure_error = (end["pressure"][selected].astype(float)
                              /p0["pressure"][selected].astype(float)-1.)
            log_density = np.log1p(density_error)
            log_pressure = np.log1p(pressure_error)
            gamma = float(np.float32(1.666667))
            adiabatic_pressure_residual = log_pressure-gamma*log_density
            entropy_log = (np.log(end["ie"][selected].astype(float)
                                  /p0["ie"][selected].astype(float))
                           +(gamma-1.)*np.log(end["volume"][selected].astype(float)
                                                /p0["volume"][selected].astype(float)))
            ell_error = l[selected]/l0[selected]-1.
            bins[f"{lo:g}:{hi:g}"] = {
                "sites": int(selected.sum()),
                "density_relative_mean_abs": weighted_abs(density_error, mass),
                "density_relative_p99_abs": float(np.quantile(np.abs(density_error), .99)),
                "pressure_relative_mean_abs": weighted_abs(pressure_error, mass),
                "pressure_relative_p99_abs": float(np.quantile(np.abs(pressure_error), .99)),
                "pressure_adiabatic_residual_log_mean_abs": weighted_abs(
                    adiabatic_pressure_residual, mass),
                "entropy_log_mean_abs": weighted_abs(entropy_log, mass),
                "pressure_over_gamma_density_rms_slope": float(
                    np.sum(mass*log_pressure*log_density)
                    /(gamma*np.sum(mass*log_density*log_density)))
                    if np.any(log_density) else None,
                "ell_relative_mean_abs": weighted_abs(ell_error, mass),
                "ell_relative_p99_abs": float(np.quantile(np.abs(ell_error), .99)),
            }
        result["rows"].append({"receipt": str(path), "time": receipt["time"],
                               "endpoint_sha256": receipt["endpoint_sha256"],
                               "gates_passed": sum(receipt["checks"].values()),
                               "gates_total": len(receipt["checks"]),
                               "total_energy_relative": receipt["total_energy_relative"],
                               "bins": bins})
    result["scope"] = ("Pointwise initial polar-ring labels. Diskwide 21-check "
                       "acceptance does not impose a separate first-ring "
                       "pressure or density threshold; this timeline "
                       "reports their growth explicitly.")
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for row in result["rows"]:
        ring = row["bins"]["0.5:0.52"]
        disk = row["bins"]["0.5:2"]
        print(row["time"], "first-ring density", ring["density_relative_mean_abs"],
              "pressure", ring["pressure_relative_mean_abs"],
              "ell", ring["ell_relative_mean_abs"],
              "diskwide density", disk["density_relative_mean_abs"])


if __name__ == "__main__":
    main()
