"""Read-only radial diagnostics for accepted fixed-polygon Kepler snapshots.

Bins follow the initial generators, so their angular momentum changes track
the same fixed-mass cells.  They are not Eulerian fluxes through a circle.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
GAMMA_MINUS_ONE = float(np.float32(1.666667)) - 1.0
DOMAIN_AREA = 512 * 2.75**2 * np.sin(2 * np.pi / 1024)
BINS = (
    (0.0, 0.2), (0.2, 0.5), (0.5, 0.75), (0.75, 1.0),
    (1.0, 1.5), (1.5, 2.0), (0.5, 2.0),
    (2.0, 2.35), (2.35, 2.5), (2.5, 2.75),
)


def weighted_abs(values, mass):
    return float(np.sum(np.abs(values) * mass, dtype=np.longdouble)
                 / np.sum(mass, dtype=np.longdouble))


def global_invariants(particles):
    x = particles["x"].astype(float) - 3.0
    y = particles["y"].astype(float) - 3.0
    vx = particles["vx"].astype(float)
    vy = particles["vy"].astype(float)
    mass = particles["mass"].astype(float)
    angular_momentum = np.sum(mass * (x * vy - y * vx), dtype=np.longdouble)
    energy = np.sum(.5 * mass * (vx * vx + vy * vy)
                    + particles["ie"].astype(float)
                    - mass / np.sqrt(x * x + y * y + .01175**2),
                    dtype=np.longdouble)
    return float(angular_momentum), float(energy)


def radial_rows(initial, current, time):
    x0 = initial["x"].astype(float) - 3.0
    y0 = initial["y"].astype(float) - 3.0
    x = current["x"].astype(float) - 3.0
    y = current["y"].astype(float) - 3.0
    r0 = np.hypot(x0, y0)
    r = np.hypot(x, y)
    mass = initial["mass"].astype(float)
    vx = current["vx"].astype(float)
    vy = current["vy"].astype(float)
    vcirc = r0 * (r0 * r0 + .01175**2)**(-.75)
    vr = (x * vx + y * vy) / np.maximum(r, 1e-30)
    phase = time * (r0 * r0 + .01175**2)**(-.75)
    orbit_error = np.hypot(x - (x0 * np.cos(phase) - y0 * np.sin(phase)),
                           y - (x0 * np.sin(phase) + y0 * np.cos(phase)))
    density_change = current["den"].astype(float) / initial["den"].astype(float) - 1.0
    pressure_change = (current["pressure"].astype(float)
                       / initial["pressure"].astype(float) - 1.0)
    entropy_log = (np.log(current["ie"].astype(float) / initial["ie"].astype(float))
                   + GAMMA_MINUS_ONE * np.log(current["volume"].astype(float)
                                              / initial["volume"].astype(float)))
    l0 = mass * (x0 * initial["vy"].astype(float)
                 - y0 * initial["vx"].astype(float))
    l1 = mass * (x * vy - y * vx)
    rows = []
    for lo, hi in BINS:
        selected = (r0 >= lo) & (r0 < hi)
        if not np.any(selected):
            raise ValueError(f"Empty initial radial bin {lo}:{hi}")
        w = mass[selected]
        speed_selected = selected & (r0 >= .05)
        wv = mass[speed_selected]
        initial_l = np.sum(l0[selected], dtype=np.longdouble)
        current_l = np.sum(l1[selected], dtype=np.longdouble)
        rows.append({
            "rlo": lo, "rhi": hi, "sites": int(np.sum(selected)),
            "mass": float(np.sum(w, dtype=np.longdouble)),
            "density_change_abs_mass_weighted": weighted_abs(density_change[selected], w),
            "pressure_change_abs_mass_weighted": weighted_abs(pressure_change[selected], w),
            "entropy_log_abs_mass_weighted": weighted_abs(entropy_log[selected], w),
            "entropy_log_signed_mass_weighted": float(
                np.sum(entropy_log[selected] * w, dtype=np.longdouble)
                / np.sum(w, dtype=np.longdouble)),
            "orbit_position_error_abs_mass_weighted": weighted_abs(orbit_error[selected], w),
            "radial_speed_rms_over_circular": float(np.sqrt(
                np.sum((vr[speed_selected] / vcirc[speed_selected])**2 * wv,
                       dtype=np.longdouble) / np.sum(wv, dtype=np.longdouble))),
            "angular_momentum_relative": float(current_l / initial_l - 1),
        })
    crossings = {}
    for boundary in (.5, 2.0, 2.35):
        across = (r0 < boundary) != (r < boundary)
        crossings[str(boundary)] = {
            "sites": int(np.sum(across)),
            "mass": float(np.sum(mass[across], dtype=np.longdouble)),
        }
    return rows, crossings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoints", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if abs(snapshot_time(args.initial, N)) > 1e-12:
        raise ValueError("Initial snapshot is not exact t=0")
    initial, _ = ordered(args.initial, N)
    initial_l, initial_e = global_invariants(initial)
    results = []
    last_time = 0.0
    for path in args.endpoints:
        time = snapshot_time(path, N)
        if time <= last_time:
            raise ValueError("Endpoints must have strictly increasing accepted times")
        current, _ = ordered(path, N)
        rows, crossings = radial_rows(initial, current, time)
        current_l, current_e = global_invariants(current)
        volume_sum = float(np.sum(current["volume"], dtype=np.longdouble))
        results.append({
            "endpoint": str(path.resolve()), "sha256": digest(path), "time": time,
            "fixed_mass_exact": bool(np.array_equal(initial["mass"], current["mass"])),
            "W0_exact": bool(np.all(initial["w2"] == 0) and np.all(current["w2"] == 0)),
            "positive_state": bool(all(np.all(current[key] > 0)
                                       for key in ("mass", "volume", "den", "ie", "pressure"))),
            "volume_sum": volume_sum,
            "area_tiles_polygon": bool(abs(volume_sum - DOMAIN_AREA) < 1e-5),
            "whole_domain_angular_momentum": current_l,
            "whole_domain_angular_momentum_change": current_l - initial_l,
            "whole_domain_angular_momentum_relative": current_l / initial_l - 1,
            "whole_domain_energy": current_e,
            "whole_domain_energy_relative_change": (current_e - initial_e) / abs(initial_e),
            "radial_bins": rows, "generator_crossings": crossings,
        })
        last_time = time
    result = {
        "initial": str(args.initial.resolve()), "initial_sha256": digest(args.initial),
        "site_count": N, "initial_whole_domain_angular_momentum": initial_l,
        "initial_whole_domain_energy": initial_e, "series": results,
        "scope": ("Accepted native snapshots; bins track initial fixed-mass generators. "
                  "Whole-domain angular momentum change includes wall torque and "
                  "numerical error; it is not an independent wall-face ledger. "
                  "No face-integrated flux or causal attribution of disk changes "
                  "is inferred."),
    }
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(rendered)
    print(rendered)


if __name__ == "__main__":
    main()
