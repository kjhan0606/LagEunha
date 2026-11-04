"""Assess frozen native volume work and linear pressure-force decomposition."""
from pathlib import Path
import argparse
import json

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("root", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
rows = []
for directory in sorted(args.root.iterdir()):
    if not directory.is_dir() or not (directory / "result.json").exists():
        continue
    result = json.loads((directory / "result.json").read_text())
    assert result["exit"] == 0 and result["accepted_steps"] == 0
    s = np.genfromtxt(directory / "unified_transport.csv", delimiter=",", names=True)
    heat = np.genfromtxt(directory / "unified_initial_heat.csv", delimiter=",", names=True)
    assert np.array_equal(s["id"], heat["id"])
    assert all(np.isfinite(s[name]).all() for name in s.dtype.names)
    volume = s["V"]
    radius = np.hypot(s["x"] - 0.5, s["y"] - 0.5)
    row = {"case": directory.name, "time": result["input_time"],
           "cells": len(s), "accepted_steps": 0,
           "max_linear_pressure_patch_error": float(s["pressure_moment_error"].max()),
           "max_dilation_moment_error": float(np.max(abs(s["Kxx"] + s["Kyy"] - 2))),
           "max_rotation_moment_error": float(np.max(abs(s["Kxy"] - s["Kyx"]))),
           "max_cpu_gpu_force_error": float(s["cpu_gpu_force_error"].max()),
           "max_native_thermal_rhs_error": float(np.max(abs(s["thermal_rhs_error"]))),
           "regions": []}
    for label, mask in [("all", np.ones(len(s), dtype=bool)),
                        ("affine_core_r_lt_015", radius < 0.15),
                        ("smooth_annulus_025_to_035", (radius > 0.25) & (radius < 0.35)),
                        ("outer_r_gt_045", radius > 0.45)]:
        v = volume[mask]

        def average(value):
            return float(np.sum(v * value[mask]) / np.sum(v))

        grad_error = np.hypot(s["pressure_gradient_error_x"], s["pressure_gradient_error_y"])
        pressure_noise = np.hypot(s["pressure_noise_force_x"], s["pressure_noise_force_y"])
        total_error = np.hypot(s["pressure_gradient_error_x"] + s["pressure_noise_force_x"],
                               s["pressure_gradient_error_y"] + s["pressure_noise_force_y"])
        observed_gap = s["div_geometric"] - s["div_fit"]
        affine_gap = s["div_geometric_affine_fit"] - s["div_fit"]
        remainder = s["div_geometric"] - s["div_geometric_affine_fit"]
        region = {"name": label, "cells": int(mask.sum()), "volume": float(v.sum()),
                  "density_L1_from_one": average(abs(s["m"] / volume - 1)),
                  "pressure_L1_from_exact_nodal_profile": average(abs(s["P"] - s["P_exact_node"])),
                  "actual_geometric_divergence_L1": average(abs(s["div_geometric"])),
                  "actual_fit_divergence_L1": average(abs(s["div_fit"])),
                  "actual_divergence_gap_L1": average(abs(observed_gap)),
                  "affine_geometry_gap_L1": average(abs(affine_gap)),
                  "curved_velocity_remainder_L1": average(abs(remainder)),
                  "divergence_decomposition_max_error": float(np.max(abs(observed_gap - affine_gap - remainder))),
                  "exact_velocity_geometric_divergence_L1": average(abs(s["div_geometric_exact_velocity"])),
                  "exact_velocity_fit_divergence_L1": average(abs(s["div_fit_exact_velocity"])),
                  "exact_analytic_gradient_geometric_divergence_L1": average(abs(s["div_geometric_exact_gradient"])),
                  "pressure_discretization_force_density_L1": average(grad_error),
                  "evolved_pressure_noise_force_density_L1": average(pressure_noise),
                  "total_pressure_force_density_error_L1": average(total_error),
                  "instantaneous_heat_rate": float(heat["heat_rate"][mask].sum())}
        row["regions"].append(region)
    rows.append(row)
    all_region = row["regions"][0]
    print(row["case"], json.dumps({key: all_region[key] for key in (
        "actual_divergence_gap_L1", "exact_velocity_geometric_divergence_L1",
        "pressure_discretization_force_density_L1", "evolved_pressure_noise_force_density_L1")}), flush=True)
with args.output.open("x") as output:
    json.dump({"scope": "Frozen exact geometric work and force decomposition. Exact pressure is a generator sample, not a cell integral. Norm components need not add; vector decompositions do.",
               "cases": rows}, output, indent=2)
    output.write("\n")
