"""Gate a native weighted Kepler continuation against the disk equilibrium.

An exact-orbit W0 Voronoi area is a useful geometry counterfactual, but it is
not the target area for an evolving power diagram. The physical density gate
therefore compares the final density with the initial equilibrium density.
Passing still does not establish local moving-face mass/energy conservation.
"""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_bounded_native_evolution_20261027 import (
    APOTHEM, DOMAIN_AREA, GAMMA_MINUS_ONE, N, energy, weighted_abs,
)
from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from audit_kepler_sparse_state_consistency_20261027 import state_metrics
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


def unique_power_edge_ratio(triangles, xy, weight):
    edges = np.concatenate((triangles[:, (0, 1)], triangles[:, (1, 2)],
                            triangles[:, (2, 0)]))
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    delta = xy[edges[:, 0]] - xy[edges[:, 1]]
    d2 = np.einsum("ij,ij->i", delta, delta)
    if not np.all(d2 > 0):
        raise ValueError("Duplicate sites in power edges")
    return float(np.max(np.abs(weight[edges[:, 0]]-weight[edges[:, 1]])/d2)), len(edges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("initial", type=Path)
    ap.add_argument("parent", type=Path)
    ap.add_argument("run", type=Path)
    ap.add_argument("--time", type=float, required=True)
    ap.add_argument("--parent-time", type=float, default=6.501)
    ap.add_argument("--expected-binary-sha256", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if (abs(snapshot_time(args.initial, N)) > 1e-12
            or abs(snapshot_time(args.parent, N)-args.parent_time) > 1e-10):
        raise ValueError("Expected t=0 reference and requested weighted parent time")
    restart = args.run / "khout.000001.dat"
    if (abs(snapshot_time(restart, N)-args.parent_time) > 1e-10
            or digest(restart) != digest(args.parent)
            or digest(args.run/"initial.dat") != digest(args.initial)):
        raise ValueError("Run lineage differs from the supplied initial/parent state")
    binary_sha256 = digest(args.run/"eunha2.exe")
    if binary_sha256 != args.expected_binary_sha256:
        raise ValueError("Run binary differs from the expected weighted build")
    candidates = [p for p in args.run.glob("khout.*.dat")
                  if abs(snapshot_time(p, N)-args.time)<1e-10]
    if len(candidates) != 1:
        raise ValueError(f"Expected one complete t={args.time} endpoint")
    endpoint = candidates[0]
    initial, _ = ordered(args.initial, N)
    parent, _ = ordered(args.parent, N)
    final, _ = ordered(endpoint, N)
    s0, s1 = state_metrics(initial), state_metrics(final)
    x0 = initial["x"].astype(float)-3.
    y0 = initial["y"].astype(float)-3.
    r0 = np.hypot(x0, y0)
    disk = (r0>=.5)&(r0<2.)
    ids = np.flatnonzero(disk)
    masses = initial["mass"][ids].astype(float)
    exact_xy, ballistic_area, image_safety = disk_geometry(initial, args.time, ids)
    volume = final["volume"][ids].astype(float)
    density_change = final["den"][ids].astype(float)/initial["den"][ids].astype(float)-1.
    w0_geometry_difference = ballistic_area/volume-1.
    position_error = np.hypot(final["x"][ids].astype(float)-exact_xy[ids, 0],
                              final["y"][ids].astype(float)-exact_xy[ids, 1])
    entropy_log = (
        np.log(final["ie"][ids].astype(float)/initial["ie"][ids].astype(float))
        + GAMMA_MINUS_ONE*np.log(volume/initial["volume"][ids].astype(float))
    )
    x = final["x"].astype(float)-3.
    y = final["y"].astype(float)-3.
    r = np.hypot(x, y)
    vx = final["vx"].astype(float)
    vy = final["vy"].astype(float)
    radial_speed = (x*vx+y*vy)/np.maximum(r, 1e-30)
    circular_speed = r0*(r0*r0+.01175**2)**(-.75)
    radial_rms = float(np.sqrt(np.average(
        (radial_speed[ids]/circular_speed[ids])**2, weights=masses)))
    l0 = float(np.sum(masses*(x0[ids]*initial["vy"][ids].astype(float)
                              - y0[ids]*initial["vx"][ids].astype(float)),
                      dtype=np.longdouble))
    l1 = float(np.sum(masses*(x[ids]*vy[ids]-y[ids]*vx[ids]), dtype=np.longdouble))
    energy_relative = abs(energy(final)-energy(initial))/abs(energy(initial))
    weight = final["w2"].astype(float)
    xy = np.column_stack((x, y))
    lookup = np.full(N, -1, dtype=int)
    lookup[ids] = np.arange(len(ids))
    power_area, power_meta, triangles, _ = power_areas(
        xy, weight, ids, lookup, return_geometry=True)
    if power_meta["hidden_active"] or not np.all(power_area > 0):
        raise ValueError("Hidden or nonpositive disk power cell")
    area_error = np.abs(volume/power_area-1.)
    face_ratio, unique_edges = unique_power_edge_ratio(triangles, xy, weight)
    log_lines = (args.run/"log").read_text(errors="replace").splitlines()
    endpoint_step = int(endpoint.name.split(".")[1])
    endpoint_marker = next((i for i, line in enumerate(log_lines)
                            if line.startswith("Time=") and
                            re.search(rf"\bstep=\s*{endpoint_step}\b", line)), None)
    if endpoint_marker is None:
        raise ValueError("Missing complete log marker for the endpoint step")
    log = "\n".join(log_lines[:endpoint_marker+1])
    ledgers = [line for line in log.splitlines() if "[RK4E]" in line]
    if not ledgers:
        raise ValueError("Missing native RK4 energy ledger")
    injection = {}
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", ledgers[-1])
        if not match:
            raise ValueError(f"Missing {key} energy ledger")
        injection[key] = float(match.group(1))
    stage_pattern = re.compile(
        r"\[GPU_VOLW2\] stage=(\d+) active=(\d+) faces=(\d+) "
        r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=([\deE+.-]+) s")
    stages = [stage_pattern.fullmatch(line) for line in log.splitlines()
              if line.startswith("[GPU_VOLW2] stage=")]
    if any(match is None for match in stages):
        raise ValueError("Malformed weighted GPU stage receipt")
    max_work = max((abs(float(match.group(5))) for match in stages), default=float("inf"))
    max_stage_seconds = max((float(match.group(6)) for match in stages), default=float("inf"))
    ordered_stages = bool(stages and len(stages)%4 == 0 and
                          all(int(match.group(1)) == i%4 and int(match.group(2)) == 49038
                              for i, match in enumerate(stages)))
    disk_density = weighted_abs(density_change, masses)
    disk_w0_geometry_difference = weighted_abs(w0_geometry_difference, masses)
    disk_entropy = weighted_abs(entropy_log, masses)
    disk_position = weighted_abs(position_error, masses)
    angular_drift = l1/l0-1.
    checks = {
        "finite_positive": bool(s0["positive"] and s1["positive"]),
        "fixed_mass_exact": bool(np.array_equal(initial["mass"], final["mass"])),
        "initial_W0_and_evolved_weights": bool(np.all(initial["w2"] == 0)
                                              and np.all(weight >= 0)
                                              and np.any(weight[ids] > 0)),
        "area_tiles_polygon": bool(abs(s0["volume_sum"]-DOMAIN_AREA)<1e-6
                                   and abs(s1["volume_sum"]-DOMAIN_AREA)<1e-5),
        "independent_power_area": bool(np.max(area_error)<1e-4),
        "weight_face_guard": bool(face_ratio<.8),
        "density_and_EOS": bool(max(s0["density_relative_max"],
                                    s0["eos_pressure_relative_max"],
                                    s1["density_relative_max"],
                                    s1["eos_pressure_relative_max"])<1e-6),
        "zero_injection": bool(all(value==0 for value in injection.values())),
        "total_energy": bool(energy_relative<1e-6),
        "disk_density_broad_safety": bool(disk_density<.05),
        "disk_density_equilibrium": bool(disk_density<1e-3),
        "disk_entropy": bool(disk_entropy<.01),
        "disk_orbit": bool(disk_position<.01),
        "disk_radial_velocity": bool(radial_rms<.05),
        "disk_angular_momentum": bool(abs(angular_drift)<.01),
        "inside_polygon_incircle": bool(np.all(r<APOTHEM)),
        "projected_weight_work": bool(ordered_stages and max_work<1e-10),
        "weighted_pressure_force_active": bool(
            "[VARPRESS_W0] gated weighted fixed-mass variational pressure active" in log
            or "[VARPRESS_W0] experimental smooth fixed-mass variational pressure active" in log),
        "no_geometry_abort": bool("[GPU_TESS_INVALID]" not in log and
                                  "[GPU_CLIP_CAPACITY]" not in log),
        "no_wall_crossing": bool("[FIXED_POLYGON_CROSSING]" not in log),
        "fixed_polygon_enabled": bool("[FIXED_POLYGON] enabled" in log),
    }
    result = {
        "initial": str(args.initial.resolve()), "parent": str(args.parent.resolve()),
        "endpoint": str(endpoint.resolve()), "initial_sha256": digest(args.initial),
        "parent_sha256": digest(args.parent), "restart_sha256": digest(restart),
        "binary_sha256": binary_sha256, "endpoint_sha256": digest(endpoint),
        "time": args.time, "parent_time": args.parent_time,
        "n": N, "disk_sites": len(ids),
        "initial_state": s0, "endpoint_state": s1,
        "total_energy_relative": energy_relative,
        "disk_density_abs_mass_weighted": disk_density,
        "disk_weighted_minus_W0_ballistic_geometry_abs_mass_weighted":
            disk_w0_geometry_difference,
        "disk_density_equilibrium_threshold": 1e-3,
        "disk_entropy_log_abs_mass_weighted": disk_entropy,
        "disk_position_error_abs_mass_weighted": disk_position,
        "disk_radial_speed_rms_over_circular": radial_rms,
        "disk_angular_momentum_relative": angular_drift,
        "independent_power_area_error_max": float(np.max(area_error)),
        "independent_power_area_error_p99": float(np.quantile(area_error, .99)),
        "unique_power_edges": unique_edges, "max_abs_weight_difference_over_d2": face_ratio,
        "weight_min": float(weight.min()), "weight_max": float(weight.max()),
        "gpu_weight_stage_count": len(stages), "gpu_weight_max_abs_stage_work": max_work,
        "gpu_weight_max_stage_seconds": max_stage_seconds,
        "rk4_injection": injection, "minimum_ballistic_image_safety_ratio": image_safety,
        "checks": checks,
        "scope": "Experimental fixed-mass weighted RK4 endpoint. Equilibrium density, "
                 "circular orbit and independent power geometry are gates. The exact-orbit "
                 "W0 Voronoi area difference is diagnostic only: a weighted cell is not "
                 "required to equal that area. Local moving-face mass/energy consistency "
                 "and GIZMO t=600 stability remain unproven.",
    }
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    if not all(checks.values()):
        raise SystemExit("Native weighted long-endpoint physics gate failed")


if __name__ == "__main__":
    main()
