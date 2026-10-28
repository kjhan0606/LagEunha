"""Gate a sparse native Kepler IC against complete periodic t=0 geometry."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import struct
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def low_index(p):
    return (p["u4if"].astype(np.uint64) & np.uint64((1 << 56) - 1)).astype(np.int64)


def read_glass(path):
    with path.open("rb") as stream:
        n, lx, ly = struct.unpack("<idd", stream.read(20))
        xy = np.column_stack((np.fromfile(stream, dtype="<f8", count=n),
                              np.fromfile(stream, dtype="<f8", count=n)))
    if (lx, ly) != (6., 6.) or len(xy) != n or path.stat().st_size != 20 + 16 * n:
        raise ValueError("Invalid native glass")
    return xy


def read_areas(path, n):
    with path.open("rb") as stream:
        count, lx, ly = struct.unpack("<idd", stream.read(20))
        area = np.fromfile(stream, dtype="<f8", count=count)
    if (count, lx, ly) != (n, 6., 6.) or len(area) != n or path.stat().st_size != 20 + 8 * n:
        raise ValueError("Invalid periodic area sidecar")
    return area


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expected-count", type=int, required=True)
    ap.add_argument("--glass-sha256", required=True)
    ap.add_argument("--area-sha256", required=True)
    ap.add_argument("--area-file", type=Path, default=Path("periodic_areas.bin"))
    ap.add_argument("--output", type=Path, default=Path("audit_t001.json"))
    args = ap.parse_args()
    run = Path.cwd()
    snapshots = sorted(run.glob("khout.*.dat"))
    if len(snapshots) != 2 or snapshots[0].name != "khout.000000.dat":
        raise ValueError("Expected exactly initial and one t=.001 snapshot")
    p0, meta0 = read_snapshot(snapshots[0])
    p1, meta1 = read_snapshot(snapshots[1])
    glass = run / "glass.bin"
    sidecar = args.area_file if args.area_file.is_absolute() else run / args.area_file
    if digest(glass) != args.glass_sha256 or digest(sidecar) != args.area_sha256:
        raise ValueError("Glass or periodic area checksum differs from approved IC")
    xy = read_glass(glass)
    area = read_areas(sidecar, args.expected_count)
    n = args.expected_count
    if len(xy) != n or len(p0) != n or len(p1) != n:
        raise ValueError("Sparse site count mismatch")
    if abs(meta0["t"]) > 1e-12 or abs(meta1["t"] - .001) > 1e-12:
        raise ValueError("Unexpected native times")
    expected_ids = np.arange(n)
    if not (np.array_equal(np.sort(low_index(p0)), expected_ids) and
            np.array_equal(np.sort(low_index(p1)), expected_ids)):
        raise ValueError("Incomplete native cell IDs")
    p0 = p0[np.argsort(low_index(p0))]
    p1 = p1[np.argsort(low_index(p1))]
    radius = np.hypot(xy[:, 0] - 3., xy[:, 1] - 3.)
    profile = np.ones_like(radius)
    inner = radius < .5
    outer = radius > 2.
    profile[inner] = (radius[inner] / .5) ** 3
    profile[outer] = (1 + (radius[outer] - 2) / .1) ** -3
    rho_target = .01 + profile
    expected_mass = rho_target * area
    disk = (radius >= .5) & (radius < 2.)
    phase = .001 * (radius * radius + .01175 ** 2) ** (-.75)
    dx, dy = xy[:, 0] - 3., xy[:, 1] - 3.
    orbit_x = 3 + dx * np.cos(phase) - dy * np.sin(phase)
    orbit_y = 3 + dx * np.sin(phase) + dy * np.cos(phase)
    orbit_error = np.hypot(p1["x"] - orbit_x, p1["y"] - orbit_y)
    p0_area = p0["volume"].astype(float)
    p1_area = p1["volume"].astype(float)
    mass = p0["mass"].astype(float)
    log = (run / "log").read_text(errors="replace")
    energy_lines = [line for line in log.splitlines() if "[RK4E]" in line]
    zero_injection = bool(energy_lines)
    if energy_lines:
        for key in ("floor_cum", "de_cum", "sfl_cum"):
            match = re.search(rf"\b{key}=([0-9.eE+-]+)", energy_lines[-1])
            zero_injection &= bool(match and float(match.group(1)) == 0.)
    finite = all(np.isfinite(q[field]).all() for q in (p0, p1) for field in
                 ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2"))
    checks = {
        "finite": bool(finite),
        "periodic_area_sidecar_tiles_box": bool(abs(np.sum(area, dtype=np.longdouble) - 36.) < 1e-8),
        "initial_native_areas_match_sidecar": bool(np.max(abs(p0_area / area - 1)) < 1e-6),
        "initial_native_areas_tile_box": bool(abs(np.sum(p0_area, dtype=np.longdouble) - 36.) < 1e-6),
        "evolved_native_areas_tile_box": bool(abs(np.sum(p1_area, dtype=np.longdouble) - 36.) < 1e-6),
        "initial_mass_matches_target_profile_and_sidecar": bool(np.max(abs(mass / expected_mass - 1)) < 1e-6),
        "initial_total_mass_matches_profile": bool(abs(np.sum(mass) / np.sum(expected_mass) - 1) < 1e-8),
        "mass_exact_across_step": bool(np.array_equal(p0["mass"], p1["mass"])),
        "W0": bool(np.all(p0["w2"] == 0) and np.all(p1["w2"] == 0)),
        "positive_state": bool(np.all(p0_area > 0) and np.all(p1_area > 0) and
                               np.all(p0["ie"] > 0) and np.all(p1["ie"] > 0)),
        "no_floor_injection": zero_injection,
        "glass_initial_positions": bool(np.max(np.hypot(p0["x"] - xy[:, 0], p0["y"] - xy[:, 1])) < 1e-12),
        "disk_pressure": bool(np.max(abs(p0["pressure"][disk] / 1e-6 - 1)) < 1e-6),
        "disk_orbit": bool(np.max(orbit_error[disk]) < 1e-6),
    }
    result = {
        "initial": str(snapshots[0]), "evolved": str(snapshots[1]),
        "t0": meta0["t"], "t1": meta1["t"], "n": n,
        "sha256": {p.name: digest(p) for p in (glass, sidecar, snapshots[0], snapshots[1], run / "eunha2.exe")},
        "initial_volume_sum": float(np.sum(p0_area, dtype=np.longdouble)),
        "evolved_volume_sum": float(np.sum(p1_area, dtype=np.longdouble)),
        "initial_volume_max": float(np.max(p0_area)),
        "expected_total_mass": float(np.sum(expected_mass, dtype=np.longdouble)),
        "initial_total_mass": float(np.sum(mass, dtype=np.longdouble)),
        "initial_mass_max": float(np.max(mass)),
        "disk_mass": float(np.sum(mass[disk], dtype=np.longdouble)),
        "initial_area_relative_error_max": float(np.max(abs(p0_area / area - 1))),
        "initial_mass_relative_error_max": float(np.max(abs(mass / expected_mass - 1))),
        "disk_orbit_error_max": float(np.max(orbit_error[disk])),
        "checks": checks,
        "scope": "Native t=.001 sparse W0 CUDA+MPI+OpenMP IC gate with independently certified periodic t=0 volumes and fixed mass; not a long-time result",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Sparse periodic mass preflight failed")


if __name__ == "__main__":
    main()
