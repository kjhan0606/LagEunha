"""Independent periodic Voronoi audit of the active smooth W=0 pressure law.

Rebuild a 3-by-3 tiling for each finite difference. The reference energy uses
only central-image cell areas, while each base generator perturbation moves all
of its periodic images. This tests seam faces and global adjoint identities in
a small mathematical mesh, not the native CUDA tessellation or long evolution.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi


GRID = 5
BOX = float(GRID)
GAMMA = 5.0 / 3.0
H = 1.0e-4
SEEDS = (17, 23, 41)
SOURCE = Path(__file__).with_name("KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27") / "Exam" / "exam_gpu.cu"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tessellate(base):
    shifts = np.array([(ix * BOX, iy * BOX)
                       for iy in (-1, 0, 1) for ix in (-1, 0, 1)])
    images = (base[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    vor = Voronoi(images, qhull_options="Qbb Qc Qz")
    first = 4 * len(base)
    areas = np.empty(len(base))
    for i in range(len(base)):
        region = vor.regions[vor.point_region[first + i]]
        if not region or -1 in region:
            raise ValueError("Central periodic cell is unbounded")
        polygon = vor.vertices[region]
        areas[i] = 0.5 * abs(
            np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
            - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
    if np.any(areas <= 0):
        raise ValueError("Nonpositive periodic area")
    return vor, images, areas


def energy(area, entropy_constant):
    return float(np.sum(entropy_constant * area ** (1.0 - GAMMA),
                        dtype=np.longdouble) / (GAMMA - 1.0))


def source_face_sums(vor, images, base_velocity, base_pressure, site):
    count = len(base_velocity)
    center = 4 * count + site
    velocity = np.tile(base_velocity, (9, 1))
    pressure = np.tile(base_pressure, 9)
    force = np.zeros(2)
    area_vector_sum = np.zeros(2)
    area_rate = 0.0
    thermal_rate = 0.0
    faces = 0
    for (a, b), vertices in zip(vor.ridge_points, vor.ridge_vertices):
        if center not in (a, b):
            continue
        neighbor = b if a == center else a
        if len(vertices) != 2 or min(vertices) < 0:
            raise ValueError("Central periodic face is unbounded")
        corners = vor.vertices[vertices]
        length = float(np.linalg.norm(corners[1] - corners[0]))
        delta = images[neighbor] - images[center]
        distance = float(np.linalg.norm(delta))
        normal = delta / distance
        area_vector = length * normal
        offset = 0.5 * (corners[0] + corners[1] - images[center] - images[neighbor])
        b_i = 0.5 * area_vector + length * offset / distance
        force += (pressure[center] - pressure[neighbor]) * b_i
        area_vector_sum += area_vector
        geometric_speed = (
            np.dot(0.5 * (velocity[center] + velocity[neighbor]), normal)
            - np.dot(velocity[neighbor] - velocity[center], offset) / distance)
        area_rate += length * geometric_speed
        thermal_rate -= pressure[center] * length * (
            geometric_speed - np.dot(velocity[center], normal))
        faces += 1
    if faces < 3:
        raise ValueError("Incomplete central polygon")
    return force, area_vector_sum, area_rate, thermal_rate, faces


def run_case(seed, uniform_pressure):
    rng = np.random.default_rng(seed)
    base = np.array([(x + 0.5, y + 0.5)
                     for y in range(GRID) for x in range(GRID)], dtype=float)
    base += rng.uniform(-0.12, 0.12, size=base.shape)
    velocity = rng.normal(scale=0.2, size=base.shape)
    pressure = (np.ones(len(base)) if uniform_pressure else
                1.0 + 0.25 * rng.normal(size=len(base)))
    if np.any(pressure <= 0):
        raise ValueError("Nonpositive test pressure")
    vor, images, area = tessellate(base)
    entropy_constant = pressure * area ** GAMMA
    force = np.empty_like(base)
    face_area_rate = np.empty(len(base))
    thermal_rate = np.empty(len(base))
    closure = np.empty(len(base))
    face_count = np.empty(len(base), dtype=int)
    for i in range(len(base)):
        force[i], area_sum, face_area_rate[i], thermal_rate[i], face_count[i] = \
            source_face_sums(vor, images, velocity, pressure, i)
        closure[i] = np.linalg.norm(area_sum)

    drift_plus = tessellate(base + H * velocity)[2]
    drift_minus = tessellate(base - H * velocity)[2]
    fd_area_rate = (drift_plus - drift_minus) / (2.0 * H)
    fd_force = np.empty_like(base)
    for i in range(len(base)):
        for axis in range(2):
            delta = np.zeros_like(base)
            delta[i, axis] = H
            ep = energy(tessellate(base + delta)[2], entropy_constant)
            em = energy(tessellate(base - delta)[2], entropy_constant)
            fd_force[i, axis] = -(ep - em) / (2.0 * H)

    return {
        "seed": seed, "uniform_pressure": uniform_pressure,
        "sites": len(base), "faces_min": int(face_count.min()),
        "periodic_area_error": float(abs(area.sum() - BOX ** 2)),
        "face_vector_closure_max": float(closure.max()),
        "force_gradient_error_max": float(np.linalg.norm(force - fd_force, axis=1).max()),
        "area_rate_error_max": float(np.abs(face_area_rate - fd_area_rate).max()),
        "thermal_work_error_max": float(np.abs(thermal_rate + pressure * fd_area_rate).max()),
        "total_area_rate_error": float(abs(face_area_rate.sum())),
        "total_force_error": float(np.linalg.norm(force.sum(axis=0))),
        "total_energy_power_error": float(abs(np.sum(force * velocity) + thermal_rate.sum())),
        "force_max": float(np.linalg.norm(force, axis=1).max()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = [run_case(seed, uniform)
             for seed in SEEDS for uniform in (False, True)]
    maxima = {
        key: max(case[key] for case in cases)
        for key in ("periodic_area_error", "face_vector_closure_max",
                    "force_gradient_error_max", "area_rate_error_max",
                    "thermal_work_error_max", "total_area_rate_error",
                    "total_force_error", "total_energy_power_error")
    }
    checks = {
        "periodic_area_partition": maxima["periodic_area_error"] < 1e-10,
        "closed_cell_faces": maxima["face_vector_closure_max"] < 1e-10,
        "force_is_negative_energy_gradient": maxima["force_gradient_error_max"] < 1e-7,
        "face_rate_matches_area_derivative": maxima["area_rate_error_max"] < 1e-7,
        "thermal_work_is_negative_P_dVdt": maxima["thermal_work_error_max"] < 1e-7,
        "total_area_rate_zero": maxima["total_area_rate_error"] < 1e-10,
        "total_force_zero": maxima["total_force_error"] < 1e-10,
        "total_energy_power_zero": maxima["total_energy_power_error"] < 1e-10,
        "uniform_pressure_force_zero": all(
            case["force_max"] < 1e-10 for case in cases if case["uniform_pressure"]),
    }
    result = {
        "source": str(SOURCE), "source_sha256": sha256(SOURCE),
        "reproducer": str(Path(__file__).resolve()),
        "reproducer_sha256": sha256(Path(__file__)),
        "grid": [GRID, GRID], "box": BOX, "gamma": GAMMA,
        "central_difference_h": H, "seeds": list(SEEDS),
        "cases": cases, "maxima": maxima, "checks": checks,
        "scope": "Small periodic smooth W=0 mathematical operator audit; no native CUDA execution, GPU seam-pair extraction, shocks, or long-time Kepler inference.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"maxima": maxima, "checks": checks}, indent=2))
    if not all(checks.values()):
        raise SystemExit("Periodic variational pressure adjoint audit failed")


if __name__ == "__main__":
    main()
