"""Independent small-mesh finite-difference audit of the active W0 pressure law.

SciPy constructs Voronoi polygons. The tested force formula is copied from
the active CUDA branch, but the reference force comes from independently
rebuilding all polygon areas after perturbing one generator. This checks an
interior geometric identity, not GPU face extraction or long-time dynamics.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi


GAMMA = 5.0 / 3.0
GRID = 7
SITES = ((2, 2), (3, 2), (3, 3), (2, 4), (4, 4))
SEEDS = (17, 23, 41)
H = 1.0e-4
QHULL_OPTIONS = "Qbb Qc Qz"
SOURCE = Path(__file__).with_name("KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27") / "Exam" / "exam_gpu.cu"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_native(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported JSON value: {type(value).__name__}")


def diagram(points):
    vor = Voronoi(points, qhull_options=QHULL_OPTIONS)
    area = np.full(len(points), np.nan)
    for k in range(len(points)):
        region = vor.regions[vor.point_region[k]]
        if region and -1 not in region:
            polygon = vor.vertices[region]
            area[k] = 0.5 * abs(
                np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
                - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
            )
    return vor, area


def internal_energy(area, bounded, entropy_constant):
    if not np.array_equal(np.isfinite(area), bounded):
        raise ValueError("Bounded-cell set changed under a small perturbation")
    if np.any(area[bounded] <= 0):
        raise ValueError("Nonpositive bounded Voronoi area")
    value = entropy_constant * area[bounded] ** (1.0 - GAMMA)
    return float(np.sum(value, dtype=np.longdouble) / (GAMMA - 1.0))


def source_face_sums(vor, points, velocity, pressure, site):
    force = np.zeros(2)
    outward_area_sum = np.zeros(2)
    geometric_area_rate = 0.0
    source_thermal_rate = 0.0
    faces = 0
    for (a, b), vertices in zip(vor.ridge_points, vor.ridge_vertices):
        if site not in (a, b):
            continue
        neighbor = b if a == site else a
        if len(vertices) != 2 or min(vertices) < 0:
            raise ValueError("Test site has an unbounded face")
        corners = vor.vertices[vertices]
        length = float(np.linalg.norm(corners[1] - corners[0]))
        delta = points[neighbor] - points[site]
        distance = float(np.linalg.norm(delta))
        if not (length > 0 and distance > 0):
            raise ValueError("Degenerate face")
        normal = delta / distance
        area_vector = length * normal
        centroid = 0.5 * (corners[0] + corners[1])
        offset = centroid - 0.5 * (points[site] + points[neighbor])
        b_i = 0.5 * area_vector + length * offset / distance
        force += (pressure[site] - pressure[neighbor]) * b_i
        outward_area_sum += area_vector
        geometric_speed = (
            np.dot(0.5 * (velocity[site] + velocity[neighbor]), normal)
            - np.dot(velocity[neighbor] - velocity[site], offset) / distance
        )
        geometric_area_rate += length * geometric_speed
        source_thermal_rate -= pressure[site] * length * (
            geometric_speed - np.dot(velocity[site], normal)
        )
        faces += 1
    if faces < 3:
        raise ValueError("Incomplete interior polygon")
    return force, outward_area_sum, geometric_area_rate, source_thermal_rate, faces


def finite_difference_force(points, bounded, entropy_constant, site):
    force = np.zeros(2)
    for axis in range(2):
        plus = points.copy()
        minus = points.copy()
        plus[site, axis] += H
        minus[site, axis] -= H
        _, area_plus = diagram(plus)
        _, area_minus = diagram(minus)
        up = internal_energy(area_plus, bounded, entropy_constant)
        um = internal_energy(area_minus, bounded, entropy_constant)
        force[axis] = -(up - um) / (2.0 * H)
    return force


def finite_difference_area_rate(points, velocity, site):
    _, plus = diagram(points + H * velocity)
    _, minus = diagram(points - H * velocity)
    if not (np.isfinite(plus[site]) and np.isfinite(minus[site])):
        raise ValueError("Drifted test site became unbounded")
    return float((plus[site] - minus[site]) / (2.0 * H))


def run_case(seed, uniform_pressure):
    rng = np.random.default_rng(seed)
    points = np.array([(x, y) for y in range(GRID) for x in range(GRID)], dtype=float)
    points += rng.uniform(-0.12, 0.12, size=points.shape)
    velocity = rng.normal(scale=0.2, size=points.shape)
    pressure = (np.ones(len(points)) if uniform_pressure
                else 1.0 + 0.25 * rng.normal(size=len(points)))
    if np.any(pressure <= 0):
        raise ValueError("Test pressure must be positive")
    vor, area = diagram(points)
    bounded = np.isfinite(area)
    entropy_constant = pressure[bounded] * area[bounded] ** GAMMA
    rows = []
    for x, y in SITES:
        site = y * GRID + x
        if not bounded[site]:
            raise ValueError("Selected site is unbounded")
        source_force, area_vector_sum, face_rate, thermal_rate, faces = source_face_sums(
            vor, points, velocity, pressure, site
        )
        fd_force = finite_difference_force(points, bounded, entropy_constant, site)
        fd_area_rate = finite_difference_area_rate(points, velocity, site)
        rows.append({
            "site": site,
            "faces": faces,
            "source_force": source_force.tolist(),
            "negative_energy_gradient": fd_force.tolist(),
            "force_absolute_error": float(np.linalg.norm(source_force - fd_force)),
            "area_vector_closure": float(np.linalg.norm(area_vector_sum)),
            "source_geometric_area_rate": face_rate,
            "finite_difference_area_rate": fd_area_rate,
            "area_rate_absolute_error": abs(face_rate - fd_area_rate),
            "thermal_work_absolute_error": abs(thermal_rate + pressure[site] * fd_area_rate),
        })
    return {"seed": seed, "uniform_pressure": uniform_pressure, "rows": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = [run_case(seed, uniform)
             for seed in SEEDS for uniform in (False, True)]
    rows = [row for case in cases for row in case["rows"]]
    maxima = {
        key: max(row[key] for row in rows)
        for key in ("force_absolute_error", "area_vector_closure",
                    "area_rate_absolute_error", "thermal_work_absolute_error")
    }
    checks = {
        "force_equals_negative_energy_gradient": maxima["force_absolute_error"] < 1e-7,
        "closed_polygon_area_vectors": maxima["area_vector_closure"] < 1e-12,
        "face_rate_equals_area_finite_difference": maxima["area_rate_absolute_error"] < 1e-7,
        "thermal_work_equals_negative_P_dVdt": maxima["thermal_work_absolute_error"] < 1e-7,
        "uniform_pressure_force_zero": all(
            float(np.linalg.norm(row["source_force"])) < 1e-12
            for case in cases if case["uniform_pressure"] for row in case["rows"]
        ),
    }
    result = {
        "source": str(SOURCE), "source_sha256": digest(SOURCE),
        "reproducer": str(Path(__file__)), "reproducer_sha256": digest(Path(__file__)),
        "grid": [GRID, GRID], "selected_sites": list(SITES),
        "seeds": list(SEEDS), "central_difference_h": H,
        "gamma": GAMMA, "cases": cases, "maxima": maxima, "checks": checks,
        "scope": "Independent finite-difference check of the smooth interior W=0 CUDA source formula; does not test native GPU polygon pairing, periodic seams, ghost boundaries, shocks, or long-time stability.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, default=json_native)
        stream.write("\n")
    print(json.dumps({"maxima": maxima, "checks": checks}, indent=2,
                     default=json_native))
    if not all(checks.values()):
        raise SystemExit("Variational pressure adjoint audit failed")


if __name__ == "__main__":
    main()
