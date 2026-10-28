"""Audit the energy term required by time-dependent power weights.

The existing W=0 variational pressure force is the position adjoint of
geometric cell-volume change. A changing power weight adds J Wdot to that
volume rate without changing the position derivative of the pressure force.
This small periodic mesh makes the resulting power identity explicit.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_kepler_varpress_periodic_adjoint_20261027 import (
    BOX, GRID, source_face_sums, tessellate,
)


def power_areas(base, weights):
    """Central-cell areas from direct clipping against 3x3 periodic images."""
    count = len(base)
    shifts = np.array([(ix * BOX, iy * BOX)
                       for iy in (-1, 0, 1) for ix in (-1, 0, 1)])
    images = (base[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    image_weights = np.tile(weights, 9)
    areas = np.empty(count)
    for i in range(count):
        center = 4 * count + i
        xi = images[center]
        polygon = np.array([xi + [-BOX, -BOX], xi + [BOX, -BOX],
                            xi + [BOX, BOX], xi + [-BOX, BOX]])
        for j, xj in enumerate(images):
            if j == center:
                continue
            normal = 2.0 * (xj - xi)
            rhs = (np.dot(xj, xj) - np.dot(xi, xi)
                   + image_weights[center] - image_weights[j])
            distance = polygon @ normal - rhs
            if np.max(distance) <= 0:
                continue
            clipped = []
            for k, end in enumerate(polygon):
                start = polygon[k - 1]
                ds, de = distance[k - 1], distance[k]
                if (ds <= 0) != (de <= 0):
                    clipped.append(start + (end - start) * ds / (ds - de))
                if de <= 0:
                    clipped.append(end)
            polygon = np.asarray(clipped)
            if len(polygon) < 3:
                raise ValueError(f"Empty power cell {i} at image {j}")
        areas[i] = 0.5 * abs(
            np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
            - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
    return areas


def jacobian_from_faces(vor, images, count):
    jacobian = np.zeros((count, count))
    for (a, b), vertices in zip(vor.ridge_points, vor.ridge_vertices):
        for center, neighbor in ((a, b), (b, a)):
            if not (4 * count <= center < 5 * count):
                continue
            if len(vertices) != 2 or min(vertices) < 0:
                raise ValueError("Unbounded central face")
            length = np.linalg.norm(np.diff(vor.vertices[vertices], axis=0))
            distance = np.linalg.norm(images[neighbor] - images[center])
            conductance = float(length / (2.0 * distance))
            i, j = center - 4 * count, neighbor % count
            jacobian[i, i] += conductance
            jacobian[i, j] -= conductance
    return jacobian


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(17)
    base = np.array([(x + 0.5, y + 0.5)
                     for y in range(GRID) for x in range(GRID)], dtype=float)
    base += rng.uniform(-0.12, 0.12, size=base.shape)
    velocity = rng.normal(scale=0.2, size=base.shape)
    pressure = 1.0 + 0.25 * rng.normal(size=len(base))
    vor, images, voronoi_area = tessellate(base)
    zero = np.zeros(len(base))
    power_area = power_areas(base, zero)
    jacobian = jacobian_from_faces(vor, images, len(base))
    fd_jacobian = np.empty_like(jacobian)
    h = 1.0e-4
    for j in range(len(base)):
        dw = np.zeros(len(base))
        dw[j] = h
        fd_jacobian[:, j] = (power_areas(base, dw)
                             - power_areas(base, -dw)) / (2.0 * h)
    force = np.empty_like(base)
    geom_volume_rate = np.empty(len(base))
    thermal_rate = np.empty(len(base))
    for i in range(len(base)):
        force[i], _, geom_volume_rate[i], thermal_rate[i], _ = \
            source_face_sums(vor, images, velocity, pressure, i)
    rate = np.linalg.lstsq(jacobian, -geom_volume_rate, rcond=1e-12)[0]
    rate -= rate.mean()
    weighted_volume_rate = geom_volume_rate + jacobian @ rate
    kinetic_pressure_power = float(np.sum(force * velocity))
    unweighted_thermal_power = float(np.sum(thermal_rate))
    weighted_thermal_power = float(-pressure @ weighted_volume_rate)
    weight_power = float(pressure @ (jacobian @ rate))
    result = {
        "reproducer": str(Path(__file__).resolve()),
        "reproducer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference": "W0 variational force in Exam/exam_gpu.cu; J_ij from power-face conductance",
        "grid": [GRID, GRID], "box": BOX, "sites": len(base), "h": h,
        "voronoi_power_area_max_difference": float(np.max(np.abs(voronoi_area - power_area))),
        "jacobian_finite_difference_max_error": float(np.max(np.abs(jacobian - fd_jacobian))),
        "jacobian_row_sum_max": float(np.max(np.abs(jacobian.sum(axis=1)))),
        "jacobian_column_sum_max": float(np.max(np.abs(jacobian.sum(axis=0)))),
        "unweighted_area_rate_sum": float(geom_volume_rate.sum()),
        "compensated_area_rate_max": float(np.max(np.abs(weighted_volume_rate))),
        "kinetic_pressure_power": kinetic_pressure_power,
        "unweighted_thermal_power": unweighted_thermal_power,
        "weighted_thermal_power": weighted_thermal_power,
        "weight_power": weight_power,
        "unweighted_total_power": kinetic_pressure_power + unweighted_thermal_power,
        "weighted_total_power": kinetic_pressure_power + weighted_thermal_power,
        "weighted_power_identity_error": abs(
            kinetic_pressure_power + weighted_thermal_power + weight_power),
        "scope": "Small periodic smooth mesh at W=0 with nonzero Wdot; first-order weighted-volume audit, not native CUDA or finite-W evolution. Weight power must be represented by another conservative degree of freedom or a modified force/work law before a weighted fixed-mass hydro claim.",
    }
    checks = {
        "power_reduces_to_voronoi": result["voronoi_power_area_max_difference"] < 1e-10,
        "weight_jacobian_matches_power_geometry": result["jacobian_finite_difference_max_error"] < 1e-7,
        "jacobian_has_constant_gauge": result["jacobian_row_sum_max"] < 1e-10,
        "jacobian_preserves_total_area": result["jacobian_column_sum_max"] < 1e-10,
        "compensated_area_rate_zero": result["compensated_area_rate_max"] < 1e-10,
        "unweighted_energy_conserved": abs(result["unweighted_total_power"]) < 1e-10,
        "weighted_power_identity": result["weighted_power_identity_error"] < 1e-10,
        "weight_power_is_nontrivial": abs(weight_power) > 1e-4,
    }
    result["checks"] = checks
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    if not all(checks.values()):
        raise SystemExit("Weight-rate energy identity gate failed")


if __name__ == "__main__":
    main()
