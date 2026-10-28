"""Finite-weight periodic power-cell pressure/volume/work adjoint audit.

This checks the algebra needed for a fixed-mass weighted RK4 candidate.
The changing weight is only globally work-neutral after projection; the
test does not claim that the face is locally material or that native CUDA
uses the same force and work terms yet.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


GRID = 5
BOX = float(GRID)
GAMMA = 5.0 / 3.0
H = 1e-4


def power_diagram(base, weight, with_faces=False):
    count = len(base)
    shifts = np.array([(ix*BOX, iy*BOX)
                       for iy in (-1, 0, 1) for ix in (-1, 0, 1)])
    image = (base[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    image_w = np.tile(weight, 9)
    areas = np.empty(count)
    face_rows = []
    for i in range(count):
        center = 4*count+i
        xi = image[center]
        polygon = np.array([xi+[-BOX, -BOX], xi+[BOX, -BOX],
                            xi+[BOX, BOX], xi+[-BOX, BOX]])
        for j, xj in enumerate(image):
            if j == center:
                continue
            normal = 2*(xj-xi)
            rhs = (np.dot(xj, xj)-np.dot(xi, xi)
                   +image_w[center]-image_w[j])
            signed = polygon @ normal-rhs
            if np.max(signed) <= 0:
                continue
            clipped = []
            for k, end in enumerate(polygon):
                start = polygon[k-1]
                ds, de = signed[k-1], signed[k]
                if (ds <= 0) != (de <= 0):
                    clipped.append(start+(end-start)*ds/(ds-de))
                if de <= 0:
                    clipped.append(end)
            polygon = np.asarray(clipped)
            if len(polygon) < 3:
                raise ValueError(f"Hidden power cell {i}")
        area = 0.5*(np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
                    -np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1)))
        if not area > 0:
            raise ValueError(f"Nonpositive oriented area {i}")
        areas[i] = area
        if not with_faces:
            continue
        for k, end in enumerate(polygon):
            start = polygon[k-1]
            centroid = 0.5*(start+end)
            length = float(np.linalg.norm(end-start))
            if length < 1e-10:
                continue
            relative = image-centroid
            score = np.abs(np.sum(relative**2, axis=1)-image_w
                           -(np.sum((xi-centroid)**2)-weight[i]))
            score[center] = np.inf
            neighbor = int(np.argmin(score))
            if score[neighbor] > 2e-7:
                raise ValueError(f"Unmatched face of cell {i}: {score[neighbor]}")
            face_rows.append((i, neighbor, length,
                              centroid[0], centroid[1]))
    return areas, image, face_rows


def energy(area, entropy):
    return float(np.sum(entropy*area**(1-GAMMA), dtype=np.longdouble)
                 /(GAMMA-1))


def operators(base, weight, velocity, pressure, rate):
    count = len(base)
    areas, image, faces = power_diagram(base, weight, True)
    image_v = np.tile(velocity, (9, 1))
    image_rate = np.tile(rate, 9)
    image_pressure = np.tile(pressure, 9)
    force = np.zeros((count, 2))
    geom_rate = np.zeros(count)
    weight_rate = np.zeros(count)
    jac = np.zeros((count, count))
    face_count = np.zeros(count, dtype=int)
    for i, j, length, cx, cy in faces:
        centroid = np.array((cx, cy))
        xi, xj = image[4*count+i], image[j]
        delta = xj-xi
        distance = np.linalg.norm(delta)
        if distance <= 0:
            raise ValueError("Coincident generators")
        g = length/(2*distance)
        b_i = length*(centroid-xi)/distance
        force[i] += (pressure[i]-image_pressure[j])*b_i
        geom_rate[i] += length*(np.dot(xj-centroid, image_v[j])
                                -np.dot(xi-centroid, velocity[i]))/distance
        weight_rate[i] += g*(rate[i]-image_rate[j])
        jac[i, i] += g
        jac[i, j%count] -= g
        face_count[i] += 1
    return areas, force, geom_rate, weight_rate, jac, face_count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20261027)
    base = np.array([(x+.5, y+.5)
                     for y in range(GRID) for x in range(GRID)], dtype=float)
    base += rng.uniform(-.12, .12, size=base.shape)
    weight = rng.normal(scale=.03, size=len(base))
    velocity = rng.normal(scale=.2, size=base.shape)
    pressure = 1+.25*rng.normal(size=len(base))
    rate = rng.normal(scale=.05, size=len(base))
    if np.any(pressure <= 0):
        raise ValueError("Nonpositive pressure")
    area, force, geom, weighted, jac, face_count = operators(
        base, weight, velocity, pressure, rate)
    entropy = pressure*area**GAMMA
    fd_combined = (power_diagram(base+H*velocity, weight+H*rate)[0]
                   -power_diagram(base-H*velocity, weight-H*rate)[0])/(2*H)
    fd_weight = (power_diagram(base, weight+H*rate)[0]
                 -power_diagram(base, weight-H*rate)[0])/(2*H)
    fd_force = np.empty_like(base)
    for i in range(len(base)):
        for axis in range(2):
            delta = np.zeros_like(base)
            delta[i, axis] = H
            ep = energy(power_diagram(base+delta, weight)[0], entropy)
            em = energy(power_diagram(base-delta, weight)[0], entropy)
            fd_force[i, axis] = -(ep-em)/(2*H)
    h = jac.T@pressure
    projected = rate-h*(np.dot(h, rate)/np.dot(h, h))
    neutral_weight_rate = jac@projected
    kinetic_power = float(np.sum(force*velocity))
    thermal_power = float(-pressure@(geom+weighted))
    neutral_thermal_power = float(-pressure@(geom+neutral_weight_rate))
    result = {
        "reproducer": str(Path(__file__).resolve()),
        "reproducer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "grid": [GRID, GRID], "box": BOX, "sites": len(base),
        "gamma": GAMMA, "central_difference_h": H,
        "minimum_faces_per_cell": int(face_count.min()),
        "area_partition_error": float(abs(area.sum()-BOX**2)),
        "force_is_negative_energy_gradient_max_error": float(
            np.max(np.linalg.norm(force-fd_force, axis=1))),
        "volume_rate_finite_difference_max_error": float(
            np.max(np.abs(geom+weighted-fd_combined))),
        "weight_rate_finite_difference_max_error": float(
            np.max(np.abs(weighted-fd_weight))),
        "jacobian_rate_error": float(np.max(np.abs(weighted-jac@rate))),
        "jacobian_row_sum_max": float(np.max(np.abs(jac.sum(axis=1)))),
        "jacobian_column_sum_max": float(np.max(np.abs(jac.sum(axis=0)))),
        "global_force_norm": float(np.linalg.norm(force.sum(axis=0))),
        "kinetic_pressure_power": kinetic_power,
        "thermal_power_with_original_rate": thermal_power,
        "unbalanced_total_power": kinetic_power+thermal_power,
        "predicted_unbalanced_power": float(-pressure@weighted),
        "unbalanced_identity_error": float(abs(kinetic_power+thermal_power+pressure@weighted)),
        "original_weight_power": float(pressure@weighted),
        "projected_weight_power": float(pressure@neutral_weight_rate),
        "projected_total_power": kinetic_power+neutral_thermal_power,
        "projected_rate_relative_change_l2": float(
            np.linalg.norm(projected-rate)/np.linalg.norm(rate)),
        "scope": "Finite nonzero W and Wdot, 25-cell periodic mathematical power diagram. Force/volume/work adjoint and one global energy-neutral rate are checked; no local material-face law, native CUDA, RK4 finite step, or long Kepler evolution is proven.",
    }
    checks = {
        "area_tiles_periodic_domain": result["area_partition_error"] < 1e-10,
        "finite_weight_force_adjoint": result["force_is_negative_energy_gradient_max_error"] < 1e-7,
        "finite_weight_geometric_rate": result["volume_rate_finite_difference_max_error"] < 1e-7,
        "finite_weight_jacobian_rate": result["weight_rate_finite_difference_max_error"] < 1e-7,
        "nonzero_rate_energy_identity": result["unbalanced_identity_error"] < 1e-10,
        "energy_neutral_projection": abs(result["projected_total_power"]) < 1e-10,
        "nontrivial_weight_work": abs(result["original_weight_power"]) > 1e-4,
        "total_force_zero": result["global_force_norm"] < 1e-10,
    }
    result["checks"] = {key: bool(value) for key, value in checks.items()}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    if not all(checks.values()):
        raise SystemExit("Finite-weight adjoint audit failed")


if __name__ == "__main__":
    main()
