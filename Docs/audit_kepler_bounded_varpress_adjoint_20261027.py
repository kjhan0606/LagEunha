"""Check the W=0 variational pressure law on a fixed convex boundary.

This independent small-mesh calculation clips every Voronoi cell to a
regular polygon. It verifies the interior-pair force and geometric area
rate against finite differences and measures the boundary-face term missing
from the current periodic-only CUDA total-energy work law. It does not
implement a native boundary or demonstrate long-time Kepler stability.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


GRID = 4
DOMAIN_SIDES = 64
CENTER = np.array([2.0, 2.0])
RADIUS = 2.9
GAMMA = 5.0 / 3.0
H = 1.0e-4
SEEDS = (17, 23)
SOURCE = (Path(__file__).with_name("KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27")
          / "Exam" / "exam_gpu.cu")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def boundary():
    angle = np.arange(DOMAIN_SIDES) * (2.0 * np.pi / DOMAIN_SIDES)
    return [(float(CENTER[0] + RADIUS * np.cos(a)),
             float(CENTER[1] + RADIUS * np.sin(a)), -k - 1)
            for k, a in enumerate(angle)]


def append_vertex(vertices, value):
    if vertices and vertices[-1][:2] == value[:2]:
        vertices[-1] = value
    else:
        vertices.append(value)


def clip(poly, normal, threshold, neighbor):
    out = []
    for e, a in enumerate(poly):
        b = poly[(e + 1) % len(poly)]
        da = a[0] * normal[0] + a[1] * normal[1] - threshold
        db = b[0] * normal[0] + b[1] * normal[1] - threshold
        inside_a, inside_b = da <= 0.0, db <= 0.0
        if inside_a:
            append_vertex(out, a)
        if inside_a != inside_b:
            t = da / (da - db)
            x = a[0] + t * (b[0] - a[0])
            y = a[1] + t * (b[1] - a[1])
            support = neighbor if inside_a else a[2]
            append_vertex(out, (x, y, support))
    if len(out) > 1 and out[0][:2] == out[-1][:2]:
        out.pop()
    if len(out) < 3:
        raise ValueError("Voronoi cell vanished during clipping")
    changed = True
    while changed:
        changed = False
        for e in range(len(out)):
            prev = (e - 1) % len(out)
            if out[e][2] >= 0 and out[e][2] == out[prev][2]:
                del out[e]
                changed = True
                break
    return out


def area(poly):
    xy = np.asarray([p[:2] for p in poly])
    return float(0.5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                               - xy[:, 1] * np.roll(xy[:, 0], -1)))


def tessellate(sites):
    domain = boundary()
    polygons = []
    areas = []
    for i, xi in enumerate(sites):
        poly = domain
        for j, xj in enumerate(sites):
            if i == j:
                continue
            normal = xj - xi
            threshold = 0.5 * (np.dot(xj, xj) - np.dot(xi, xi))
            poly = clip(poly, normal, threshold, j)
        value = area(poly)
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Invalid bounded Voronoi area")
        polygons.append(poly)
        areas.append(value)
    return polygons, np.asarray(areas)


def energy(areas, entropy_constant):
    return float(np.sum(entropy_constant * areas ** (1.0 - GAMMA),
                        dtype=np.longdouble) / (GAMMA - 1.0))


def face_operator(sites, velocity, pressure, polygons):
    force = np.zeros_like(sites)
    volume_rate = np.zeros(len(sites))
    interior_relative = np.zeros(len(sites))
    boundary_relative = np.zeros(len(sites))
    wall_traction = np.zeros(2)
    closure = np.zeros_like(sites)
    boundary_count = 0
    reciprocal = {}
    for i, poly in enumerate(polygons):
        for e, a in enumerate(poly):
            b = poly[(e + 1) % len(poly)]
            segment = np.asarray(b[:2]) - np.asarray(a[:2])
            length = float(np.linalg.norm(segment))
            if length <= 0:
                continue
            outward = np.array([segment[1], -segment[0]]) / length
            closure[i] += length * outward
            j = a[2]
            if j < 0:
                # The fixed domain face has zero laboratory-frame speed.
                boundary_relative[i] -= length * np.dot(velocity[i], outward)
                wall_traction += pressure[i] * length * outward
                boundary_count += 1
                continue
            delta = sites[j] - sites[i]
            distance = float(np.linalg.norm(delta))
            normal = delta / distance
            if np.dot(outward, normal) < 1.0 - 1.0e-9:
                raise ValueError("Interior face support has wrong orientation")
            centroid = 0.5 * (np.asarray(a[:2]) + np.asarray(b[:2]))
            offset = centroid - 0.5 * (sites[i] + sites[j])
            bi = 0.5 * length * normal + length * offset / distance
            force[i] += (pressure[i] - pressure[j]) * bi
            speed = (0.5 * np.dot(velocity[i] + velocity[j], normal)
                     - np.dot(velocity[j] - velocity[i], offset) / distance)
            volume_rate[i] += length * speed
            interior_relative[i] += length * (speed - np.dot(velocity[i], normal))
            reciprocal[(i, j)] = (length, centroid)
    pair_error = 0.0
    for (i, j), (length, centroid) in reciprocal.items():
        if (j, i) not in reciprocal:
            raise ValueError(f"Nonreciprocal bounded face {i},{j}")
        reverse_length, reverse_centroid = reciprocal[(j, i)]
        pair_error = max(pair_error, abs(length - reverse_length),
                         float(np.linalg.norm(centroid - reverse_centroid)))
    return (force, volume_rate, interior_relative, boundary_relative,
            wall_traction, float(np.max(np.linalg.norm(closure, axis=1))),
            boundary_count, pair_error)


def run_case(seed, uniform_pressure):
    rng = np.random.default_rng(seed)
    sites = np.array([(x + 0.5, y + 0.5)
                      for y in range(GRID) for x in range(GRID)], dtype=float)
    sites += rng.uniform(-0.12, 0.12, size=sites.shape)
    velocity = rng.normal(scale=0.2, size=sites.shape)
    pressure = (np.ones(len(sites)) if uniform_pressure else
                1.0 + 0.25 * rng.normal(size=len(sites)))
    polygons, areas = tessellate(sites)
    entropy_constant = pressure * areas ** GAMMA
    (force, volume_rate, interior_relative, boundary_relative, wall,
     closure, boundary_count, pair_error) = face_operator(
         sites, velocity, pressure, polygons)
    domain_area = area(boundary())
    plus = tessellate(sites + H * velocity)[1]
    minus = tessellate(sites - H * velocity)[1]
    finite_difference_rate = (plus - minus) / (2.0 * H)
    finite_difference_force = np.empty_like(sites)
    for i in range(len(sites)):
        for axis in range(2):
            delta = np.zeros_like(sites)
            delta[i, axis] = H
            ep = energy(tessellate(sites + delta)[1], entropy_constant)
            em = energy(tessellate(sites - delta)[1], entropy_constant)
            finite_difference_force[i, axis] = -(ep - em) / (2.0 * H)
    # Current periodic CUDA path omits boundary faces. In a bounded domain,
    # each omitted static edge must contribute -P*(-v_i.n*A) to total E_dot.
    power_without_boundary = float(np.sum(force * velocity)
                                   - np.dot(pressure, interior_relative))
    power_with_boundary = float(np.sum(force * velocity)
                                - np.dot(pressure,
                                         interior_relative + boundary_relative))
    return {
        "seed": seed, "uniform_pressure": uniform_pressure,
        "sites": len(sites), "boundary_faces": boundary_count,
        "area_partition_error": float(abs(np.sum(areas) - domain_area)),
        "face_vector_closure_max": closure,
        "reciprocal_face_error_max": pair_error,
        "force_gradient_error_max": float(np.linalg.norm(
            force - finite_difference_force, axis=1).max()),
        "area_rate_error_max": float(np.abs(
            volume_rate - finite_difference_rate).max()),
        "relative_sweep_error_max": float(np.abs(
            volume_rate - interior_relative - boundary_relative).max()),
        "total_area_rate_error": float(abs(np.sum(volume_rate))),
        "wall_momentum_balance_error": float(np.linalg.norm(
            force.sum(axis=0) + wall)),
        "energy_power_without_boundary": power_without_boundary,
        "energy_power_with_boundary": power_with_boundary,
        "force_max": float(np.linalg.norm(force, axis=1).max()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = [run_case(seed, uniform) for seed in SEEDS
             for uniform in (False, True)]
    maximum = lambda key: max(abs(case[key]) for case in cases)
    checks = {
        "bounded_area_partition": maximum("area_partition_error") < 1e-10,
        "reciprocal_interior_faces": maximum("reciprocal_face_error_max") < 1e-10,
        "closed_polygons": maximum("face_vector_closure_max") < 1e-10,
        "pressure_force_is_energy_gradient": maximum("force_gradient_error_max") < 1e-7,
        "geometric_area_rate": maximum("area_rate_error_max") < 1e-7,
        "boundary_relative_sweep_closes": maximum("relative_sweep_error_max") < 1e-10,
        "fixed_domain_area_rate_zero": maximum("total_area_rate_error") < 1e-10,
        "wall_traction_accounts_for_momentum": maximum("wall_momentum_balance_error") < 1e-10,
        "total_energy_power_closes_with_boundary": maximum("energy_power_with_boundary") < 1e-10,
        "unmodified_periodic_work_detectably_fails":
            max(abs(case["energy_power_without_boundary"]) for case in cases) > 1e-3,
        "uniform_pressure_force_zero": all(
            case["force_max"] < 1e-10 for case in cases if case["uniform_pressure"]),
    }
    result = {
        "source": str(SOURCE), "source_sha256": digest(SOURCE),
        "reproducer": str(Path(__file__).resolve()),
        "reproducer_sha256": digest(Path(__file__)),
        "grid": [GRID, GRID], "boundary_sides": DOMAIN_SIDES,
        "boundary_center": CENTER.tolist(), "boundary_circumradius": RADIUS,
        "gamma": GAMMA, "central_difference_h": H,
        "cases": cases, "checks": checks,
        "scope": "Small fixed convex-domain mathematical W0 audit; no native CUDA boundary, fluid reflection, multi-rank scaling, or long-time Kepler result",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"checks": checks, "cases": cases}, indent=2))
    if not all(checks.values()):
        raise SystemExit("Bounded W0 variational operator audit failed")


if __name__ == "__main__":
    main()
