"""Compare generator-centered and cell-averaged Kepler gravity in sparse cores.

Only sites near the origin are tessellated. Two neighborhood radii are checked
to verify that omitted distant generators do not affect the core cells.
This is an IC/source-quadrature diagnostic, not a native hydrodynamic run.
"""

import argparse
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import Voronoi, cKDTree


CENTER = np.array([3.0, 3.0])
EPS = 0.01175


def glass(path):
    with path.open("rb") as stream:
        n, lx, ly = struct.unpack("<idd", stream.read(20))
        x = np.fromfile(stream, dtype="<f8", count=n)
        y = np.fromfile(stream, dtype="<f8", count=n)
    if (not 0 < n <= 384*384 or (lx, ly) != (6., 6.) or
            len(x) != n or len(y) != n or path.stat().st_size != 20+16*n):
        raise ValueError(f"Bad glass: {path}")
    return np.column_stack((x, y))


def polygon_for_sites(xy, wanted, neighborhood_radius):
    select = np.flatnonzero(np.linalg.norm(xy-CENTER, axis=1) < neighborhood_radius)
    local = xy[select]
    lookup = {int(global_id): int(local_id) for local_id, global_id in enumerate(select)}
    vor = Voronoi(local, qhull_options="Qbb Qc Qz")
    polygons = []
    for global_id in wanted:
        local_id = lookup[int(global_id)]
        region = vor.regions[vor.point_region[local_id]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded core cell {global_id} in radius {neighborhood_radius}")
        polygons.append(vor.vertices[region])
    return polygons, len(select)


def area(poly):
    return .5*abs(np.dot(poly[:, 0], np.roll(poly[:, 1], -1))-
                  np.dot(poly[:, 1], np.roll(poly[:, 0], -1)))


def integrate_polygon(poly, order):
    roots, weights = np.polynomial.legendre.leggauss(order)
    u = .5*(roots+1.)
    w = .5*weights
    uu, vv = np.meshgrid(u, u, indexing="ij")
    ww = np.outer(w, w)
    values = np.zeros(5, dtype=np.float64)
    a = poly[0]
    for k in range(1, len(poly)-1):
        b, c = poly[k], poly[k+1]
        det = abs(np.linalg.det(np.column_stack((b-a, c-a))))
        xy = a + uu[..., None]*(b-a) + ((1.-uu)*vv)[..., None]*(c-a)
        offset = xy-CENTER
        r2 = np.sum(offset*offset, axis=-1)
        accel = -offset/(r2+EPS*EPS)[..., None]**1.5
        rho = .01+(np.sqrt(r2)/.5)**3
        jac_w = ww*det*(1.-uu)
        values[0] += np.sum(jac_w)
        values[1:3] += np.sum(jac_w[..., None]*accel, axis=(0, 1))
        values[3] += np.sum(jac_w*rho)
        values[4] += np.sum(jac_w*rho*np.sum(accel*offset, axis=-1))
    return values


def audit(name, path, positions, match_tolerance):
    xy = glass(path)
    distances, wanted = cKDTree(xy).query(positions)
    if np.max(distances) > match_tolerance or len(set(wanted.tolist())) != len(positions):
        raise ValueError(f"Expected core-ring angular counterparts are absent from {path}")
    short, nshort = polygon_for_sites(xy, wanted, .35)
    long, nlong = polygon_for_sites(xy, wanted, .5)
    area_short = np.array([area(p) for p in short])
    area_long = np.array([area(p) for p in long])
    if np.max(abs(area_short/area_long-1.)) > 1e-9:
        raise ValueError(f"Voronoi core polygons depend on neighborhood: {path}")
    entries = []
    for site, poly in zip(xy[wanted], long):
        q64 = integrate_polygon(poly, 64)
        q128 = integrate_polygon(poly, 128)
        if abs(q128[0]/area(poly)-1.) > 1e-10:
            raise ValueError("Quadrature does not recover polygon area")
        offset = site-CENTER
        r = np.linalg.norm(offset)
        outward = offset/r
        a_site = -offset/(r*r+EPS*EPS)**1.5
        a_mean = q128[1:3]/q128[0]
        a64 = q64[1:3]/q64[0]
        rho_site = .01+(r/.5)**3
        entries.append({
            "site_xy": site.tolist(), "site_radius": float(r),
            "area": float(q128[0]),
            "area_max_vertex_radius": float(np.max(np.linalg.norm(poly-CENTER, axis=1))),
            "generator_inward_acceleration": float(-np.dot(a_site, outward)),
            "cell_mean_inward_acceleration": float(-np.dot(a_mean, outward)),
            "cell_mean_over_generator_radial_acceleration":
                float(np.dot(a_mean, outward)/np.dot(a_site, outward)),
            "cell_mean_tangential_acceleration":
                float(outward[0]*a_mean[1]-outward[1]*a_mean[0]),
            "cell_mean_density_target_over_site_density_target":
                float((q128[3]/q128[0])/rho_site),
            "quadrature_64_vs_128_acceleration_max_abs": float(np.max(abs(a64-a_mean))),
        })
    ratio = np.array([e["cell_mean_over_generator_radial_acceleration"] for e in entries])
    return {
        "name": name, "glass": str(path), "cells": len(xy),
        "max_distance_to_profile_sampled_core_site": float(np.max(distances)),
        "local_sites_radius_0p35": nshort, "local_sites_radius_0p5": nlong,
        "max_polygon_area_change_between_local_radii": float(np.max(abs(area_short-area_long))),
        "radial_acceleration_ratio_min_median_max":
            [float(np.min(ratio)), float(np.median(ratio)), float(np.max(ratio))],
        "core_cells": entries,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("narrow", type=Path)
    parser.add_argument("wide", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    narrow_xy = glass(args.narrow)
    rr = np.linalg.norm(narrow_xy-CENTER, axis=1)
    rmin = np.min(rr)
    core_sites = narrow_xy[np.abs(rr-rmin)<1e-9]
    if len(core_sites) != 8 or not (.1 < rmin < .102):
        raise ValueError("Expected eight symmetric profile-sampled core-ring sites")
    result = {
        "softening": EPS,
        "core_ring_sites": len(core_sites),
        "core_ring_radius": float(rmin),
        "cases": [
            audit("full_q2_nearby", args.reference, core_sites, .02),
            audit("narrow", args.narrow, core_sites, 1e-12),
            audit("wide", args.wide, core_sites, 1e-8),
        ],
        "scope": "Initial 2D Voronoi cell average of prescribed Plummer acceleration under piecewise-constant cell density. Native source uses generator-centered acceleration; this diagnostic does not prove later hydrodynamic instability or require per-step subcell integration.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"core_ring_sites":result["core_ring_sites"],
                      "core_ring_radius":result["core_ring_radius"],
                      "cases":[{"name":c["name"], "cells":c["cells"],
                                "ratio":c["radial_acceleration_ratio_min_median_max"],
                                "area_error":c["max_polygon_area_change_between_local_radii"]}
                               for c in result["cases"]]},indent=2))


if __name__ == "__main__":
    main()
