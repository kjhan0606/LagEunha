"""Build a locally refined W=0 Kepler glass from the full-q2 t=0 snapshot.

The q2 disk is ring-sampled, not a Cartesian lattice. This generator keeps
every original site at r>=r_match bitwise unchanged and replaces the inner
sites with staggered polar rings. The target spacing factor is two at
r<=r_refine and transitions smoothly to one at r_match. Periodic areas and
the fixed-mass rebase are supplied by the existing sidecar/native workflow.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np

from plot_kepler_native_density_long_20261026 import complete_ordered


BOX = 6.0
CENTER = 3.0
BASE_COUNT = 384 * 384


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def spacing_factor(r, r_refine, r_match):
    if r <= r_refine:
        return 2.0
    if r >= r_match:
        return 1.0
    u = (r - r_refine) / (r_match - r_refine)
    smooth = 3.0*u*u - 2.0*u*u*u
    return 2.0 - smooth


def ring_sites(radii, dr_base):
    sites = [(CENTER, CENTER)]
    ring_counts = []
    for j, r in enumerate(radii, start=1):
        factor = spacing_factor(float(r), R_REFINE, R_MATCH)
        count = max(6, int(np.rint(4.0*np.pi*float(r)*factor/dr_base)))
        phase = 0.5 * (j % 2)
        theta = 2.0*np.pi*(np.arange(count, dtype=np.float64) + phase)/count
        sites.extend(zip(CENTER + r*np.cos(theta), CENTER + r*np.sin(theta)))
        ring_counts.append(count)
    return np.asarray(sites, dtype=np.float64), ring_counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_snapshot", type=Path)
    parser.add_argument("output_glass", type=Path)
    parser.add_argument("--r-refine", type=float, default=0.75)
    parser.add_argument("--r-match", type=float, default=1.0)
    parser.add_argument("--nominal-grid-side", type=int, default=448)
    args = parser.parse_args()
    global R_REFINE, R_MATCH
    R_REFINE, R_MATCH = args.r_refine, args.r_match
    if not (0 < R_REFINE < R_MATCH < 2.0):
        raise ValueError("Require 0 < r_refine < r_match < 2")
    manifest_path = args.output_glass.with_suffix(".json")
    if args.output_glass.exists() or manifest_path.exists():
        raise FileExistsError("Refusing to overwrite the refined glass or manifest")
    if args.nominal_grid_side <= 0:
        raise ValueError("Nominal grid side must be positive")

    p, meta = complete_ordered(args.base_snapshot)
    n0 = len(p)
    if n0 != BASE_COUNT or float(meta["t"]) != 0:
        raise ValueError("Expected the complete full-q2 384^2 t=0 snapshot")
    if (float(meta["Lx"]), float(meta["Ly"])) != (BOX, BOX):
        raise ValueError("Expected the periodic 6x6 q2 domain")
    if not np.all(p["w2"] == 0):
        raise ValueError("Base q2 snapshot must have exactly zero weights")
    x0 = p["x"].astype(np.float64)
    y0 = p["y"].astype(np.float64)
    if not np.all(np.isfinite(x0) & np.isfinite(y0) &
                  (x0 >= 0) & (x0 < BOX) & (y0 >= 0) & (y0 < BOX)):
        raise ValueError("Invalid base q2 positions")

    radius0 = np.hypot(x0 - CENTER, y0 - CENTER)
    original_radii = np.unique(np.round(radius0[(radius0 >= 0.1) &
                                                (radius0 < 1.5)], 10))
    dr_candidates = np.diff(original_radii)
    if len(dr_candidates) < 10:
        raise ValueError("Could not infer the q2 radial-ring spacing")
    dr_base = float(np.median(dr_candidates))
    if float(np.std(dr_candidates) / dr_base) > 2e-5:
        raise ValueError("Base q2 radii are not a uniform polar-ring sequence")

    # Equal increments of integral(f(r)/dr_base dr) give radial spacing
    # dr_base/f(r); f transitions smoothly from 2 to 1.
    radial_grid = np.linspace(0.0, R_MATCH, 20001)
    factors = np.array([spacing_factor(float(r), R_REFINE, R_MATCH)
                        for r in radial_grid])
    phase_coordinate = np.zeros_like(radial_grid)
    phase_coordinate[1:] = np.cumsum(
        0.5*(factors[1:] + factors[:-1]) * np.diff(radial_grid) / dr_base)
    ring_indices = np.arange(1, int(np.floor(phase_coordinate[-1])) + 1,
                             dtype=np.float64)
    new_radii = np.interp(ring_indices, phase_coordinate, radial_grid)
    new_radii = new_radii[new_radii < R_MATCH]
    inner_xy, ring_counts = ring_sites(new_radii, dr_base)

    keep = radius0 >= R_MATCH
    exterior_xy = np.column_stack((x0[keep], y0[keep]))
    xy = np.vstack((inner_xy, exterior_xy))
    if len(xy) >= args.nominal_grid_side**2:
        raise ValueError("Refined site count exceeds nominal NX*NY capacity")
    if not np.all((xy >= 0.0) & (xy < BOX)):
        raise ValueError("Refined positions fall outside the periodic box")

    # Verify the outer geometry is preserved as an exact coordinate multiset.
    expected_outer = np.column_stack((x0[keep], y0[keep]))
    actual_outer = xy[len(inner_xy):]
    if not np.array_equal(actual_outer, expected_outer):
        raise AssertionError("Exterior q2 coordinates changed")

    # Periodic nearest-neighbor distances provide an early geometry guard.
    from scipy.spatial import cKDTree
    nearest = cKDTree(xy, boxsize=BOX).query(xy, k=2, workers=1)[0][:, 1]
    if not np.all(np.isfinite(nearest)) or float(nearest.min()) < 0.004:
        raise ValueError(f"Refined geometry has a close pair: dmin={nearest.min():.9g}")

    output_glass = args.output_glass
    with output_glass.open("xb") as stream:
        stream.write(struct.pack("<idd", len(xy), BOX, BOX))
        stream.write(np.asarray(xy[:, 0], dtype="<f8").tobytes())
        stream.write(np.asarray(xy[:, 1], dtype="<f8").tobytes())

    r_new = np.hypot(xy[:, 0] - CENTER, xy[:, 1] - CENTER)
    bands = [(0.0, 0.5), (0.5, R_REFINE), (R_REFINE, 0.85),
             (0.85, 0.93), (0.93, R_MATCH), (R_MATCH, 2.0),
             (2.0, 2.75), (2.75, 5.0)]
    radial_counts = [{"lo": lo, "hi": hi,
                      "base": int(np.sum((radius0 >= lo) & (radius0 < hi))),
                      "refined": int(np.sum((r_new >= lo) & (r_new < hi)))}
                     for lo, hi in bands]
    report = {
        "base_snapshot": str(args.base_snapshot),
        "base_snapshot_sha256": sha256_file(args.base_snapshot),
        "base_time": float(meta["t"]),
        "base_count": n0,
        "output_glass": str(output_glass),
        "output_glass_sha256": sha256_file(output_glass),
        "output_count": int(len(xy)),
        "nominal_NX_NY": [args.nominal_grid_side, args.nominal_grid_side],
        "domain": [BOX, BOX],
        "r_refine": R_REFINE,
        "r_match": R_MATCH,
        "base_radial_spacing": dr_base,
        "refined_ring_count": int(len(new_radii)),
        "first_last_refined_ring_radius": [float(new_radii[0]), float(new_radii[-1])],
        "ring_count_min_max": [int(min(ring_counts)), int(max(ring_counts))],
        "min_periodic_site_separation": float(nearest.min()),
        "median_periodic_site_separation": float(np.median(nearest)),
        "exterior_site_count": int(np.sum(keep)),
        "exterior_coordinates_bitwise_preserved": True,
        "radial_counts": radial_counts,
        "initial_weights": "all zero in source q2 state; native glass stores positions only",
        "mass_and_area": "Use the validated periodic area sidecar and native fixed-mass rebase; this geometry-only generator assigns no masses.",
        "scope": "Locally refined W=0 Voronoi IC candidate only; no native stability claim",
    }
    with manifest_path.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
