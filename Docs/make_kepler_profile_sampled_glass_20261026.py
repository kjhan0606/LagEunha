"""Experimental profile-sampled Kepler glass with the q=2 disk kept intact.

Allocate polar sites inside r<0.5 and 2<r<2.5 according to the Hopkins
surface-density profile at the measured median q=2 disk-cell mass. Outside
the polar region use a coarse periodic Cartesian background. This supplies a
geometry-only candidate for the native sparse-glass path, not an accepted IC.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree

DISK_MASS = 0.00024276668751554097
Q = 2.0
BOX = 6.0
CENTER = 3.0


def rho(r):
    return .01 + np.where(r < .5, (r/.5)**3,
                           np.where(r <= 2., 1., (1.+(r-2.)/.1)**-3))


def variable_rings(lo, hi, phase_offset):
    grid = np.linspace(lo, hi, 20001)
    rate = np.sqrt(rho(grid)/(Q*DISK_MASS))
    coordinate = np.empty_like(grid)
    coordinate[0] = 0.
    coordinate[1:] = np.cumsum(.5*(rate[1:]+rate[:-1])*np.diff(grid))
    ring_coordinates = np.arange(.5, coordinate[-1], 1.)
    radius = np.interp(ring_coordinates, coordinate, grid)
    points = []
    counts = []
    for j, r in enumerate(radius):
        ds = np.sqrt(DISK_MASS/(Q*float(rho(np.array(r)))))
        count = max(6, int(np.rint(2*np.pi*r/ds)))
        theta = 2*np.pi*(np.arange(count)+.5*((j+phase_offset)%2))/count
        points.append(np.column_stack((CENTER+r*np.cos(theta),
                                       CENTER+r*np.sin(theta))))
        counts.append(count)
    return np.vstack(points), radius, np.array(counts, dtype=int)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_q2_glass", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--inner-buffer-start", type=float, default=.5,
                        help="Keep original q=2 core sites from this radius to 0.5")
    parser.add_argument("--outer-buffer-end", type=float, default=2.,
                        help="Keep original q=2 outer sites from 2 to this radius")
    args = parser.parse_args()
    if not (0 < args.inner_buffer_start <= .5 and
            2. <= args.outer_buffer_end < 2.5):
        raise ValueError("Invalid original-q=2 transition buffers")
    raw = args.source_q2_glass.read_bytes()
    n, lx, ly = struct.unpack_from("<idd", raw)
    if (n, lx, ly) != (384*384, BOX, BOX) or len(raw) != 20+16*n:
        raise ValueError("Expected complete q=2 384² glass")
    x = np.frombuffer(raw, dtype="<f8", count=n, offset=20)
    y = np.frombuffer(raw, dtype="<f8", count=n, offset=20+8*n)
    r = np.hypot(x-CENTER, y-CENTER)
    disk = np.column_stack((x[(r>=.5)&(r<2.)], y[(r>=.5)&(r<2.)]))
    inner_buffer = np.column_stack((x[(r>=args.inner_buffer_start)&(r<.5)],
                                    y[(r>=args.inner_buffer_start)&(r<.5)]))
    outer_buffer = np.column_stack((x[(r>=2.)&(r<args.outer_buffer_end)],
                                    y[(r>=2.)&(r<args.outer_buffer_end)]))
    core, core_r, core_n = variable_rings(0., args.inner_buffer_start, 0)
    edge, edge_r, edge_n = variable_rings(args.outer_buffer_end, 2.5, len(core_r))
    ncart = int(np.rint(BOX/np.sqrt(DISK_MASS/.01)))
    hcart = BOX/ncart
    axis = (np.arange(ncart)+.5)*hcart
    gx, gy = np.meshgrid(axis, axis)
    background = np.column_stack((gx.ravel(), gy.ravel()))
    bg_r = np.hypot(background[:, 0]-CENTER, background[:, 1]-CENTER)
    background = background[bg_r >= 2.55]
    xy = np.vstack((core, inner_buffer, disk, outer_buffer, edge, background))
    if len(xy) >= n:
        raise ValueError("Profile sampling did not reduce the total site count")
    nearest = cKDTree(xy).query(xy, k=2)[0][:, 1]
    if not np.all(np.isfinite(nearest)) or float(np.min(nearest)) < .004:
        raise ValueError("Generator collision")
    if not ((xy > 0) & (xy < BOX)).all():
        raise ValueError("Out-of-box site")
    with args.output.open("xb") as stream:
        stream.write(struct.pack("<idd", len(xy), BOX, BOX))
        stream.write(np.asarray(xy[:, 0], dtype="<f8").tobytes())
        stream.write(np.asarray(xy[:, 1], dtype="<f8").tobytes())
    result = {
        "source": str(args.source_q2_glass),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "output": str(args.output),
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "target_disk_cell_mass": DISK_MASS,
        "site_counts": {"core": len(core), "inner_q2_buffer": len(inner_buffer),
                        "disk": len(disk), "outer_q2_buffer": len(outer_buffer),
                        "outer_polar": len(edge),
                        "outer_cartesian": len(background), "total": len(xy)},
        "buffer_radii": [args.inner_buffer_start, args.outer_buffer_end],
        "core_rings": len(core_r),
        "core_first_last_radius": [float(core_r[0]), float(core_r[-1])],
        "core_first_last_count": [int(core_n[0]), int(core_n[-1])],
        "outer_rings": len(edge_r),
        "outer_first_last_radius": [float(edge_r[0]), float(edge_r[-1])],
        "outer_first_last_count": [int(edge_n[0]), int(edge_n[-1])],
        "background_cartesian_side": ncart,
        "background_cartesian_spacing": hcart,
        "minimum_radius": float(np.min(np.hypot(xy[:, 0]-CENTER, xy[:, 1]-CENTER))),
        "minimum_pair_distance": float(np.min(nearest)),
        "disk_generators_bitwise_preserved": True,
        "scope": "Geometry candidate only; Voronoi areas, mass distribution, and native stability must be audited",
    }
    with args.output.with_suffix(".json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
