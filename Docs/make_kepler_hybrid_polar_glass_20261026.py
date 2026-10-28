"""Build an opt-in polar inner disk for the user's native glass-file IC path.

The square periodic box stays Cartesian outside radius 2.35. Particle count,
box, and Hopkins profile resolution stay unchanged. The native executable must
set EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1 to assign exact initial cell masses.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree

N = 384
L = 6.
CX = CY = 3.
R_INTERFACE = 2.35


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--aspect', type=float, choices=(1., 2.), required=True)
    parser.add_argument('--no-center', action='store_true',
                        help='Start rings at a/2 and let their cells meet at the origin')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    h = L / N
    radial_spacing = h * np.sqrt(args.aspect)
    target_azimuthal_spacing = h / np.sqrt(args.aspect)
    ids = np.arange(N * N, dtype=np.int64)
    x = (ids % N + .5) * h
    y = (ids // N + .5) * h
    r0 = np.hypot(x - CX, y - CY)
    selected = np.flatnonzero(r0 < R_INTERFACE)
    theta0 = np.mod(np.arctan2(y[selected] - CY, x[selected] - CX), 2 * np.pi)
    selected = selected[np.lexsort((theta0, r0[selected]))]
    if args.no_center:
        radius = (np.arange(int(np.floor((R_INTERFACE - .5*radial_spacing) /
                                         radial_spacing)) + 1) + .5) * radial_spacing
    else:
        radius = np.arange(1, int(np.floor((R_INTERFACE - .5 * radial_spacing) /
                                           radial_spacing)) + 1) * radial_spacing
    raw_count = 2 * np.pi * radius / target_azimuthal_spacing
    desired_count = len(selected) - (0 if args.no_center else 1)
    scaled = raw_count * desired_count / np.sum(raw_count)
    count = np.floor(scaled).astype(int)
    remainder = desired_count - int(np.sum(count))
    if remainder:
        order = np.argsort(-(scaled - count))
        count[order[:remainder]] += 1
    if np.any(count < 3) or (0 if args.no_center else 1) + int(np.sum(count)) != len(selected):
        raise ValueError('Polar ring count allocation failed')
    points = [] if args.no_center else [np.array([[CX, CY]])]
    for j, (r, nphi) in enumerate(zip(radius, count)):
        angle = 2 * np.pi * (np.arange(nphi) + .5 * (j % 2)) / nphi
        points.append(np.column_stack((CX + r * np.cos(angle),
                                       CY + r * np.sin(angle))))
    polar = np.vstack(points)
    x[selected] = polar[:, 0]
    y[selected] = polar[:, 1]
    xy = np.column_stack((x, y))
    nearest = cKDTree(xy).query(xy, k=2)[0][:, 1]
    r = np.hypot(x - CX, y - CY)
    ring = (r >= .5) & (r <= 2.)
    if np.min(nearest) < .35 * h or not ((xy >= 0) & (xy < L)).all():
        raise ValueError('Generator collision or out-of-box position')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as stream:
        stream.write(struct.pack('<idd', N * N, L, L))
        stream.write(x.astype('<f8').tobytes())
        stream.write(y.astype('<f8').tobytes())
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    result = dict(file=str(args.output), sha256=digest, n=N*N, box=L,
                  h=h, polar_interface_radius=R_INTERFACE,
                  polar_generators=len(selected), polar_rings=len(radius),
                  central_generator=not args.no_center,
                  minimum_radius=float(radius[0] if args.no_center else 0.),
                  radial_over_azimuthal_spacing=args.aspect,
                  target_radial_spacing=float(radial_spacing),
                  target_azimuthal_spacing=float(target_azimuthal_spacing),
                  first_last_ring_radii=[float(radius[0]), float(radius[-1])],
                  first_last_ring_counts=[int(count[0]), int(count[-1])],
                  minimum_any_pair_distance=float(np.min(nearest)),
                  ring_0p5_to_2_cells=int(np.sum(ring)),
                  ring_min_pair_distance=float(np.min(nearest[ring])),
                  requires_native_env='EUNHA_KEPLER_MASS_FROM_CELL_VOLUME=1')
    manifest = args.output.with_suffix('.json')
    with manifest.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
