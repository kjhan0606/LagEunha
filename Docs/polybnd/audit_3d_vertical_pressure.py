"""Read-only vertical pressure-gradient audit of an extruded vortex.

The derivative is affine, so its absolute maximum over a convex cell is
attained at a vertex. This audits the reconstructed polynomial, not an
enclosure of the continuum solution or a complete force budget.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from moment_state_3d import decode
from pressure_reconstruction_3d import prepare
from dense_record_3d import load


def gradient_z(metadata, points):
    z = (np.asarray(points) - metadata['center']) / metadata['scale']
    c = metadata['coefficient']
    # Features: x,y,z,xx,xy,xz,yy,yz,zz.
    return (c[2] + c[5]*z[:, 0] + c[7]*z[:, 1] + 2*c[8]*z[:, 2]) / metadata['scale']


def inspect(path, library):
    path = Path(path)
    original = path.read_bytes()
    if path.suffix == '.npz':
        record, saved = load(path)
        t = record['t'] + record['dt']
    else:
        saved = json.loads(original)
        t = saved['t']
    sites, weights, raw = (np.asarray(saved[k]) for k in ('sites', 'weights', 'raw'))
    initial = json.loads((path.parent/'initial_state.json').read_text())
    np.testing.assert_array_equal(raw[:, 0], initial['mass'])
    polys, _, _ = Geometry3D(library).build(sites, weights)
    states = [decode(poly, row) for poly, row in zip(polys, raw)]
    pressure, diagnostic = prepare(states, polys)
    peaks = np.array([np.max(abs(gradient_z(p, v))) for p, v in zip(pressure, polys)])
    worst = int(np.argmax(peaks))
    assert path.read_bytes() == original
    return dict(path=str(path), sha256=hashlib.sha256(original).hexdigest(), t=t,
                cells=len(raw), max_vertical_pressure_gradient=float(peaks[worst]),
                worst_cell=worst, worst_site=sites[worst].tolist(),
                max_mean_vertical_velocity=float(np.max(abs(raw[:, 15]/raw[:, 0]))),
                pressure_diagnostic=diagnostic)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('library')
    parser.add_argument('paths', nargs='+')
    args = parser.parse_args()
    report = dict(rows=[inspect(p, args.library) for p in args.paths],
                  full_validation=False, scope='Reconstructed vertical pressure gradient on preserved states.')
    Path('report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
