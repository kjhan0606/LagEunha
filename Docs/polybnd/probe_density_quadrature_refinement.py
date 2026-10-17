"""Refit one saved cell at independent quadrature orders; no evolution."""
import argparse, json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from kepler_disk_reference import polar_rule


def triangle_rule(poly, order):
    z, w = np.polynomial.legendre.leggauss(order)
    z = (z+1)/2; w = w/2
    center = poly.mean(0); points = []; weights = []
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        u = a-center; v = b-center
        determinant = abs(u[0]*v[1]-u[1]*v[0])
        for r, wr in zip(z, w):
            points.extend(center+r*u+((1-r)*z)[:, None]*v)
            weights.extend(determinant*(1-r)*wr*w)
    return np.asarray(points), np.asarray(weights)


def run(fixture):
    data = json.loads(Path(fixture).read_text())
    poly = np.asarray(data['polygon']); raw = np.asarray(data['raw'])
    points, weights = polar_rule(poly, 32, radii=(.25, .5))
    old = base.quadrature; rows = []
    try:
        for order in (6, 8, 12, 16, 24):
            base.quadrature = lambda p: triangle_rule(p, order)
            density = base.decode(poly, raw)['density']
            dm = weights*density.value(points)
            native = base.features(density.points, density.origin, density.scale)
            reference = base.features(points, density.origin, density.scale)
            error = dm @ reference/raw[0] - density.weighted_mass @ native/raw[0]
            rows.append(dict(order=order, native_fit_residual=density.residual,
                independent_mass_relative=float(dm.sum()/raw[0]-1),
                independent_centered_moment_difference=float(np.max(abs(error)))))
    finally:
        base.quadrature = old
    report = dict(cell=data['cell'], rows=rows, scope='one saved cell refit, no trajectory validation')
    Path('refinement_report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('fixture')
    run(p.parse_args().fixture)
