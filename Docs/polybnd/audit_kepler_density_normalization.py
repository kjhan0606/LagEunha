"""Independent continuum moment audit of a saved Kepler density closure.

No time evolution, mass renormalization, or modification of saved states.
The native fit residual is distinguished from independently integrated moments.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from audit_kepler_energy_step import configure, atomic
from kepler_disk_reference import polar_rule, frozen


def run(seed):
    root = Path(seed)
    data = (root / 'source.json').read_bytes()
    cfg = (root / 'configuration.json').read_bytes()
    identity = dict(source_sha256=hashlib.sha256(data).hexdigest(),
                    configuration_sha256=hashlib.sha256(cfg).hexdigest())
    configure(json.loads(cfg))
    cp = json.loads(data)
    mass = np.asarray(cp['mass'])
    a = np.asarray(cp['state'])[:-4].reshape(-1, 20)
    cells, _ = base.geometry(a[:, :2], a[:, 2], cp['t'])
    mask = frozen(cp['t'], a[:, :2])
    saved = dict(identity=identity, next_cell=0, rows=[])
    if Path('progress.json').exists():
        saved = json.loads(Path('progress.json').read_text())
        if saved['identity'] != identity:
            raise ValueError('input identity changed')
    for i in range(saved['next_cell'], len(cells)):
        if not mask[i]:
            raw = np.r_[mass[i], a[i, 3:]]
            d = base.decode(cells[i], raw)['density']
            points, weights = polar_rule(cells[i], 16, radii=(.25, .5))
            dm = weights * d.value(points)
            native = base.features(d.points, d.origin, d.scale)
            independent = base.features(points, d.origin, d.scale)
            error = dm @ independent / mass[i] - d.weighted_mass @ native / mass[i]
            saved['rows'].append(dict(cell=i,
                mass_relative=float(dm.sum()/mass[i]-1),
                centered_moment_difference=float(np.max(abs(error))),
                native_mass_relative=float(d.weighted_mass.sum()/mass[i]-1),
                native_fit_residual=d.residual))
        saved['next_cell'] = i+1
        if (i+1) % 32 == 0 or i+1 == len(cells):
            atomic('progress.json', saved)
            print('density audit', i+1, '/', len(cells), flush=True)
    worst = max(saved['rows'], key=lambda r: abs(r['mass_relative']))
    i = worst['cell']; raw = np.r_[mass[i], a[i, 3:]]
    d = base.decode(cells[i], raw)['density']
    refinement = []
    for order in (16, 24, 32):
        x, w = polar_rule(cells[i], order, radii=(.25, .5))
        refinement.append(dict(order=order, mass_relative=float(w @ d.value(x)/mass[i]-1)))
    atomic('worst_cell.json', dict(cell=i, polygon=cells[i].tolist(), raw=raw.tolist(),
        theta=d.theta.tolist(), logZ=d.logZ, refinement=refinement))
    report = dict(identity=identity, t=cp['t'], worst_mass=worst,
        worst_centered_moment=max(saved['rows'], key=lambda r:r['centered_moment_difference']),
        refinement=refinement, full_validation=False)
    atomic('report.json', report)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seed')
    run(parser.parse_args().seed)
