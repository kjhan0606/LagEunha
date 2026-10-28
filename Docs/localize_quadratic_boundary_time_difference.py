"""Localize two-step-size differences on identical initial cells at two boundaries."""
import hashlib
import json
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
roots = {
    'r055': [base / f'Codex_continuous_quadratic_evolution_{tag}_20261018'
             for tag in ('dt1e4', 'dt5e5')],
    'r060': [base / f'Codex_quadratic_r060_{tag}_20261018'
             for tag in ('dt1e4', 'dt5e5')],
}
reference = None
results = {}
for label, folders in roots.items():
    states, sources = [], []
    for folder, step in zip(folders, (1, 2)):
        summary = json.loads((folder / 'summary.json').read_text())
        assert summary['status'] == 'completed'
        with np.load(folder / 'initial.npz') as initial:
            identity = {key: initial[key].copy() for key in ('state', 'mass', 'ids', 'factors')}
        if reference is None:
            reference = identity
        assert all(np.array_equal(reference[key], identity[key]) for key in reference)
        path = folder / f'step{step:04d}.npz'
        with np.load(path) as saved:
            assert abs(float(saved['time']) - 1e-4) < 1e-15
            assert np.array_equal(saved['ids'], reference['ids'])
            assert np.array_equal(saved['mass'], reference['mass'])
            states.append(saved['state'][:-4].reshape(-1, 20).copy())
        sources.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    mass, ids = reference['mass'], reference['ids']
    radius = np.linalg.norm(reference['state'][:-4].reshape(-1, 20)[:, :2], axis=1)
    dv = (states[1][:, 9:11] - states[0][:, 9:11]) / mass[:, None]
    square = np.sum(dv * dv, axis=1)
    assert np.isfinite(square).all()
    active = radius < .55
    total = float(np.sum(mass[active] * square[active]))
    rows = []
    for lo, hi in ((0, .4), (.4, .48), (.48, .52), (.52, .55)):
        mask = (radius >= lo) & (radius < hi)
        contribution = float(np.sum(mass[mask] * square[mask]))
        rows.append(dict(radius=[lo, hi], cells=int(mask.sum()),
                         velocity_mass_RMS=float(np.sqrt(np.average(square[mask], weights=mass[mask]))),
                         fraction_of_common_squared_difference=contribution / total if total else None))
    worst = np.flatnonzero(active)[np.argsort(square[active])[-12:][::-1]]
    results[label] = dict(sources=sources, regions=rows,
        largest_velocity_differences=[dict(id=int(ids[i]), radius=float(radius[i]),
            velocity_difference=float(np.sqrt(square[i]))) for i in worst])
print(json.dumps(dict(scope=__doc__, time=1e-4, results=results,
    limits='Two temporal resolutions only. Differences are not errors against an exact solution. '
           'Boundary sensitivity is not a causal proof or a convergence-order measurement.'), indent=2))
