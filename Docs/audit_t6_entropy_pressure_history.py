"""Read-only pressure decomposition using original and inherited entropy."""
import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha

base = Path('/gpfs/kjhan/LagForce')
initial = base / 'Codex_canonical_initial_20261017/run/khout.000000.dat'
manifest = json.loads((base / 'Codex_adaptive_qr_physical_to_t10_20261018/manifest.json').read_text())
source = manifest['cases']['adaptive_from_t2']
seed = Path(source['source'])
assert sha(seed) == source['source_sha256']
paths = {
    'adaptive': base / 'Codex_adaptive_qr_physical_to_t10_20261018/adaptive_from_t2/khout.001783.dat',
    'legacy': base / 'Codex_smooth_qr_physical_to_t10_20261018/legacy_continue/khout.001783.dat',
}
ic, t0 = read(initial)
s, ts = read(seed)
assert t0 == 0 and np.array_equal(ic['id'], s['id'])
gamma = float(np.float32(1.666667))
r0 = np.linalg.norm(np.column_stack((ic['x']-3, ic['y']-3)), axis=1)
masks = {'inner_background': (r0 < .5) & (ic['den'] < .5),
         'inner_shear_band': (r0 >= .3) & (r0 < .5),
         'disk_core': (r0 >= .6) & (r0 < 1.9) & (ic['den'] > .5)}
result = {'initial': str(initial), 'initial_sha256': sha(initial),
          'seed': str(seed), 'seed_sha256': sha(seed), 'seed_time': ts,
          'scope': 'Algebraic pressure-error decomposition on each actual geometry. Signed projections sum to one but are not causal fractions or force/work budgets. No state or entropy is modified.',
          'cases': {}}
for name, path in paths.items():
    a, t = read(path)
    assert np.array_equal(a['id'], ic['id'])
    for z in (ic, s, a):
        assert np.array_equal(z['mass'], ic['mass']) and np.all(z['w2'] == 0)
        assert np.all(np.isfinite(z['pressure'])) and np.all(z['pressure'] > 0)
        assert np.all(np.isfinite(z['volume'])) and np.all(z['volume'] > 0)
    p0iso = ic['pressure'] * (ic['volume']/a['volume'])**gamma
    p2iso = s['pressure'] * (s['volume']/a['volume'])**gamma
    parts = {'initial_entropy_compression': p0iso/ic['pressure'] - 1,
             'inherited_entropy_before_t2': (p2iso-p0iso)/ic['pressure'],
             'entropy_after_t2': (a['pressure']-p2iso)/ic['pressure']}
    error = a['pressure']/ic['pressure'] - 1
    residual = float(np.max(np.abs(sum(parts.values())-error)))
    assert residual < 1e-11
    row = {'snapshot': str(path), 'sha256': sha(path), 'time': t,
           'decomposition_max_abs_residual': residual, 'regions': {}}
    for region, mask in masks.items():
        def mean(x):
            return float(np.average(x[mask], weights=ic['mass'][mask]))
        variance = mean(error**2)
        stats = {key: {'RMS': mean(val**2)**.5,
                       'signed_projection_on_total': mean(val*error)/variance}
                 for key, val in parts.items()}
        closure = sum(v['signed_projection_on_total'] for v in stats.values())
        assert abs(closure-1) < 1e-12
        row['regions'][region] = {'cells': int(mask.sum()),
                                  'total_pressure_error_RMS': variance**.5,
                                  'components': stats}
    result['cases'][name] = row
assert result['cases']['adaptive']['time'] == result['cases']['legacy']['time']
output = base / 'Codex_adaptive_qr_physical_to_t10_20261018/t6_entropy_pressure_history.json'
if output.exists():
    assert json.loads(output.read_text()) == result
else:
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
print(json.dumps(result, indent=2))
