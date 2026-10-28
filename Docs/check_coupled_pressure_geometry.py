"""Check independent Voronoi volumes for the saved pressure-response audit."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi

sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha

parser = argparse.ArgumentParser()
parser.add_argument('--step', type=int, required=True)
args = parser.parse_args()
base = Path('/gpfs/kjhan/LagForce')
root = base / 'Codex_coupled_onset_to_t4_20261018'
seed, _ = read(base / 'Codex_coupled_onset_seed_20261018/reference/khout.000184.dat')
out = {}
for name in ['plus', 'minus']:
    path = root / name / f'khout.{args.step:06d}.dat'
    a, time = read(path)
    assert np.array_equal(a['id'], seed['id'])
    assert np.array_equal(a['mass'], seed['mass'])
    assert np.all(a['w2'] == 0)
    x = np.column_stack((a['x'], a['y']))
    mesh = Voronoi(x)
    selected = np.flatnonzero(np.hypot(x[:, 0]-3, x[:, 1]-3) < 2.4)
    errors = []
    for i in selected:
        region = mesh.regions[mesh.point_region[i]]
        assert len(region) >= 3 and -1 not in region
        q = mesh.vertices[region] - x[i]
        volume = abs(np.sum(q[:, 0]*np.roll(q[:, 1], -1)
                            - q[:, 1]*np.roll(q[:, 0], -1))) / 2
        errors.append(abs(volume/a['volume'][i]-1))
    error = float(max(errors))
    assert error < 1e-10
    out[name] = {'time': time, 'checked_cells': len(selected),
                 'max_relative_volume_error': error, 'snapshot_sha256': sha(path)}
assert abs(out['plus']['time']-out['minus']['time']) < 1e-12
out['scope'] = 'Interior unweighted geometry and fixed masses; not a force or full energy-budget check.'
with (root / f'pressure_geometry_check_{args.step:06d}.json').open('x') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
