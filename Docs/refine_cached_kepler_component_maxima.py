"""q32 only for extra component maxima absent from saved q32 selection."""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce')
sys.path[:0] = [str(base/'Codex_steady_ale_reference_modules_20261018'), str(base)]
from steady_kepler_ale_reference import field_moments, polygon_integral
from audit_canonical_pair_history import gravity_coefficient
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference-dir', type=Path, required=True)
parser.add_argument('--output-prefix', type=Path, required=True)
args = parser.parse_args()
root = args.reference_dir
report = json.loads((root/'report.json').read_text())
components = json.loads((root/'components.json').read_text())
data = np.load(root/'reference.npz')
assert components['time'] == report['time'] == float(data['time'])
for path, digest in report['source_sha256'].items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
ids_to_refine = sorted({row['maximum_id'] for regions in components['numerical_minus_q16'].values()
                       for row in regions.values() if row['cells'] and not row['maximum_cell_q32_refined']})
indices = np.array([int(np.flatnonzero(data['ids'] == i)[0]) for i in ids_to_refine], dtype=int)
assert not data['refined_mask'][indices].any()
path = base/'Codex_weight_remap_moments_centered_20261017/audit.py'
spec = importlib.util.spec_from_file_location('saved_potential', path)
potential = importlib.util.module_from_spec(spec)
spec.loader.exec_module(potential)
physical = report['physical']


def integrand(points):
    fields, _ = field_moments(points, gravity_coefficient, potential.potential,
        physical['pressure'], physical['gamma'], physical['rin'], physical['rout'], physical['floor'])
    return np.column_stack([fields, np.ones(len(points))])


values, estimates = [], []
for index in indices:
    lo, hi = data['polygon_offsets'][index:index+2]
    poly = data['polygon_vertices'][lo:hi]
    value, error = polygon_integral(poly, integrand, physical['split_radii'], 32)
    assert np.isfinite(value).all() and value[0] > 0 and value[19] > 0
    values.append(value)
    estimates.append(float(error))
npz = Path(str(args.output_prefix)+'.npz')
meta_path = Path(str(args.output_prefix)+'.json')
assert not npz.exists() and not meta_path.exists()
np.savez(npz, ids=data['ids'], indices=indices, values=np.array(values).reshape(-1, 20),
         checkpoint_sha256=report['checkpoint_sha256'])
meta = dict(scope=__doc__, indices=indices.tolist(), ids=ids_to_refine,
    maximum_quad_vec_estimate=max(estimates, default=0.), time=report['time'],
    source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                   (Path(__file__), root/'reference.npz', root/'report.json', root/'components.json', npz)},
    limits='Only original q16 maximum-error cells of the added derived components. No geometry/fit/RHS or evolution; base reference arrays preserved.')
meta_path.write_text(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
