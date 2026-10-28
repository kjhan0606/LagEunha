"""Velocity components and origin angular momentum from cached raw integrals.

No geometry rebuild, quadrature, fit, RHS, or evolution. Missing q32 coverage
at a new maximum is recorded explicitly for the bundled reference review.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference-dir', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--supplemental-q32', type=Path)
args = parser.parse_args()
root = args.reference_dir
report = json.loads((root/'report.json').read_text())
data = np.load(root/'reference.npz')
extra_sources = []
if args.supplemental_q32:
    supplement = np.load(args.supplemental_q32)
    assert np.array_equal(supplement['ids'], data['ids'])
    assert str(supplement['checkpoint_sha256']) == report['checkpoint_sha256']
    data = {key: data[key].copy() for key in data.files}
    indices = supplement['indices']
    assert not data['refined_mask'][indices].any()
    data['reference_q32_selected'][indices] = supplement['values']
    data['refined_mask'][indices] = True
    extra_sources.append(args.supplemental_q32)
checkpoint_path = Path(report['checkpoint'])
assert hashlib.sha256(checkpoint_path.read_bytes()).hexdigest() == report['checkpoint_sha256']
checkpoint = np.load(checkpoint_path)
ids, mass = data['ids'], data['mass']
assert np.array_equal(checkpoint['ids'], ids) and np.array_equal(checkpoint['mass'], mass)
a = checkpoint['state'][:-4].reshape(-1, 20)
assert np.array_equal(np.column_stack([mass, a[:, 3:]]), data['raw'])
assert float(checkpoint['time']) == float(data['time']) == report['time']
axes = a[:, :2]/np.linalg.norm(a[:, :2], axis=1)[:, None]
assert np.isfinite(axes).all()
azimuthal = np.column_stack([-axes[:, 1], axes[:, 0]])
r = data['initial_radius']
frozen = data['frozen']
regions = dict(active=~frozen, frozen=frozen, inner=r < .4,
    inner_transition=(r >= .4)&(r < .48), contact_band=(r >= .48)&(r < .52),
    outer_interior=(r >= .52)&(r < .55), active_edge=(r >= .55)&(r < .60),
    contact_active=data['contact_mask']&~frozen,
    boundary_active=data['boundary_mask']&~frozen)


def fields(raw, area):
    m = raw[:, 0]
    velocity = raw[:, 7:9]/m[:, None]
    thermal = raw[:, 17]-.5*(raw[:, 13]+raw[:, 16])
    result = dict(velocity_radial=np.sum(velocity*axes, axis=1),
        velocity_azimuthal=np.sum(velocity*azimuthal, axis=1),
        angular_momentum_origin=raw[:, 10]-raw[:, 11],
        specific_thermal=thermal/m, mean_density=m/area)
    assert all(np.isfinite(x).all() for x in result.values())
    return result


def statistics(delta, mask, coverage=False):
    ix = np.flatnonzero(mask)
    if len(ix) == 0:
        return dict(cells=0)
    k = ix[np.argmax(abs(delta[ix]))]
    result = dict(cells=len(ix), mass_weighted_RMS=float(np.sqrt(np.average(delta[ix]**2, weights=mass[ix]))),
        maximum_absolute=float(abs(delta[k])), maximum_signed=float(delta[k]), maximum_id=int(ids[k]))
    if coverage:
        result['maximum_cell_q32_refined'] = bool(data['refined_mask'][k])
    return result


numerical = fields(data['raw'], data['areas'])
reference16 = fields(data['reference_q16'][:, :18], data['reference_q16'][:, 19])
reference32 = fields(data['reference_q32_selected'][:, :18], data['reference_q32_selected'][:, 19])
errors = {}
refinement = {}
for name in numerical:
    errors[name] = {region: statistics(numerical[name]-reference16[name], mask, True)
                    for region, mask in regions.items()}
    refinement[name] = {region: statistics(reference32[name]-reference16[name], mask & data['refined_mask'])
                        for region, mask in regions.items()}
totals = {region: dict(numerical=float(numerical['angular_momentum_origin'][mask].sum()),
    exact_q16=float(reference16['angular_momentum_origin'][mask].sum()),
    exact_q32_selected_else_q16=float(reference32['angular_momentum_origin'][mask].sum()))
    for region, mask in regions.items()}
out = dict(scope=__doc__, time=report['time'], checkpoint=str(checkpoint_path),
    numerical_minus_q16=errors, q32_minus_q16_on_selected_cells=refinement,
    angular_momentum_origin_region_sums=totals,
    velocity_basis='Common to numerical/reference within each endpoint: radial from current generator position about Kepler origin; azimuthal=(-er_y,er_x). Not material-centroid orbit scoring.',
    source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in (Path(__file__), root/'reference.npz', root/'report.json', checkpoint_path, *extra_sources)},
    supplemental_q32=str(args.supplemental_q32) if args.supplemental_q32 else None,
    limits='q32 data exist only on original selected mask. New maxima outside it remain unrefined. Region sums compare different fields on same current cells, not a closed-system conservation test. No long-time or order verdict.')
with args.output.open('x') as stream:
    json.dump(out, stream, indent=2)
print(json.dumps(dict(status='completed', output=str(args.output)), indent=2))
