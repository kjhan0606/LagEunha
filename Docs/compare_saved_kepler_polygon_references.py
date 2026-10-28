"""Equal-time contact endpoints and each endpoint's own Eulerian reference.

Saved-state analysis only. Moving-polygon raw differences are path differences,
not pointwise errors on a common volume and not a time-order measurement.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--coarse-root', type=Path, required=True)
parser.add_argument('--fine-root', type=Path, required=True)
parser.add_argument('--coarse-reference', type=Path, required=True)
parser.add_argument('--fine-reference', type=Path, required=True)
parser.add_argument('--components-name', default='components.json')
args = parser.parse_args()
roots = (args.coarse_root, args.fine_root)
refs = (args.coarse_reference, args.fine_reference)
manifests = [json.loads((root/'manifest.json').read_text()) for root in roots]
summaries = [json.loads((root/'summary.json').read_text()) for root in roots]
assert all(s['status'] == 'completed' for s in summaries)
for key in ('sha256', 'input', 'end', 'active_radius', 'density_closure',
            'shared_face_velocity', 'integrator'):
    assert manifests[0][key] == manifests[1][key], key
assert manifests[0]['dt'] == 2*manifests[1]['dt']
initial = [np.load(root/'initial.npz') for root in roots]
for key in ('state', 'factors', 'mass', 'ids', 'time', 'weight_rate_seed'):
    assert np.array_equal(initial[0][key], initial[1][key]), key
reports = [json.loads((ref/'report.json').read_text()) for ref in refs]
arrays = [np.load(ref/'reference.npz') for ref in refs]
assert reports[0]['physical'] == reports[1]['physical']
assert reports[0]['time'] == reports[1]['time'] == .001
reference_hashes = []
checkpoints = []
admissibility = []
for root, rootref, manifest, summary, report, data in zip(roots, refs, manifests, summaries, reports, arrays):
    path = Path(report['checkpoint'])
    steps = round((manifest['end']-manifest.get('start_time', 0.))/manifest['dt'])
    assert path.resolve() == (root/f'step{steps:04d}.npz').resolve()
    assert Path(report['initial']).resolve() == (root/'initial.npz').resolve()
    assert hashlib.sha256((root/'initial.npz').read_bytes()).hexdigest() == report['initial_sha256']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == report['checkpoint_sha256']
    checkpoint = np.load(path)
    checkpoints.append(checkpoint)
    assert float(checkpoint['time']) == float(data['time']) == summary['completed_time'] == manifest['end'] == report['time']
    raw = np.column_stack([checkpoint['mass'], checkpoint['state'][:-4].reshape(-1, 20)[:, 3:]])
    assert np.array_equal(raw, data['raw'])
    for key in ('mass', 'ids', 'factors'):
        assert np.array_equal(checkpoint[key], data[key]), key
    state_cells = checkpoint['state'][:-4].reshape(-1, 20)
    initial_cells = initial[0]['state'][:-4].reshape(-1, 20)
    frozen_mask = data['frozen']
    assert np.array_equal(checkpoint['mass'], initial[0]['mass'])
    assert np.array_equal(state_cells[frozen_mask], initial_cells[frozen_mask])
    assert np.array_equal(checkpoint['factors'][frozen_mask], initial[0]['factors'][frozen_mask])
    factor_R = checkpoint['factors']@checkpoint['factors'].transpose(0, 2, 1)
    thermal = raw[:, 17]-.5*(raw[:, 13]+raw[:, 16])
    assert np.isfinite(state_cells).all() and np.isfinite(factor_R).all()
    minimum_thermal = float(thermal.min())
    minimum_R = float(np.linalg.eigvalsh(factor_R).min())
    assert minimum_thermal > 0 and minimum_R > 0
    admissibility.append(dict(mass_bitwise_unchanged=True, frozen_state_bitwise_unchanged=True,
        frozen_factors_bitwise_unchanged=True, minimum_thermal=minimum_thermal,
        minimum_factor_R_eigenvalue=minimum_R,
        limits='State admissibility and frozen storage preservation; moving frozen polygons are not thereby an exact physical reservoir.'))
    # Area is independently recalculated from the polygons retained by the
    # reference job; checkpoint alone does not contain a polygon array.
    vertices = data['polygon_vertices']
    offsets = data['polygon_offsets']
    areas = []
    for lo, hi in zip(offsets[:-1], offsets[1:]):
        p = vertices[lo:hi]-vertices[lo]
        areas.append(.5*np.sum(p[:, 0]*np.roll(p[:, 1], -1)-p[:, 1]*np.roll(p[:, 0], -1)))
    assert np.array_equal(np.array(areas), data['areas'])
    reference_hashes.append(hashlib.sha256((rootref/'reference.npz').read_bytes()).hexdigest())
assert reports[0]['source_sha256'] == reports[1]['source_sha256']
for path, expected in reports[0]['source_sha256'].items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, path
for key in ('ids', 'mass', 'frozen', 'initial_radius'):
    assert np.array_equal(arrays[0][key], arrays[1][key]), key
ids = arrays[0]['ids']
mass = arrays[0]['mass']
r = arrays[0]['initial_radius']
frozen = arrays[0]['frozen']
regions = dict(active=~frozen, frozen=frozen, inner=r < .4,
               inner_transition=(r >= .4)&(r < .48), contact_band=(r >= .48)&(r < .52),
               outer_interior=(r >= .52)&(r < .55), active_edge=(r >= .55)&(r < .60))
fields = []
gamma = reports[0]['physical']['gamma']
P0 = reports[0]['physical']['pressure']
initial_sites = initial[0]['state'][:-4].reshape(-1, 20)[:, :2]
radial = initial_sites/np.linalg.norm(initial_sites, axis=1)[:, None]
azimuthal = np.column_stack([-radial[:, 1], radial[:, 0]])
for data, checkpoint in zip(arrays, checkpoints):
    raw = data['raw']
    state = checkpoint['state'][:-4].reshape(-1, 20)
    thermal = raw[:, 17]-.5*(raw[:, 13]+raw[:, 16])
    pressure = (gamma-1)*thermal/data['areas']
    fields.append(dict(raw18_mixed_units_auxiliary=raw, mass=raw[:, 0], q=raw[:, 1:3],
                       Q=raw[:, 3:7], momentum=raw[:, 7:9], D=raw[:, 9:13], Pi=raw[:, 13:17],
                       energy=raw[:, 17], centroid=raw[:, 1:3]/mass[:, None],
                       velocity=raw[:, 7:9]/mass[:, None],
                       velocity_radial_initial_basis=np.sum(raw[:, 7:9]/mass[:, None]*radial, axis=1),
                       velocity_azimuthal_initial_basis=np.sum(raw[:, 7:9]/mass[:, None]*azimuthal, axis=1),
                       ell=raw[:, 10]-raw[:, 11]-(raw[:, 1]*raw[:, 8]-raw[:, 2]*raw[:, 7])/mass,
                       angular_momentum_origin=raw[:, 10]-raw[:, 11],
                       thermal=thermal, pressure=pressure, pressure_relative_P0=pressure/P0,
                       generator=state[:, :2], relative_weight=state[:, 2]-state[0, 2],
                       factor_R=data['factors']@data['factors'].transpose(0, 2, 1),
                       areas=data['areas']))
differences = {}
for name in fields[0]:
    delta = fields[1][name]-fields[0][name]
    norm = np.linalg.norm(delta.reshape(len(ids), -1), axis=1)
    differences[name] = {}
    for region, mask in regions.items():
        ix = np.flatnonzero(mask)
        if not len(ix):
            differences[name][region] = dict(cells=0)
            continue
        k = ix[np.argmax(norm[ix])]
        differences[name][region] = dict(cells=len(ix), maximum=float(norm[k]), maximum_id=int(ids[k]),
            mass_weighted_RMS=float(np.sqrt(np.average(norm[ix]**2, weights=mass[ix]))))
component_reports = []
component_hashes = []
for ref, report in zip(refs, reports):
    path = ref/args.components_name
    components = json.loads(path.read_text())
    assert components['time'] == report['time']
    assert Path(components['checkpoint']).resolve() == Path(report['checkpoint']).resolve()
    for source, expected in components['source_sha256'].items():
        assert hashlib.sha256(Path(source).read_bytes()).hexdigest() == expected
    component_reports.append(components)
    component_hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
print(json.dumps(dict(scope=__doc__, time=reports[0]['time'], same_initial_arrays=True,
    run_manifests=[str(root/'manifest.json') for root in roots],
    reference_reports=[str(ref/'report.json') for ref in refs],
    reference_npz_sha256=reference_hashes,
    endpoint_fine_minus_coarse=differences,
    physical_errors=reports,
    state_admissibility=admissibility,
    physical_error_components=component_reports, component_report_sha256=component_hashes,
    common_velocity_basis='Endpoint path difference projections use common initial generator radial/azimuthal axes. Each endpoint exact-reference component report uses its own current generator axes.',
    weight_gauge=dict(convention='W_i minus W of first common ID', reference_id=int(ids[0])),
    limits='No two-level order claim, long-term success or candidate ranking. Each endpoint reference uses its own moving polygon; frozen reservoir is an open boundary. Contact/boundary masks are endpoint-specific and only individually reported, not used as common comparison masks. raw18_mixed_units_auxiliary is not a physical error norm. Reference refinement assessment belongs to Astra.'), indent=2))
