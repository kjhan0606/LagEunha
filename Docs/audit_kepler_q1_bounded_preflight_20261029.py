"""Independent Q1 native initial geometry/RHS and tiny accepted-step audit."""
import argparse
import csv
import json
from pathlib import Path
import re
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import energy
from audit_kepler_wall_native_band_probe_20261029 import polygon, reconstruct


def audit(root):
    cfg = json.loads((root/'preflight_stage.json').read_text())
    ic = json.loads((root/'ic_v1/initial_audit.json').read_text())
    a, ma = ordered_global(Path(cfg['initial']))
    ids = a['u4if'] & ID_MASK
    p, mp = ordered_global(root/'k1_preflight/khdiag.step000001.wall_rhs_probe.dat', ids)
    endpoints = [f for f in (root/'tiny_preflight').glob('khout.*.dat') if f.name != 'khout.000001.dat']
    if len(endpoints) != 1:
        raise ValueError('Missing unique tiny accepted endpoint')
    b, mb = ordered_global(endpoints[0], ids)
    independent = np.load(root/'ic_v1/independent_areas.npy')
    volume_error = float(np.max(abs(p['volume'].astype(float)/independent-1)))
    with (root/'k1_preflight/rk4_stage_r0.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if sorted(int(v['id']) for v in rows) != sorted(cfg['IDs']) or any(v['rhs_valid'] != '1' for v in rows):
        raise ValueError('Incomplete selected native K1 RHS data')
    xy = np.column_stack((p['x'], p['y']));tree = cKDTree(xy)
    selected = []
    for row in rows:
        global_id = int(row['id']);i = int(np.searchsorted(ids, global_id))
        poly, area, residual = polygon(i, xy, tree)
        vertices = np.array([v[:2] for v in poly])
        lengths = np.linalg.norm(np.roll(vertices,-1,axis=0)-vertices,axis=1)
        radius = float(np.max(np.linalg.norm(vertices,axis=1)))
        resolution = 64*np.finfo(float).eps*max(1.,radius)
        omitted_length = float(np.sum(lengths[lengths<=resolution]))
        # Degenerate co-circular edges can have roundoff-scale length and an
        # undefined numerical normal. Bound their total possible contribution
        # instead of testing that normal. No native face is removed or changed.
        nearest = float(tree.query(xy[i],k=2)[0][1])
        pmax = float(np.max(p['pressure']))
        vmax = float(np.max(np.hypot(p['vx'],p['vy'])))
        acceleration_bound = 2*pmax*omitted_length*(.5+2*radius/nearest)/float(p['mass'][i])
        work_bound = 4*pmax*vmax*omitted_length*(1+2*radius/nearest)
        acc, edot, udot, closure, _ = reconstruct(i, xy, p, poly,edge_resolution=resolution)
        native = np.array([float(row['ax_hydro']), float(row['ay_hydro'])])
        native_edot = float(row['die_stored'])
        native_udot = native_edot-float(p['mass'][i])*float(np.array([p['vx'][i],p['vy'][i]])@native)
        selected.append({'id': global_id, 'native_pressure_acceleration_norm': float(np.linalg.norm(native)),
                         'independent_acceleration_error': float(np.linalg.norm(acc-native)),
                         'independent_energy_work_error': abs(edot-native_edot),
                         'independent_thermal_work_error': abs(udot-native_udot),
                         'independent_volume_relative_error': abs(area/float(p['volume'][i])-1),
                         'roundoff_edge_resolution': resolution, 'omitted_roundoff_edge_length': omitted_length,
                         'omitted_edge_acceleration_bound': acceleration_bound, 'omitted_edge_work_bound': work_bound,
                         'all_site_support_residual': residual, 'face_closure': closure})
    text = (root/'tiny_preflight/log').read_text(errors='replace')
    ledgers = [line for line in text.splitlines() if '[RK4E]' in line]
    zero = bool(ledgers) and all(float(re.search(r'\b'+key+r'=(\S+)', line)[1]) == 0.
                               for line in ledgers for key in ['floor_cum','de_cum','sfl_cum'])
    drift = (energy(b, 'baseline')-energy(a, 'baseline'))/abs(energy(a, 'baseline'))
    g = cfg['predeclared_gates']
    area = float(np.sum(b['volume'].astype(float), dtype=np.longdouble))
    checks = {'independent_IC_gate': ic['passed'],
              'both_native_exits_zero': all((root/n/'exit_status.txt').read_text().strip() == '0' for n in ['k1_preflight','tiny_preflight']),
              'expected_population_and_IDs': len(p) == len(b) == 90352,
              'K1_no_time_or_coordinate_advance': mp['t'] == ma['t'] == 0. and all(np.array_equal(a[k],p[k]) for k in ['x','y','mass','ie']),
              'accepted_exact_target': mb['t'] == 1e-9,
              'fixed_mass_W0': np.array_equal(a['mass'],b['mass']) and np.all(b['w2'] == 0) and np.all(p['w2'] == 0),
              'all_owned_finite_positive': all(np.isfinite(q[k]).all() and np.all(q[k] > 0) for q in [p,b] for k in ['mass','ie','volume','den','pressure']),
              'all_native_K1_volumes_match_independent': volume_error < g['max_native_volume_relative_to_independent'],
              'bounded_domain_area': abs(area-ic['domain_area']) < 1e-5,
              'independent_selected_force_work': all(v['independent_acceleration_error']+v['omitted_edge_acceleration_bound'] < g['max_independent_acceleration_error'] and max(v['independent_energy_work_error'],v['independent_thermal_work_error'])+v['omitted_edge_work_bound'] < g['max_independent_work_error'] for v in selected),
              'initial_uniform_pressure_force_small': max(v['native_pressure_acceleration_norm'] for v in selected) < g['max_initial_pressure_acceleration'],
              'zero_injections': zero, 'tiny_step_energy': abs(drift) < g['restart_relative_energy_drift'],
              'no_unresolved_wall_or_state_abort': '[WALL_EVENT_STOP]' not in text and '[RK4_INVALID_STATE]' not in text}
    return {'passed': bool(all(checks.values())), 'checks': {k: bool(v) for k,v in checks.items()},
            'native_K1_volume_max_relative_to_independent': volume_error, 'selected_RHS': selected,
            'accepted_endpoint': str(endpoints[0]), 'accepted_time': mb['t'], 'relative_energy_change': drift,
            'endpoint_area': area, 'scope': 'New bounded Q1 native initial-state validation only. The exact t=3,6.5,10 and GIZMO long-time physics gates remain required.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path);parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args();result = audit(args.root)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2);stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit(1)
