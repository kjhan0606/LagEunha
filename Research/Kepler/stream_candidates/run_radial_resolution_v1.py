"""Execute the approved IC comparison with the unchanged native CUDA binary."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

import numpy as np
from audit_snapshot import read_snapshot

EXPECTED_BINARY = 'f6f51995a76861b8438ee6649185b860ceb671c3209097005c7fdf6e13cb250f'
EDGES = [0, .047, .1, .2, .5, 1, 2, 4, 8, 18]
MASK = np.uint64((1 << 56)-1)


def ordered(path):
    p, meta = read_snapshot(path)
    return p[np.argsort(p['u4if'] & MASK)], meta


def audit_state(path, initial_path, gamma):
    p, meta = ordered(path)
    initial, _ = ordered(initial_path)
    if not np.array_equal(p['u4if'] & MASK, initial['u4if'] & MASK):
        raise ValueError('native IDs changed')
    mass_fixed = np.array_equal(p['mass'], initial['mass'])
    weights_zero = bool(np.all(p['w2'] == 0))
    x0 = np.column_stack((initial['x']-12, initial['y']-12))
    x = np.column_stack((p['x']-12, p['y']-12))
    v = np.column_stack((p['vx'], p['vy']))
    r0 = np.linalg.norm(x0, axis=1)
    omega = (r0*r0+.047**2)**(-.75)
    theta = omega*meta['t']
    cs, sn = np.cos(theta), np.sin(theta)
    exact = np.column_stack((cs*x0[:, 0]-sn*x0[:, 1], sn*x0[:, 0]+cs*x0[:, 1]))
    exact_v = omega[:, None]*np.column_stack((-exact[:, 1], exact[:, 0]))
    uad = initial['ie'].astype(float)*(initial['volume'].astype(float)/p['volume'])**(gamma-1)
    ratio = p['ie']/uad
    position_error = np.linalg.norm(x-exact, axis=1)/np.maximum(r0, 1e-30)
    velocity_error = np.linalg.norm(v-exact_v, axis=1)/np.maximum(r0*omega, 1e-30)
    phase_error = np.arctan2(exact[:, 0]*x[:, 1]-exact[:, 1]*x[:, 0], np.sum(exact*x, axis=1))
    row = dict(snapshot=str(path), time=meta['t'], n=len(p), mass_bitwise_fixed=mass_fixed,
               weights_identically_zero=weights_zero, minimum_U=float(p['ie'].min()),
               maximum_density_mass_volume_relative=float(np.max(abs(p['den']-p['mass']/p['volume'])/p['den'])),
               radial=[])
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        sel = (r0 >= lo) & (r0 < hi)
        if not sel.any():
            continue
        m = p['mass'][sel]
        avg = lambda a: float(np.sum(m*a[sel])/m.sum())
        angular = p['mass']*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])
        angular0 = initial['mass']*(x0[:, 0]*initial['vy']-x0[:, 1]*initial['vx'])
        row['radial'].append(dict(Rlo=lo, Rhi=hi, n=int(sel.sum()),
            mass=float(m.sum()), entropy_signed=avg(ratio-1), entropy_L1=avg(abs(ratio-1)),
            U_over_Uad_min=float(ratio[sel].min()), U_over_Uad_max=float(ratio[sel].max()),
            orbital_position_L1=avg(position_error), orbital_velocity_L1=avg(velocity_error),
            principal_phase_error_L1=avg(abs(phase_error)),
            radial_drift_L1=avg(abs(np.linalg.norm(x, axis=1)-r0)/np.maximum(r0, 1e-30)),
            angular_momentum_relative_change=float((angular[sel].sum()-angular0[sel].sum())/angular0[sel].sum())))
    valid = mass_fixed and weights_zero and np.all(np.isfinite(ratio)) and np.all(p['ie'] > 0)
    # Fixed engineering continuation gates, not GIZMO comparison thresholds.
    physical = [a for a in row['radial'] if a['Rlo'] < 8]
    valid = valid and all(a['entropy_L1'] <= .1 and a['U_over_Uad_min'] >= .1
                         and a['U_over_Uad_max'] <= 10 and a['orbital_position_L1'] <= .01
                         and a['orbital_velocity_L1'] <= .01 for a in physical)
    row['eligible_for_longer_continuation'] = bool(valid)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--case', choices=['uniform', 'graded'], required=True)
    parser.add_argument('--version', default='v2')
    parser.add_argument('--batch', default='native_v2')
    parser.add_argument('--reuse-preflight', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    base = root/('radial_resolution_20261021_'+args.version)
    config = json.loads((base/'manifest.json').read_text())
    if args.case not in config['cases']:
        raise ValueError('requested IC is not complete')
    own_ic = config['cases'][args.case]
    if own_ic['max_relative_mass_quadrature_change'] > 1e-8 or abs(own_ic['total_area']-576) > 1e-8:
        raise ValueError('requested IC area/quadrature gate failed')
    historical = json.loads((root/'core_refinement/manifest.json').read_text())['cases']['coarse_integrated']['mass']
    if abs(own_ic['total_mass']-historical)/historical > 1e-8:
        raise ValueError('total mass differs from the preserved original continuum IC')
    binary = root/'raw_vertex_dedup/v4/eunha2_raw_vertex_dedup.exe'
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    if digest != EXPECTED_BINARY:
        raise ValueError('frozen native binary SHA differs')
    dest = base/args.batch/args.case
    dest.mkdir(parents=True, exist_ok=False)
    initial_path = base/args.case/'khout.000001.dat'
    if hashlib.sha256(initial_path.read_bytes()).hexdigest() != own_ic['IC_sha256']:
        raise ValueError('requested IC hash differs')
    result = dict(scope='matched IC resolution only, unchanged native RK4/old limiter; no intrinsic-spin/remap claim',
                  case=args.case, binary_sha256=digest, started=time.time(), stages=[],
                  IC_sha256=own_ic['IC_sha256'], paired_IC_preparation_complete_at_launch=bool(config.get('IC_preparation_pass')),
                  requested_long_target=10, first_interval_target=.2,
                  physical_continuation_gates=dict(entropy_L1=.1, U_over_Uad_min=.1,
                      U_over_Uad_max=10, orbital_position_L1=.01, orbital_velocity_L1=.01),
                  output_policy='initial one-step geometry plus .05/.1/.2 checkpoints; failure start/stage retained; no per-step snapshots')
    def save():
        (dest/'manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    save()

    def execute(label, input_path, stop, ranks=1, faces=False, audit=False):
        d = dest/label
        d.mkdir()
        step = int(input_path.name.split('.')[1])
        (d/input_path.name).symlink_to(input_path)
        (d/'params.dat').symlink_to(base/args.case/'params.dat')
        (d/'eunha2.exe').symlink_to(binary)
        env = os.environ.copy()
        for key in list(env):
            if key.startswith(('GFS_', 'EUNHA_', 'HYDRO_', 'KEPLER_', 'SEDOV_', 'LAGMFM_')):
                env.pop(key)
        env.update(SEDOV_PHASE1='1', EUNHA_KEPLER_CX='12', EUNHA_KEPLER_CY='12',
                   EUNHA_KEPLER_RIN='2', EUNHA_KEPLER_ROUT='8', EUNHA_KEPLER_EPS='.047',
                   EUNHA_KEPLER_DTETA='.06', EUNHA_DBG_LO='1000000000', EUNHA_DBG_HI='0',
                   GFS_FLOOR_LOG='1', GFS_E_AUDIT='1', GFS_KEPLER_HALO_CELLS='8')
        env.update(OMP_NUM_THREADS='2' if ranks == 1 else '1', OMP_PROC_BIND='false',
                   NATIVE_MPI_RANKS=str(ranks), I_MPI_PIN='0', I_MPI_HYDRA_BOOTSTRAP='fork',
                   KEPLER_RESTART_STEP=str(step), EUNHA_IC='kepler', HYDRO_TSTOP=str(stop),
                   EUNHA_DUMP_DT=str(stop), GFS_EXACT_STOP='1', GFS_FAIL_NONPOSITIVE_U='1',
                   GFS_RK4_FAILURE_SNAPSHOT='1' if ranks == 1 else '0',
                   GFS_GPU_FINAL_VOLUME='1', GFS_GPU_RAW_FACES='1', GFS_VELOCITY_LIMITER='0',
                   GFS_GPU_FACE_DUMP_CALL='1' if faces else '0', GFS_GPU_TESS_DUMP='0')
        for key in ('I_MPI_PMI_LIBRARY', 'PMI_RANK', 'PMI_SIZE', 'PMI_FD', 'PMI_PORT',
                    'GFS_DIAG_DT_MAX', 'OMP_PLACES', 'GOMP_CPU_AFFINITY', 'KMP_AFFINITY',
                    'GFS_AFFINE_UNLIMITED_VELOCITY', 'GFS_AFFINE_UNMERGED_FACES',
                    'GFS_RK4_AUDIT_STEP_MIN', 'GFS_RK4_AUDIT_STEP_MAX', 'GFS_RK4_AUDIT_STAGE',
                    'GFS_RK4_AUDIT_IDS', 'GFS_SSPRK3', 'GFS_SSPRK_RETRY'):
            env.pop(key, None)
        if audit:
            env.update(GFS_RK4_AUDIT_STEP_MIN=str(step), GFS_RK4_AUDIT_STEP_MAX=str(step),
                       GFS_RK4_AUDIT_STAGE='k1_rhs')
        (d/'environment.json').write_text(json.dumps({k: v for k, v in env.items()
            if k.startswith(('GFS_', 'EUNHA_', 'HYDRO_', 'KEPLER_', 'OMP_', 'I_MPI_', 'NATIVE_'))}, indent=2)+'\n')
        start = time.monotonic()
        code = subprocess.run([sys.executable, str(Path(__file__).with_name('run_native.py'))],
                              cwd=d, env=env).returncode
        row = dict(label=label, requested_endpoint=stop, ranks=ranks, exit_code=code,
                   seconds=time.monotonic()-start)
        result['stages'].append(row)
        if code:
            row['classification'] = 'FAILED; no restart from invalid stage'
            starts = list(d.glob('khdiag.step*.start.dat'))
            if starts:
                row['normal_step_start'] = str(starts[-1])
                row['normal_start_physical_audit'] = audit_state(starts[-1], initial_path, config['gamma'])
            result['status'] = 'NATIVE_FAILURE'
            save()
            return None
        candidates = [p for p in d.glob('khout.*.dat') if not p.is_symlink()]
        path = max(candidates, key=lambda p: int(p.name.split('.')[1]))
        _, meta = read_snapshot(path)
        if meta['t'] != stop:
            raise ValueError('native exact endpoint differs')
        row['checkpoint'] = str(path)
        text = (d/'log').read_text()
        if 'halo_search_clipped=1' in text:
            raise ValueError('halo search clipped')
        floors = [float(x) for x in re.findall(r'floor_(?:inj|cum)=\s*([^\s]+)', text)]
        if any(x != 0 for x in floors):
            raise ValueError('native floor injection')
        row['floor_injection_zero'] = True
        save()
        return path

    if args.reuse_preflight:
        d = args.reuse_preflight.resolve()
        old_timing = json.loads((d/'timing.json').read_text())
        if (old_timing['exit_code'] != 0 or (d/'khout.000001.dat').resolve() != initial_path.resolve()
                or (d/'eunha2.exe').resolve() != binary.resolve()):
            raise ValueError('reused preflight provenance failed')
        first = d/'khout.000002.dat'
        if read_snapshot(first)[1]['t'] != 1e-7:
            raise ValueError('reused preflight time differs')
        result['stages'].append(dict(label='preflight_mpi1', ranks=1, exit_code=0,
            reused_native_output=str(d), checkpoint=str(first),
            note='reclassify face audit; do not repeat native step'))
    else:
        first = execute('preflight_mpi1', initial_path, 1e-7, faces=True, audit=True)
    if first is None:
        return 2
    d = first.parent
    states = np.genfromtxt(d/'rk4_stage_r0.tsv', names=True, dtype=None, encoding=None, delimiter='\t')
    states = states[np.argsort(states['id'])]
    expected = np.load(base/args.case/'geometry.npz')['volume']
    error = abs(states['V']-expected)/expected
    if error.max() >= 2e-6 or not np.all(states['U_raw'] > 0):
        raise ValueError('native first volume/positive-state gate failed')
    subprocess.run([sys.executable, str(Path(__file__).with_name('audit_face_reciprocity_v2.py')),
                    '--call', '1', '--output', str(dest/'legacy_reciprocity.json'), str(d)], check=True)
    reciprocity = next(iter(json.loads((dest/'legacy_reciprocity.json').read_text()).values()))
    subprocess.run([sys.executable, str(Path(__file__).with_name('audit_radial_faces_v1.py')),
                    str(d), str(dest/'explicit_face_audit.json')], check=True)
    explicit = json.loads((dest/'explicit_face_audit.json').read_text())
    if not explicit['pass']:
        raise ValueError('native shared-face geometry gate failed')
    first_rhs = []
    radius = np.load(base/args.case/'geometry.npz')['radius']
    faces = np.loadtxt(d/'face_dump_r0_call1.tsv', skiprows=1, usecols=(4, 5, 7, 10, 11))
    initial, _ = ordered(initial_path)
    ii, jj = faces[:, 0].astype(int), faces[:, 1].astype(int)
    length = np.hypot(faces[:, 3], faces[:, 4])
    usable = length > 0
    generator_jump = np.zeros(len(faces))
    generator_jump[usable] = ((initial['vx'][jj[usable]]-initial['vx'][ii[usable]])*faces[usable, 3]
        +(initial['vy'][jj[usable]]-initial['vy'][ii[usable]])*faces[usable, 4])/length[usable]
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        sel = (radius >= lo) & (radius < hi)
        if sel.any():
            fs = sel[ii]
            first_rhs.append(dict(Rlo=lo, Rhi=hi, n=int(sel.sum()),
                signed_Udot=float(states['Udot'][sel].sum()),
                positive_Udot=float(np.maximum(states['Udot'][sel], 0).sum()),
                Udot_over_U_min=float(np.min(states['Udot'][sel]/states['U_raw'][sel])),
                Udot_over_U_max=float(np.max(states['Udot'][sel]/states['U_raw'][sel])),
                max_pstar_over_initial_pressure=float(np.max(faces[fs, 2])/config['pressure']),
                face_ends_pstar_over_initial_pressure_gt_100=int(np.sum(faces[fs, 2] > 100*config['pressure'])),
                generator_normal_jump_rms=float(np.sqrt(np.mean(generator_jump[fs]**2)))))
    result['preflight'] = dict(native_volume_relative_max=float(error.max()),
                               reciprocity=reciprocity, explicit_tiny_face_audit=explicit, first_RHS=first_rhs,
                               reconstructed_face_velocity_jump='not present in frozen v4 face packet; generator jump is separately labeled')
    save()
    if args.case == 'graded':
        parallel = execute('preflight_mpi4', initial_path, 1e-7, ranks=4)
        if parallel is None:
            return 2
        a, ma = ordered(first)
        b, mb = ordered(parallel)
        errors = {k: float(np.max(abs(a[k].astype(float)-b[k].astype(float))))
                  for k in a.dtype.names if k != 'u4if'}
        tolerance_failures = []
        normalized = {}
        for k in errors:
            typ = a.dtype.fields[k][0]
            multiplier = 8 if typ.itemsize == 4 else 256
            scale = np.maximum(abs(a[k].astype(float)), abs(b[k].astype(float)))
            floor = 24 if k in ('x', 'y') else (1 if k in ('vx', 'vy') or k.startswith('stress.') else 0)
            tolerance = multiplier*np.finfo(typ).eps*np.maximum(scale, floor)
            difference = abs(a[k].astype(float)-b[k].astype(float))
            normalized[k] = float(np.max(np.divide(difference, tolerance,
                out=np.zeros_like(difference), where=tolerance > 0)))
            if np.any(difference > tolerance):
                tolerance_failures.append(k)
        result['mpi1_mpi4'] = dict(all_recorded_fields_bitwise_equal=not any(errors.values()),
                                   max_absolute_errors=errors, metadata_equal=ma == mb,
                                   normalized_roundoff_comparison=normalized,
                                   comparison_failures=tolerance_failures,
                                   criterion='8eps float fields;256eps double fields; scale floors x/y24,v/gradient1; fixed before observing comparison')
        save()
        if tolerance_failures or ma != mb or not np.array_equal(a['u4if'] & MASK, b['u4if'] & MASK):
            raise ValueError('new-layout MPI1/4 agreement failed; inspect before evolution')
    current = first
    for stop in (.05, .1, .2):
        current = execute(f't{stop:g}', current, stop)
        if current is None:
            return 2
        audit = audit_state(current, initial_path, config['gamma'])
        result['stages'][-1]['physical_audit'] = audit
        # Complete the initial .2 diagnosis while positive, even when a
        # local accuracy gate already rules out promotion to expensive t10.
        save()
    result['status'] = ('SHORT_INTERVAL_COMPLETE_REQUIRES_HALF_STEP_CONTROL'
        if result['stages'][-1]['physical_audit']['eligible_for_longer_continuation']
        else 'SHORT_INTERVAL_PHYSICAL_ACCURACY_FAILED')
    result['long_goal_achieved'] = False
    save()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
