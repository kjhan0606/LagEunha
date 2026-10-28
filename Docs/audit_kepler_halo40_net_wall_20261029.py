"""Gate one dedicated-Ada halo40 canonical/combined dumps0 ABBA benchmark."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys

Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
PARENT = Q1/'halo40_reuse_alias_t2654918_ada_20261029'
ORIGINAL = Q1/'vertex_halo_pair_t2654918_ada_20261029'
STAGE = Q1/'halo40_net_wall_ada_t2654918_20261029'
CASES = [('old_a', 'qualifier_old'), ('combined_b', 'combined40'),
         ('combined_b2', 'combined40'), ('old_a2', 'qualifier_old')]
GPU = 'GPU-42b3367f-1df7-2b2d-7961-7781c87b022d'
BASE_SHA = '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'
COMBINED_SHA = 'a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def performance_summary(wall):
    if set(wall) != {name for name, _ in CASES} or not all(
            math.isfinite(value) and value > 0 for value in wall.values()):
        raise ValueError('Expected four positive finite solver wall measurements')
    old = [wall['old_a'], wall['old_a2']]
    new = [wall['combined_b'], wall['combined_b2']]
    conservative = min(old)/max(new)
    return {'solver_wall_seconds': wall,
            'canonical_median_seconds': statistics.median(old),
            'combined_median_seconds': statistics.median(new),
            'median_speedup': statistics.median(old)/statistics.median(new),
            'conservative_observed_speedup': conservative,
            'canonical_repeat_max_over_min': max(old)/min(old),
            'combined_repeat_max_over_min': max(new)/min(new),
            'clear_improvement': conservative > 1.05,
            'predeclared_clarity_rule': 'fastest canonical / slowest combined > 1.05',
            'scope': 'This GPU and interval only; no A100 throughput or long-time stability claim'}


def stage():
    result = json.loads((PARENT/'comparison.json').read_text())
    if not result['passed'] or (PARENT/'manual_exit_status.txt').read_text().strip() != '0':
        raise ValueError('Independent exact halo40 combination must finish and pass first')
    subprocess.run(['sha256sum', '-c', 'audit_inputs.sha256'], cwd=PARENT,
                   stdout=subprocess.DEVNULL, check=True)
    for source, expected in [(ORIGINAL/'qualifier_old', BASE_SHA), (PARENT/'combined40', COMBINED_SHA)]:
        if sha(source/'eunha2.exe') != expected:
            raise ValueError('Canonical or combined executable changed')
    STAGE.mkdir()
    helpers = ['audit_kepler_padding_reuse_pair_20261028.py',
               'audit_kepler_polar_long_gate_20261027.py',
               'audit_kepler_geometry_force_dump_20261029.py',
               'audit_kepler_vertex_halo_manual_20261029.py']
    for name in helpers:
        shutil.copy2(PARENT/name, STAGE/name)
    self_name = Path(__file__).name
    shutil.copy2(Path(__file__), STAGE/self_name)
    subprocess.run([sys.executable, '-c',
        'import audit_kepler_halo40_net_wall_20261029; import audit_kepler_padding_reuse_pair_20261028'],
        cwd=STAGE, check=True)
    inputs = ['eunha2.exe', 'params.dat', 'khout.000000.dat', 'khout.023008.dat',
              'audit_kepler_polar_long_gate_20261027.py', 'audit_kepler_native_radial_20261026.py']
    hashes = {}
    for name, source_name in CASES:
        source_root = ORIGINAL if source_name == 'qualifier_old' else PARENT
        source = source_root/source_name
        run = STAGE/name
        run.mkdir()
        for item in inputs:
            original = (source/item).resolve()
            (run/item).symlink_to(original)
            if original not in hashes:
                hashes[original] = sha(original)
        script = (source/'run.sh').read_text()
        script, n1 = re.subn(r'HYDRO_TSTOP=[^ ]+', 'HYDRO_TSTOP=2.66', script)
        script, n2 = re.subn(r'--time \S+ --output audit_endpoint_gate.json',
                            '--time 2.66 --output audit_endpoint_gate.json', script)
        script, n3 = re.subn(r'GFS_GEOM_FORCE_DUMP_MAXCALL=\d+',
                            'GFS_GEOM_FORCE_DUMP_MAXCALL=0', script)
        if (n1, n2, n3) != (1, 1, 1):
            raise ValueError('Unexpected frozen native wrapper')
        script = script.replace(str(source_root), str(STAGE))
        (run/'run.sh').write_text(script)
        subprocess.run(['bash', '-n', str(run/'run.sh')], check=True)
        (run/'inputs.sha256').write_text(''.join(
            f'{hashes[(run/item).resolve()]}  {item}\n' for item in inputs)
            +f'{sha(run/"run.sh")}  run.sh\n')
    batch = '#!/usr/bin/env bash\nset -euo pipefail\n'
    for name, _ in CASES:
        batch += f'cd {STAGE/name}\nbash run.sh > runner.out 2> runner.err\n'
    batch += f'cd {STAGE}\n/home/kjhan/miniconda3/bin/python3 {STAGE/self_name} --audit > comparison.log\n'
    (STAGE/'run_cases.sh').write_text(batch)
    wrapper = (PARENT/'run_manual.sh').read_text().replace(str(PARENT), str(STAGE))
    (STAGE/'run_manual.sh').write_text(wrapper)
    for name in ['run_cases.sh', 'run_manual.sh']:
        subprocess.run(['bash', '-n', str(STAGE/name)], check=True)
    files = helpers+[self_name, 'run_cases.sh', 'run_manual.sh']
    (STAGE/'audit_inputs.sha256').write_text(''.join(f'{sha(STAGE/name)}  {name}\n' for name in files))
    receipt = {'observed_at': datetime.datetime.now().astimezone().isoformat(),
        'root': str(STAGE), 'prerequisite': str(PARENT/'comparison.json'),
        'prerequisite_sha256': sha(PARENT/'comparison.json'), 'prerequisite_passed': True,
        'case_order': [name for name, _ in CASES], 'start': 2.6549183918680868, 'target': 2.66,
        'expected_steps_per_case': 59, 'geometry_force_dumps': 0,
        'halo': 40, 'vertex_certificate': 0, 'small_halo_promoted': False,
        'native_endpoint_and_physics_audits_retained': True,
        'GPU_UUID': GPU, 'OMP_threads': 2, 'MPI_ranks': 1, 'timeout_seconds': 2400,
        'output_estimated_GB': .34,
        'predeclared_clarity_rule': 'fastest canonical / slowest combined > 1.05',
        'production_promoted': False, 'existing_Slurm_jobs_modified': False,
        'source_sha256': {name: sha(STAGE/name) for name in files}}
    (STAGE/'staging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (Path(__file__).parent/'KEPLER_HALO40_NET_WALL_ADA_STAGE_2026-10-29.json').write_text(
        json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


def audit():
    root = Path(__file__).resolve().parent
    import audit_kepler_padding_reuse_pair_20261028 as common
    common.ROOT = root
    common.RANKS.update({name: 1 for name, _ in CASES})
    pairs = [common.compare('old_a', 'combined_b', exact=True),
             common.compare('combined_b', 'combined_b2', exact=True),
             common.compare('old_a', 'old_a2', exact=True)]
    checks = {'all_native_fields_traces_physics_exact': all(p['passed'] for p in pairs)}
    manifest = json.loads((root/'manual_execution.json').read_text())
    wall = {}
    for name, source in CASES:
        _, _, metrics = common.load(name)
        wall[name] = metrics['wall_seconds']
        context = json.loads((root/name/'launch_context.json').read_text())
        identity = ['execution_mode', 'hostname', 'GPU_UUID', 'CUDA_VISIBLE_DEVICES',
                    'manual_session_id', 'supervisor_pid']
        checks[name+'_same_actual_execution'] = (all(context.get(k) == manifest.get(k) for k in identity)
            and all(manifest.get(k) not in (None, '') for k in identity)
            and manifest.get('hostname') == 'LagEunha' and manifest.get('GPU_UUID') == GPU
            and manifest.get('execution_mode') == 'manual_dedicated_Ada'
            and not manifest.get('SLURM_JOB_ID') and not context.get('SLURM_JOB_ID'))
        expected = (0, 40, 0, 0) if source == 'qualifier_old' else (0, 40, 1, 1)
        keys = ['GFS_GPU_VERTEX_HALO_CERT', 'GFS_KEPLER_HALO_CELLS',
                'GFS_REUSE_TESS_PADDING', 'GFS_GPU_LOCAL_PERIODIC_SYNC']
        checks[name+'_flags_exact'] = all(context.get(k) == str(v) for k, v in zip(keys, expected))
        checks[name+'_binary_exact'] = sha(root/name/'eunha2.exe') == (BASE_SHA if source == 'qualifier_old' else COMBINED_SHA)
        checks[name+'_full_target_and_steps'] = abs(metrics['time']-2.66) < 1e-12 and metrics['steps'] == 59
        checks[name+'_no_geometry_dumps'] = not any((root/name).glob('geom_force_rank*.bin'))
        if source != 'qualifier_old':
            text = (root/name/'log').read_text()
            checks[name+'_all_combined_branches_executed'] = (
                bool(re.search(r'\[TESS_PADDING_REUSE\] rank=0 .*reuse=1 refresh_validation=1', text))
                and '[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh' in text
                and '[HALO_SYNC]' not in text and 'alias invariants failed' not in text
                and '[KEPLER_HALO] blend cells=40 cellsize=0.0625 width=2.5' in text)
    numerical = all(checks.values())
    result = {'scope': 'Four independent same-checkpoint native solver wall measurements in ABBA order',
              'checks': checks, 'numerical_passed': numerical, 'pairs': pairs,
              'performance': performance_summary(wall), 'production_promoted': False,
              'timing_scope': 'Bash time around mpirun includes native startup/evolution/endpoint I/O. '
                              'Post-run Python physics auditors are outside that interval. '
                              'Large geometry dumps are disabled; endpoint output remains enabled.',
              'auditor_sha256': sha(Path(__file__))}
    with (root/'net_wall_comparison.json').open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not numerical:
        raise SystemExit('Numerical or execution-context gate failed; evidence retained')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--stage', action='store_true')
    mode.add_argument('--audit', action='store_true')
    args = parser.parse_args()
    if args.stage:
        stage()
    else:
        audit()
