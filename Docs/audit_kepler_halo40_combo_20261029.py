"""Exact halo40 reuse+alias fallback; retain the failed halo20 result unchanged."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
PARENT = Q1/'vertex_halo_pair_t2654918_ada_20261029'
STAGE = Q1/'halo40_reuse_alias_t2654918_ada_20261029'
CANDIDATE = 'a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc'
GPU = 'GPU-42b3367f-1df7-2b2d-7961-7781c87b022d'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage():
    prior = json.loads((PARENT/'comparison.json').read_text())
    if not prior['qualifier_passed'] or not prior['pairs'][0]['passed']:
        raise ValueError('Missing canonical qualifier and exact halo40 reference')
    if (PARENT/'manual_resume_exit_status.txt').read_text().strip() != '1':
        raise ValueError('Expected terminated original seven-case execution')
    if sha(PARENT/'control40/eunha2.exe') != CANDIDATE:
        raise ValueError('Unexpected reference executable')
    STAGE.mkdir()
    helpers = ['audit_kepler_padding_reuse_pair_20261028.py',
               'audit_kepler_polar_long_gate_20261027.py',
               'audit_kepler_geometry_force_dump_20261029.py',
               'audit_kepler_vertex_halo_manual_20261029.py']
    for name in helpers:
        shutil.copy2(PARENT/name, STAGE/name)
    name = Path(__file__).name
    shutil.copy2(Path(__file__), STAGE/name)
    subprocess.run([sys.executable, '-c', 'import audit_kepler_halo40_combo_20261029; '
                    'import audit_kepler_padding_reuse_pair_20261028'], cwd=STAGE, check=True)
    (STAGE/'control40').symlink_to(PARENT/'control40', target_is_directory=True)
    source = PARENT/'vertex20_reuse_local_sync'
    run = STAGE/'combined40'
    run.mkdir()
    inputs = ['eunha2.exe', 'params.dat', 'khout.000000.dat', 'khout.023008.dat',
              'audit_kepler_polar_long_gate_20261027.py', 'audit_kepler_native_radial_20261026.py']
    for item in inputs:
        (run/item).symlink_to((source/item).resolve())
    script = (source/'run.sh').read_text()
    for old, new in [('GFS_KEPLER_HALO_CELLS=20', 'GFS_KEPLER_HALO_CELLS=40'),
                     ('GFS_GPU_VERTEX_HALO_CERT=1', 'GFS_GPU_VERTEX_HALO_CERT=0')]:
        if script.count(old) != 1:
            raise ValueError('Unexpected native wrapper')
        script = script.replace(old, new)
    script = script.replace(str(PARENT), str(STAGE))
    (run/'run.sh').write_text(script)
    (run/'inputs.sha256').write_text(''.join(f'{sha(run/item)}  {item}\n' for item in inputs)
                                   +f'{sha(run/"run.sh")}  run.sh\n')
    wrapper = (PARENT/'run_manual.sh').read_text().replace(str(PARENT), str(STAGE))
    (STAGE/'run_manual.sh').write_text(wrapper)
    batch = f'''#!/usr/bin/env bash
set -euo pipefail
cd {run}
bash run.sh > runner.out 2> runner.err
cd {STAGE}
/home/kjhan/miniconda3/bin/python3 {STAGE/name} --audit > comparison.log
'''
    (STAGE/'run_cases.sh').write_text(batch)
    for item in ['run_manual.sh', 'run_cases.sh', 'combined40/run.sh']:
        subprocess.run(['bash', '-n', str(STAGE/item)], check=True)
    files = helpers+[name, 'run_manual.sh', 'run_cases.sh']
    (STAGE/'audit_inputs.sha256').write_text(''.join(f'{sha(STAGE/item)}  {item}\n' for item in files))
    receipt = {'observed_at': datetime.datetime.now().astimezone().isoformat(),
               'root': str(STAGE), 'reference_reused': str(PARENT/'control40'),
               'reference_replayed': False, 'new_native_cases': 1,
               'checkpoint': 2.6549183918680868, 'target': 2.66,
               'candidate_sha256': CANDIDATE, 'halo': 40,
               'vertex_certificate': 0, 'padding_reuse': 1, 'local_alias': 1,
               'geometry_force_dumps': 4, 'output_estimated_GB': .31,
               'halo20_seven_case_result_preserved_as_failed': True,
               'physical_and_exact_geometry_gates_relaxed': False,
               'existing_Slurm_and_main_modified': False,
               'source_sha256': {item: sha(STAGE/item) for item in files}}
    (STAGE/'staging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (Path(__file__).parent/'KEPLER_HALO40_COMBO_ADA_STAGE_2026-10-29.json').write_text(
        json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


def audit():
    root = Path(__file__).resolve().parent
    import audit_kepler_padding_reuse_pair_20261028 as common
    from audit_kepler_geometry_force_dump_20261029 import compare_dumps
    common.ROOT = root
    common.RANKS.update(control40=1, combined40=1)
    state = common.compare('control40', 'combined40', exact=True)
    geometry = [compare_dumps(root/'control40'/f'geom_force_rank000_call{i:03d}.bin',
                              root/'combined40'/f'geom_force_rank000_call{i:03d}.bin') for i in range(1, 5)]
    context = json.loads((root/'combined40/launch_context.json').read_text())
    baseline = json.loads((root/'control40/launch_context.json').read_text())
    manifest = json.loads((root/'manual_execution.json').read_text())
    text = (root/'combined40/log').read_text()
    keys = ['GFS_GPU_VERTEX_HALO_CERT', 'GFS_KEPLER_HALO_CELLS',
            'GFS_REUSE_TESS_PADDING', 'GFS_GPU_LOCAL_PERIODIC_SYNC']
    checks = {'full_59_step_state_and_trace_exact': state['passed'] and all(m['steps'] == 59 for m in state['metrics']),
              'four_geometry_and_force_calls_exact': all(g['passed'] for g in geometry),
              'requested_flags_exact': all(context.get(k) == str(v) for k, v in zip(keys, (0, 40, 1, 1))),
              'same_actual_reference_host_GPU': all(c.get('hostname') == 'LagEunha' and c.get('GPU_UUID') == GPU
                                                    for c in [context, baseline, manifest]),
              'honest_current_execution': all(context.get(k) == manifest.get(k) and manifest.get(k) not in (None, '')
                  for k in ['manual_session_id', 'supervisor_pid', 'CUDA_VISIBLE_DEVICES'])
                  and not context.get('SLURM_JOB_ID') and not manifest.get('SLURM_JOB_ID'),
              'correct_native_binary': sha(root/'combined40/eunha2.exe') == CANDIDATE,
              'padding_reuse_executed': bool(re.search(r'\[TESS_PADDING_REUSE\] rank=0 .*reuse=1 refresh_validation=1', text)),
              'local_alias_executed': '[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh' in text,
              'no_host_alias_fallback': '[HALO_SYNC]' not in text and 'alias invariants failed' not in text,
              'halo40_executed': '[KEPLER_HALO] blend cells=40 cellsize=0.0625 width=2.5' in text}
    result = {'passed': all(checks.values()), 'checks': checks, 'state': state,
              'four_force_calls': geometry, 'production_promoted': False,
              'scope': 'Exact halo40 reuse+alias implementation check, not halo20 acceptance, net speed or long-time stability',
              'auditor_sha256': sha(Path(__file__))}
    with (root/'comparison.json').open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    if not result['passed']:
        raise SystemExit('Halo40 combination gate failed; evidence retained')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stage', action='store_true')
    group.add_argument('--audit', action='store_true')
    args = parser.parse_args()
    if args.stage:
        stage()
    else:
        audit()
