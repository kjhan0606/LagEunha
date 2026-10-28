"""Dedicated-Ada context adapter; native state/face/force gates are unchanged."""
import datetime
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import uuid

MODE = 'manual_dedicated_Ada'
FLAG_KEYS = ['GFS_GPU_VERTEX_HALO_CERT', 'GFS_KEPLER_HALO_CELLS',
             'GFS_REUSE_TESS_PADDING', 'GFS_GPU_LOCAL_PERIODIC_SYNC']


def context():
    mask = os.environ['CUDA_VISIBLE_DEVICES']
    gpu = subprocess.check_output(['nvidia-smi', '-i', mask, '--query-gpu=uuid',
                                   '--format=csv,noheader'], text=True, timeout=10).strip()
    session = os.environ['KEPLER_MANUAL_SESSION_ID']
    uuid.UUID(session)
    if socket.gethostname() != 'LagEunha' or not re.fullmatch(r'GPU-[0-9a-f-]{36}', gpu):
        raise RuntimeError('Expected actual dedicated Ada host and one GPU UUID')
    if mask != gpu or os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Expected explicit GPU UUID mask and honest non-Slurm context')
    return {'execution_mode': MODE, 'hostname': socket.gethostname(),
            'GPU_UUID': gpu, 'CUDA_VISIBLE_DEVICES': mask,
            'manual_session_id': session,
            'supervisor_pid': int(os.environ['KEPLER_MANUAL_SUPERVISOR_PID']),
            'SLURM_JOB_ID': os.environ.get('SLURM_JOB_ID', ''),
            'observed_at': datetime.datetime.now().astimezone().isoformat(),
            **{key: os.environ.get(key, '') for key in FLAG_KEYS}}


def launch_checks(root):
    from audit_kepler_vertex_halo_pair_20261029 import EXPECTED
    manifest_path = root/'manual_execution.json'
    checks = {'actual_manual_execution_record_present': manifest_path.is_file()}
    if not manifest_path.is_file():
        return checks
    manifest = json.loads(manifest_path.read_text())
    checks['manual_execution_mode_recorded'] = manifest.get('execution_mode') == MODE
    contexts = []
    for name, expected in EXPECTED.items():
        path = root/name/'launch_context.json'
        checks[name+'_launch_context_present'] = path.is_file()
        if not path.is_file():
            continue
        c = json.loads(path.read_text())
        contexts.append(c)
        checks[name+'_requested_flags_exact'] = all(c.get(k) == str(v) for k, v in zip(FLAG_KEYS, expected))
    checks['all_seven_cases_recorded'] = len(contexts) == len(EXPECTED)
    keys = ['execution_mode', 'hostname', 'GPU_UUID', 'CUDA_VISIBLE_DEVICES',
            'manual_session_id']
    checks['one_actual_host_GPU_session'] = bool(contexts) and all(
        manifest.get(key) not in (None, '') and all(c.get(key) == manifest[key] for c in contexts)
        for key in keys)
    checks['no_fabricated_Slurm_allocation'] = not manifest.get('SLURM_JOB_ID') and all(
        not c.get('SLURM_JOB_ID') for c in contexts)
    checks['dedicated_Ada_identity'] = (manifest.get('hostname') == 'LagEunha'
        and bool(re.fullmatch(r'GPU-[0-9a-f-]{36}', manifest.get('GPU_UUID', '')))
        and manifest.get('CUDA_VISIBLE_DEVICES') == manifest.get('GPU_UUID'))
    supervisors = {manifest.get('supervisor_pid')}
    continuation = root/'manual_continuation.json'
    if continuation.is_file():
        resumed = json.loads(continuation.read_text())
        checks['continuation_has_same_actual_host_GPU_session'] = all(
            resumed.get(key) == manifest.get(key) for key in keys) and not resumed.get('SLURM_JOB_ID')
        supervisors.add(resumed.get('supervisor_pid'))
    checks['all_actual_supervisors_recorded'] = bool(contexts) and all(
        isinstance(c.get('supervisor_pid'), int) and c['supervisor_pid'] > 0
        and c['supervisor_pid'] in supervisors for c in contexts)
    return checks


if __name__ == '__main__':
    if sys.argv[1:] in (['--record'], ['--record-execution'], ['--record-continuation']):
        record = context()
        execution = sys.argv[1:] != ['--record']
        continuation = sys.argv[1:] == ['--record-continuation']
        if execution:
            apps = subprocess.check_output(['nvidia-smi', '-i', record['GPU_UUID'],
                '--query-compute-apps=pid', '--format=csv,noheader'], text=True, timeout=10).strip()
            if apps:
                raise RuntimeError('Ada GPU already has compute processes; refusing overlap')
        if continuation:
            original = json.loads(Path('manual_execution.json').read_text())
            for key in ('execution_mode', 'hostname', 'GPU_UUID', 'manual_session_id'):
                if record[key] != original[key]:
                    raise RuntimeError('Continuation changed benchmark identity')
            try:
                os.kill(original['supervisor_pid'], 0)
            except ProcessLookupError:
                pass
            else:
                raise RuntimeError('Original supervisor still exists; refusing overlap')
        output = Path('manual_continuation.json' if continuation else
                      ('manual_execution.json' if execution else 'launch_context.json'))
        with output.open('x') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')
    else:
        import audit_kepler_vertex_halo_pair_20261029 as native
        native.launch_checks = launch_checks
        native.main()
