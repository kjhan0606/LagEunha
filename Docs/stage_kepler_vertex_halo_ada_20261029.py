"""Stage the frozen seven native comparisons on the previously authorized Ada host."""
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

DOCS = Path(__file__).resolve().parent
Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
SOURCE = Q1/'vertex_halo_pair_t2654918_a100_20261029'
ROOT = Q1/'vertex_halo_pair_t2654918_ada_20261029'
CASES = ['qualifier_old', 'qualifier_new', 'control40', 'vertex40',
         'vertex20', 'vertex20_reuse', 'vertex20_reuse_local_sync']
GPU = 'GPU-42b3367f-1df7-2b2d-7961-7781c87b022d'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    subprocess.run(['sha256sum', '-c', 'audit_inputs.sha256'], cwd=SOURCE, check=True,
                   stdout=subprocess.DEVNULL)
    expected = {'qualifier_old': '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0',
                'qualifier_new': 'a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc'}
    for name, value in expected.items():
        if sha((SOURCE/name/'eunha2.exe').resolve()) != value:
            raise ValueError('Unexpected frozen native binary')
    ROOT.mkdir()
    helpers = ['audit_kepler_vertex_halo_pair_20261029.py',
               'audit_kepler_geometry_force_dump_20261029.py',
               'audit_kepler_padding_reuse_pair_20261028.py',
               'audit_kepler_polar_long_gate_20261027.py']
    for name in helpers:
        shutil.copy2(SOURCE/name, ROOT/name)
    adapter = 'audit_kepler_vertex_halo_manual_20261029.py'
    shutil.copy2(DOCS/adapter, ROOT/adapter)
    subprocess.run([sys.executable, '-c',
        'import audit_kepler_vertex_halo_manual_20261029; import audit_kepler_vertex_halo_pair_20261029'],
        cwd=ROOT, check=True)
    inputs = ['eunha2.exe', 'params.dat', 'khout.000000.dat', 'khout.023008.dat',
              'audit_kepler_polar_long_gate_20261027.py', 'audit_kepler_native_radial_20261026.py']
    hashes = {}
    for case in CASES:
        run = ROOT/case
        run.mkdir()
        for name in inputs:
            source = (SOURCE/case/name).resolve()
            (run/name).symlink_to(source)
            if source not in hashes:
                hashes[source] = sha(source)
        script = (SOURCE/case/'run.sh').read_text()
        old = 'export I_MPI_HYDRA_BOOTSTRAP=slurm'
        if script.count(old) != 1:
            raise ValueError('Unexpected MPI bootstrap declaration')
        script = script.replace(old, 'unset I_MPI_HYDRA_BOOTSTRAP')
        old = str(SOURCE/'record_kepler_vertex_case_context_20261029.py')
        if script.count(old) != 1:
            raise ValueError('Unexpected native context recorder')
        script = script.replace(old, str(ROOT/adapter)+' --record')
        (run/'run.sh').write_text(script)
        subprocess.run(['bash', '-n', str(run/'run.sh')], check=True)
        (run/'inputs.sha256').write_text(''.join(
            f'{hashes[(run/name).resolve()]}  {name}\n' for name in inputs)
            +f'{sha(run/"run.sh")}  run.sh\n')
    batch = (SOURCE/'job_with_combinations.slurm').read_text()
    batch = '\n'.join(line for line in batch.splitlines() if not line.startswith('#SBATCH'))+'\n'
    batch = batch.replace(str(SOURCE), str(ROOT))
    batch = batch.replace('audit_kepler_vertex_halo_pair_20261029.py', adapter)
    (ROOT/'run_cases.sh').write_text(batch)
    wrapper = f'''#!/usr/bin/env bash
set -euo pipefail
cd {ROOT}
[[ "$(hostname)" == LagEunha ]]
[[ -z "${{SLURM_JOB_ID:-}}" ]]
export CUDA_VISIBLE_DEVICES={GPU}
export KEPLER_MANUAL_SESSION_ID=$(/home/kjhan/miniconda3/bin/python3 -c 'import uuid;print(uuid.uuid4())')
export KEPLER_MANUAL_SUPERVISOR_PID=$$
sha256sum -c audit_inputs.sha256
/home/kjhan/miniconda3/bin/python3 {ROOT/adapter} --record-execution
exec bash run_cases.sh
'''
    (ROOT/'run_manual.sh').write_text(wrapper)
    for name in ['run_cases.sh', 'run_manual.sh']:
        subprocess.run(['bash', '-n', str(ROOT/name)], check=True)
    files = helpers+[adapter, 'run_cases.sh', 'run_manual.sh']
    (ROOT/'audit_inputs.sha256').write_text(''.join(f'{sha(ROOT/name)}  {name}\n' for name in files))
    receipt = {'observed_at': datetime.datetime.now().astimezone().isoformat(),
        'root': str(ROOT), 'source_frozen_batch': str(SOURCE), 'cases': CASES,
        'host': 'LagEunha', 'GPU_UUID': GPU, 'execution_mode': 'manual_dedicated_Ada',
        'parent_checkpoint_sha256': hashes[(ROOT/'control40/khout.023008.dat').resolve()],
        'canonical_binary_sha256': expected['qualifier_old'],
        'candidate_binary_sha256': expected['qualifier_new'],
        'changes_from_A100_cases': ['local MPI bootstrap', 'actual manual host/GPU/session recorder'],
        'physical_geometry_force_and_timestep_gates_changed': False,
        'outputs_estimated_GB': 1.8, 'evolution_processes_at_once': 1, 'OMP_threads': 2,
        'execution_timeout_seconds': 2400, 'existing_Slurm_jobs_modified': False,
        'source_hashes': {name: sha(ROOT/name) for name in files}}
    (ROOT/'staging_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    (DOCS/'KEPLER_VERTEX_HALO_ADA_STAGE_2026-10-29.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
