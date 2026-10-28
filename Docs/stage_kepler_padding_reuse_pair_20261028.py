"""Prepare one focused job to qualify and measure the padding-reuse change."""
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

DOCS = Path(__file__).resolve().parent
Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
ROOT = Q1 / 'tess_padding_reuse_pair_20261028'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    from audit_kepler_polar_long_gate_20261027 import ordered
    build = json.loads((ROOT / 'build_receipt.json').read_text())
    assert build['GPU_and_RK4_object_bytes_unchanged']
    assert digest(ROOT / 'eunha2.reuse.exe') == build['candidate_sha256']
    parent = json.loads((Q1 / 'halo40_from_t2654918_a100_20261028/accepted_start_audit/audit_gate.json').read_text())
    assert all(parent['checks'].values())
    checkpoint = Q1 / 'diag_clip_swap_t2_to65_a100_20261028/khdiag.step023008.start_geometry.dat'
    state, start = ordered(checkpoint)
    assert abs(start - 2.6549183918680868) < 1e-12
    assert digest(checkpoint) == 'c34a846dedaf2e01a8e5548de2aa091f98db53f761424ce6cd931107ff0fc4a4'
    source = Q1 / 'mpi_scale_w0_to155_a100_pair_20261028/np2'
    source_script = (source / 'run.sh').read_text()
    old_binary = Q1 / 'diag_clip_swap_ab_t2_a100_pcie_20261028/candidate_binary.exe'
    assert digest(old_binary) == build['baseline_binary_sha256']
    cases = [('qualifier_old', 1, 0, start + 1e-5, old_binary),
             ('qualifier_new', 1, 0, start + 1e-5, ROOT / 'eunha2.reuse.exe'),
             ('np1_control', 1, 0, 2.66, ROOT / 'eunha2.reuse.exe'),
             ('np2_control', 2, 0, 2.66, ROOT / 'eunha2.reuse.exe'),
             ('np2_reuse', 2, 1, 2.66, ROOT / 'eunha2.reuse.exe'),
             ('np1_reuse', 1, 1, 2.66, ROOT / 'eunha2.reuse.exe')]
    for name, ranks, reuse, target, binary in cases:
        run = ROOT / name
        run.mkdir()
        files = {'eunha2.exe': binary, 'params.dat': source / 'params.dat',
                 'khout.000000.dat': source / 'khout.000000.dat', 'khout.023008.dat': checkpoint,
                 'audit_kepler_polar_long_gate_20261027.py': source / 'audit_kepler_polar_long_gate_20261027.py',
                 'audit_kepler_native_radial_20261026.py': source / 'audit_kepler_native_radial_20261026.py'}
        for dest, original in files.items():
            (run / dest).symlink_to(original)
        (run / 'inputs.sha256').write_text(''.join(f'{digest(path)}  {dest}\n' for dest, path in files.items()))
        script = source_script.replace('mpirun -np 2', f'mpirun -np {ranks}')
        script = script.replace('params.dat 13596', 'params.dat 23008')
        script = script.replace('HYDRO_TSTOP=1.55', f'HYDRO_TSTOP={target:.17g}')
        script = script.replace('--time 1.55', f'--time {target:.17g}')
        script = script.replace('GFS_KEPLER_HALO_CELLS=28', 'GFS_KEPLER_HALO_CELLS=40')
        script = script.replace('t1p55', 'endpoint')
        script = script.replace("TIMEFORMAT=", f'export GFS_REUSE_TESS_PADDING={reuse}\nTIMEFORMAT=', 1)
        (run / 'run.sh').write_text(script)
        subprocess.run(['bash', '-n', str(run / 'run.sh')], check=True)
    for name in ('audit_kepler_padding_reuse_pair_20261028.py', 'audit_kepler_polar_long_gate_20261027.py'):
        ast.parse((DOCS / name).read_text())
        shutil.copy2(DOCS / name, ROOT / name)
    job = f'''#!/usr/bin/env bash
#SBATCH --job-name=KepPadReusePair
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:2
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --chdir={ROOT}
#SBATCH --output={ROOT}/slurm-%j.out
#SBATCH --error={ROOT}/slurm-%j.err
set -euo pipefail
for name in qualifier_old qualifier_new np1_control np2_control np2_reuse np1_reuse; do
  cd {ROOT}/$name
  bash run.sh > runner.out 2> runner.err
  if [[ "$name" == qualifier_new ]]; then
    /home/kjhan/miniconda3/bin/python3 {ROOT}/audit_kepler_padding_reuse_pair_20261028.py --qualifier > {ROOT}/qualifier.log
  fi
done
/home/kjhan/miniconda3/bin/python3 {ROOT}/audit_kepler_padding_reuse_pair_20261028.py > {ROOT}/comparison.log
'''
    (ROOT / 'job.slurm').write_text(job)
    subprocess.run(['bash', '-n', str(ROOT / 'job.slurm')], check=True)
    receipt = {'purpose': 'Single targeted validation job after poor MPI scaling; no speculative queued chain',
               'checkpoint': str(checkpoint), 'start_time': start, 'restart_step': 23008,
               'endpoint': 2.66, 'halo_cells': 40, 'ranks': [1, 2],
               'qualifier': 'Old canonical executable versus new executable with reuse disabled; all ten fields exact',
               'operator_gate': 'Same new executable with reuse off/on; exact ten fields and accepted-step trace per rank count',
               'physics_gate': 'Existing 15 checks per run, strict stage positivity and zero injection in every ledger',
               'storage_policy': 'Six exact endpoint snapshots, about 510 MB plus binaries; existing parent files linked',
               'build': build, 'cases': [{'name': n, 'ranks': r, 'reuse': e, 'target': t} for n, r, e, t, _ in cases],
               'native_simulator': 'Compiled CUDA+MPI+OpenMP eunha2; Python only stages and audits',
               'existing_long_run_untouched': True}
    (ROOT / 'staging_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (DOCS / 'KEPLER_PADDING_REUSE_STAGE_2026-10-28.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'root': str(ROOT), 'job': str(ROOT / 'job.slurm'), 'start': start, 'stop': 2.66, 'cases': len(cases)}, indent=2))


if __name__ == '__main__':
    main()
