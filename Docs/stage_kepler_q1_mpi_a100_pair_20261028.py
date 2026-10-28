"""Run the existing W0 MPI rank comparison on available A100 GPUs."""
import hashlib
import json
from pathlib import Path
import subprocess

Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
SOURCE = Q1 / 'mpi_scale_w0_t15_to155_h200_20261028/np2_after438363_20261028'
ROOT = Q1 / 'mpi_scale_w0_to155_a100_pair_20261028'
DOCS = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    if ROOT.exists():
        raise FileExistsError(ROOT)
    inputs = {}
    for line in (SOURCE / 'inputs.sha256').read_text().splitlines():
        expected, name = line.split()
        if digest(SOURCE / name) != expected:
            raise ValueError(f'Source input hash mismatch: {name}')
        inputs[name] = expected
    source_script = (SOURCE / 'job.slurm').read_text()
    if source_script.count('mpirun -np 2') != 1:
        raise ValueError('Expected one native MPI launch')
    ROOT.mkdir()
    for ranks in (1, 2):
        run = ROOT / f'np{ranks}'
        run.mkdir()
        for name in inputs:
            (run / name).symlink_to(SOURCE / name)
        (run / 'inputs.sha256').write_text((SOURCE / 'inputs.sha256').read_text())
        runner = '\n'.join(line for line in source_script.splitlines()
                           if not line.startswith('#SBATCH')) + '\n'
        runner = runner.replace('mpirun -np 2', f'mpirun -np {ranks}')
        (run / 'run.sh').write_text(runner)
        subprocess.run(['bash', '-n', str(run / 'run.sh')], check=True)
    comparison = '''import json, re, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs')
from audit_kepler_polar_long_gate_20261027 import ordered
root=Path(__file__).resolve().parent
audits=[json.loads((root/f'np{n}/audit_t1p55_gate.json').read_text()) for n in (1,2)]
states=[ordered(Path(g['endpoint'])) for g in audits]
a,ta=states[0]; b,tb=states[1]
checks={'both_native_physics_gates':all(all(g['checks'].values()) for g in audits),
        'same_exact_time':abs(ta-1.55)<1e-12 and abs(tb-1.55)<1e-12}
differences={}
for f in ('x','y','vx','vy','mass','volume','ie','den','pressure','w2'):
    av=a[f].astype(float); bv=b[f].astype(float); delta=bv-av
    scale=max(float(np.sqrt(np.mean(av*av))),1e-300)
    differences[f]={'bitwise_equal':bool(np.array_equal(a[f],b[f])),
                    'different_cells':int(np.count_nonzero(a[f]!=b[f])),
                    'max_absolute_difference':float(np.max(abs(delta))),
                    'rms_difference':float(np.sqrt(np.mean(delta*delta))),
                    'normalized_rms_difference':float(np.sqrt(np.mean(delta*delta)))/scale}
    if f in ('mass','w2'): checks[f+'_exact']=differences[f]['bitwise_equal']
    if f in ('x','y'): checks[f+'_max_abs_below_1e_6']=differences[f]['max_absolute_difference']<1e-6
    if f in ('vx','vy'): checks[f+'_max_abs_below_1e_5']=differences[f]['max_absolute_difference']<1e-5
    if f in ('volume','ie','den','pressure'):
        relative=float(np.max(abs(delta)/np.maximum(abs(av),1e-300)))
        differences[f]['max_cell_relative_difference']=relative
        checks[f+'_normalized_rms_below_1e_6']=differences[f]['normalized_rms_difference']<1e-6
        checks[f+'_max_cell_relative_below_1e_3']=relative<1e-3
times=[]
for n in (1,2):
    matches=re.findall(r'\\[RUN_TIME\\] wall=([0-9.]+)',(root/f'np{n}/runner.err').read_text())
    assert len(matches)==1
    times.append(float(matches[0]))
result={'scope':'short same-GPU-model MPI rank consistency and throughput; no long-time claim',
        'time':ta,'wall_seconds_np1_np2':times,'wall_speedup_np1_over_np2':times[0]/times[1],
        'checks':checks,'field_differences':differences,
        'predeclared_tolerances':{'position_max_abs':1e-6,'velocity_max_abs':1e-5,
                                 'positive_scalar_normalized_rms':1e-6,'positive_scalar_max_cell_relative':1e-3}}
(root/'pair_audit.json').write_text(json.dumps(result,indent=2)+'\\n')
print(json.dumps(result,indent=2))
if not all(checks.values()): raise SystemExit('MPI rank consistency gate failed')
'''
    (ROOT / 'compare.py').write_text(comparison)
    wrapper = f'''#!/usr/bin/env bash
#SBATCH --job-name=KepMPIA100Pair
#SBATCH --partition=a100,a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --chdir={ROOT}
#SBATCH --output={ROOT}/slurm-%j.out
#SBATCH --error={ROOT}/slurm-%j.err
set -euo pipefail
for ranks in 1 2; do
  cd {ROOT}/np${{ranks}}
  bash run.sh > runner.out 2> runner.err
done
/home/kjhan/miniconda3/bin/python3 {ROOT}/compare.py > {ROOT}/pair_audit.log
'''
    (ROOT / 'job.slurm').write_text(wrapper)
    subprocess.run(['bash', '-n', str(ROOT / 'job.slurm')], check=True)
    receipt = {'root': str(ROOT), 'source': str(SOURCE), 'inputs_sha256': inputs,
               'ranks': [1, 2], 'gpu_model': 'A100', 'threads_per_rank': 2,
               'target_time': 1.55, 'initial_job': str(ROOT / 'job.slurm'),
               'purpose': 'Measure existing native CUDA+MPI evolution without waiting for H200 partition QoS',
               'original_H200_job': 438394, 'original_job_preserved': True,
               'output_policy': 'one final native checkpoint per rank branch; about 170 MB total',
               'pair_tolerances_declared_before_execution': True}
    (ROOT / 'staging_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (DOCS / 'KEPLER_Q1_MPI_A100_PAIR_STAGE_2026-10-28.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
