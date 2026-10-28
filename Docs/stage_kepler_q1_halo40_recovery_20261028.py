"""Stage a guarded continuation from the accepted t=2.654918 RK4 start."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

DOCS = Path(__file__).resolve().parent
Q1 = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
PARENT = Q1 / 'diag_clip_swap_t2_to65_a100_20261028'
ROOT = Q1 / 'halo40_from_t2654918_a100_20261028'
START = 2.6549183918680868
STEP = 23008
PYTHON = '/home/kjhan/miniconda3/bin/python3'
EXPECTED_BINARY = '6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    if ROOT.exists():
        raise FileExistsError(ROOT)
    sys.path.insert(0, str(DOCS))
    from audit_kepler_polar_long_gate_20261027 import ordered
    import numpy as np
    checkpoint = PARENT / 'khdiag.step023008.start_geometry.dat'
    baseline = Q1 / 'dteta_no_cap_from_t1p5_to65_20261028/khdiag.step023008.start_geometry.dat'
    binary = (PARENT / 'eunha2.exe').resolve()
    if digest(binary) != EXPECTED_BINARY:
        raise ValueError('Unexpected CUDA candidate binary')
    candidate, tc = ordered(checkpoint)
    reference, tr = ordered(baseline)
    if abs(tc - START) > 1e-12 or abs(tr - START) > 1e-12:
        raise ValueError('Unexpected accepted step-start time')
    fields = ('x', 'y', 'vx', 'vy', 'mass', 'volume', 'ie', 'den', 'pressure', 'w2')
    parity = {f: bool(np.array_equal(candidate[f], reference[f])) for f in fields}
    differences = {f: {'different_cells': int(np.count_nonzero(candidate[f] != reference[f])),
                       'max_absolute_difference': float(np.max(np.abs(candidate[f].astype(float) - reference[f].astype(float))))}
                   for f in fields}
    # The measured cross-device difference is five vx entries, at most two
    # float64 ULPs (1.11e-16); retain exact thermal/geometric/mass equality.
    velocity_roundoff_only = all(bool(np.all(
        np.abs(candidate[f] - reference[f]) <= 2 * np.abs(np.spacing(reference[f]))))
        for f in ('vx', 'vy'))
    if not all(parity[f] for f in fields if f not in ('vx', 'vy')) or not velocity_roundoff_only:
        raise ValueError(f'Accepted state differs beyond velocity roundoff: {differences}')

    ROOT.mkdir()
    preflight = ROOT / 'accepted_start_audit'
    preflight.mkdir()
    (preflight / 'khout.023008.dat').symlink_to(checkpoint)
    (preflight / 'log').symlink_to(PARENT / 'log')
    subprocess.run([PYTHON, str(DOCS / 'audit_kepler_polar_long_gate_20261027.py'),
                    str(Q1 / 'khout.000000.dat'), str(preflight), '--time', str(START),
                    '--output', str(preflight / 'audit_gate.json')], check=True,
                   stdout=(preflight / 'audit.log').open('x'))
    for name in ('audit_kepler_polar_long_gate_20261027.py',
                 'audit_kepler_native_radial_20261026.py'):
        shutil.copyfile(DOCS / name, ROOT / name)
    inputs = {
        'eunha2.exe': binary,
        'params.dat': (PARENT / 'params.dat').resolve(),
        'khout.000000.dat': (Q1 / 'khout.000000.dat').resolve(),
    }
    hashes = {name: digest(path) for name, path in inputs.items()}
    targets = [3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0]
    stages = [ROOT / ('to' + str(t).replace('.', 'p')) for t in targets]
    for i, (target, stage) in enumerate(zip(targets, stages)):
        stage.mkdir()
        for name, source in inputs.items():
            (stage / name).symlink_to(source)
        (stage / 'inputs.sha256').write_text(''.join(f'{sha}  {name}\n' for name, sha in hashes.items()))
        manifest = {
            'purpose': 'Continue the existing accepted CUDA+MPI+OpenMP W0 trajectory with a wider certified periodic search region',
            'single_change': 'GFS_KEPLER_HALO_CELLS=28 -> 40; force and timestep parameters retained',
            'initial_time': START if i == 0 else targets[i - 1],
            'target_time': target,
            'checkpoint': str(checkpoint) if i == 0 else None,
            'checkpoint_step': STEP if i == 0 else None,
            'checkpoint_sha256': digest(checkpoint) if i == 0 else None,
            'parent_gate': str(preflight / 'audit_gate.json') if i == 0 else str(stages[i - 1] / 'audit_gate.json'),
            'next_job': str(stages[i + 1] / 'job.slurm') if i + 1 < len(stages) else None,
            'binary_sha256': hashes['eunha2.exe'],
            'initial_sha256': hashes['khout.000000.dat'],
            'physics': 'fixed mass, W=0, no generator steering, no added viscosity or spin, original RK4 force/energy bookkeeping',
            'output_policy': '85 MB native checkpoints every 0.5 code-time units and exact segment endpoint; about 1.5 GB through t=10',
        }
        (stage / 'stage.json').write_text(json.dumps(manifest, indent=2) + '\n')
        script = f'''#!/usr/bin/env bash
#SBATCH --job-name=KepQ1H40T{target:g}
#SBATCH --partition=a100,a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --chdir={stage}
#SBATCH --output={stage}/slurm-%j.out
#SBATCH --error={stage}/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
PY={PYTHON}
restart_step=$($PY {ROOT}/prepare_restart.py)
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${{LD_LIBRARY_PATH:-}}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP={target} EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=40
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=0.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset EUNHA_GPU_DT_DIAG EUNHA_GPU_TESS_PROFILE GFS_GPU_TESS_PROFILE EUNHA_GPU_FACE_AREA_AUDIT
unset GFS_VARPRESS_GEOM_AUDIT GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
rc=$?
set -e
printf 'EXIT:%s\\n' "$rc" >> log
[[ "$rc" -eq 0 ]] || exit "$rc"
$PY {ROOT}/audit_kepler_polar_long_gate_20261027.py khout.000000.dat "$PWD" --time {target} --output audit_gate.json > audit_gate.log
endpoint=$($PY -c 'import json; print(json.load(open("audit_gate.json"))["endpoint"])')
$PY {ROOT}/audit_kepler_native_radial_20261026.py khout.000000.dat "$endpoint" --mode baseline --output audit_radial.json > audit_radial.log
$PY {ROOT}/check_radial_and_next.py
next_job=$($PY -c 'import json; print(json.load(open("stage.json"))["next_job"] or "")')
if [[ -n "$next_job" ]]; then
  sbatch --parsable --dependency=afterok:${{SLURM_JOB_ID}} "$next_job" > next_job_id.txt
fi
'''
        (stage / 'job.slurm').write_text(script)
        subprocess.run(['bash', '-n', str(stage / 'job.slurm')], check=True)

    (ROOT / 'prepare_restart.py').write_text('''import hashlib, json, re, struct
from pathlib import Path
m=json.loads(Path('stage.json').read_text())
gate=json.loads(Path(m['parent_gate']).read_text())
assert all(gate['checks'].values()), 'parent physics gate failed'
assert abs(gate['time']-m['initial_time']) < 1e-10, 'parent time mismatch'
for output in ('log','audit_gate.json','audit_radial.json','restart_receipt.json'):
    assert not Path(output).exists(), 'existing run output: '+output
checkpoint=Path(m['checkpoint']) if m['checkpoint'] else Path(gate['endpoint'])
assert checkpoint.is_absolute()
with checkpoint.open('rb') as stream:
    magic,n,stride,stress=struct.unpack('<4i',stream.read(16))
    stream.seek(-52,2)
    time=struct.unpack('<d',stream.read(8))[0]
assert (magic,n,stride,stress)==(-1,147456,576,1)
assert checkpoint.stat().st_size==16+n*stride+52
assert abs(time-m['initial_time']) < 1e-10
h=hashlib.sha256()
with checkpoint.open('rb') as stream:
    for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
if m['checkpoint_sha256']: assert h.hexdigest()==m['checkpoint_sha256']
step=m['checkpoint_step'] if m['checkpoint_step'] else int(re.fullmatch(r'khout\\.(\\d+)\\.dat',checkpoint.name)[1])
link=Path(f'khout.{step:06d}.dat')
assert not link.exists() and not link.is_symlink()
link.symlink_to(checkpoint)
Path('restart_receipt.json').write_text(json.dumps({'checkpoint':str(checkpoint),'sha256':h.hexdigest(),'time':time,'step':step},indent=2)+'\\n')
print(step)
''')
    (ROOT / 'check_radial_and_next.py').write_text('''import json
from pathlib import Path
gate=json.loads(Path('audit_gate.json').read_text())
checks={}
for row in gate['bins']:
    if (row['rlo'],row['rhi']) in ((.5,.75),(.75,1.),(1.,2.)):
        label=f"{row['rlo']}_{row['rhi']}"
        checks[label+'_density_below_1pct']=row['density_abs_mass_weighted'] < .01
        checks[label+'_entropy_log_below_0p01']=row['entropy_log_abs_mass_weighted'] < .01
        checks[label+'_position_error_below_0p01']=row['position_error_abs_mass_weighted'] < .01
Path('continuation_gate.json').write_text(json.dumps({'time':gate['time'],'checks':checks},indent=2)+'\\n')
assert len(checks)==9 and all(checks.values()), 'radial continuation gate failed'
''')
    receipt = {
        'root': str(ROOT), 'accepted_time': START, 'restart_step': STEP,
        'accepted_start_physical_fields_bitwise_equal_to_baseline': parity,
        'accepted_start_field_differences': differences,
        'velocity_differences_within_two_float64_ulps': velocity_roundoff_only,
        'accepted_checkpoint_sha256': digest(checkpoint),
        'binary_sha256': EXPECTED_BINARY, 'halo_cells_before': 28, 'halo_cells_after': 40,
        'targets': targets, 'initial_job': str(stages[0] / 'job.slurm'),
        'scheduler': 'one A100; a100 or a100_pcie; submit next segment only after full and radial gates pass',
        'no_long_time_success_claim': True,
    }
    (ROOT / 'staging_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (DOCS / 'KEPLER_Q1_HALO40_RECOVERY_STAGE_2026-10-28.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
