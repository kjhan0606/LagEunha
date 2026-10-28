"""Resume the accepted weighted step-start checkpoint with the v5 build."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS = Path(__file__).resolve().parent
OLD = ROOT / 'gpu_energy_neutral_zero_pressure_v4_fresh_jitter64_to65_a10_v1'
BUILD = ROOT / 'gpu_energy_neutral_capacity64_build_v5'
RUN = ROOT / 'gpu_energy_neutral_capacity64_v5_recover_t005_a10_v1'
RECEIPT = DOCS / 'KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/accepted_step_start_audit.json'
CHECKPOINT = OLD / 'khdiag.step000024.start_geometry.dat'
CHECKPOINT_SHA = '35c9e11ff60471e7a68a991d16acd5b9571926cceb8e17c8d310141dc9039c7e'
INITIAL_SHA = '8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'
MASK_SHA = 'b560aca44c81e4560f8371c6d979e1afbe853f8cd18b592d0fc19a58becc8f4c'
OLD_BINARY_SHA = '3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc'
BINARY_SHA = 'feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'
START = 0.045616533734626756


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    receipt = json.loads(RECEIPT.read_text())
    if (len(receipt['checks']) != 21 or not all(receipt['checks'].values())
            or receipt['endpoint_sha256'] != CHECKPOINT_SHA
            or receipt['binary_sha256'] != OLD_BINARY_SHA
            or receipt['initial_sha256'] != INITIAL_SHA
            or abs(receipt['time'] - START) > 1e-14
            or digest(CHECKPOINT) != CHECKPOINT_SHA):
        raise ValueError('The accepted step-start receipt or checkpoint changed')
    template = (DOCS / 'kepler_q2_zero_pressure_v4_fresh_jitter64_to65_20261027.slurm').read_text()
    begin = template.index("/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'")
    end = template.index('\nexport GFS_GPU_VOLW2=', begin)
    restart = '''/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'
import hashlib
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

run = Path.cwd()
receipt = json.loads((run / 'parent_audit.json').read_text())
source = Path("%s")
expected = "%s"
if (len(receipt['checks']) != 21 or not all(receipt['checks'].values())
        or receipt['endpoint_sha256'] != expected
        or receipt['binary_sha256'] != "%s"
        or receipt['initial_sha256'] != "%s"
        or abs(receipt['time'] - %r) > 1e-14
        or digest(source) != expected
        or abs(snapshot_time(source, 55351) - %r) > 1e-10):
    raise ValueError('Accepted parent checkpoint failed runtime verification')
target = run / 'khout.000001.dat'
if target.exists():
    raise FileExistsError(target)
shutil.copy2(source, target)
if digest(target) != expected:
    raise ValueError('Copied checkpoint differs from accepted state')
with (run / 'restart_origin.json').open('x') as output:
    json.dump({'source_job': 437456, 'source': str(source),
               'sha256': expected, 'time': %r}, output, indent=2)
    output.write('\\n')
PY_RESTART''' % (CHECKPOINT, CHECKPOINT_SHA, OLD_BINARY_SHA, INITIAL_SHA, START, START, START)
    job = template[:begin] + restart + template[end:]
    replacements = {
        'LagQ2W64Fresh65v4': 'LagQ2W64Recover005v5',
        '--time=03:00:00': '--time=00:20:00',
        'HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.825': 'HYDRO_TSTOP=.05 EUNHA_DUMP_DT=.05',
        '--parent-time .001 --time 6.5': '--parent-time 0.045616533734626756 --time .05',
        OLD_BINARY_SHA + ' \\\n    --output audit_t65.json': BINARY_SHA + ' \\\n    --output audit_t005.json',
    }
    for old, new in replacements.items():
        if job.count(old) != 1:
            raise ValueError(f'Unexpected job template: {old!r}')
        job = job.replace(old, new)
    origins = {'eunha2.exe': BUILD / 'eunha2_capacity64_v5.exe',
               'params.dat': OLD / 'params.dat', 'initial.dat': OLD / 'initial.dat',
               'disk_active_mask.bin': OLD / 'disk_active_mask.bin', 'parent_audit.json': RECEIPT}
    for name in ('audit_kepler_q2_weighted_equilibrium_20261027.py',
                 'audit_kepler_bounded_native_evolution_20261027.py',
                 'audit_kepler_sparse_native_ballistic_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'audit_kepler_sparse_state_consistency_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        origins[name] = DOCS / name
    for name, expected in {'eunha2.exe': BINARY_SHA, 'initial.dat': INITIAL_SHA,
                           'disk_active_mask.bin': MASK_SHA}.items():
        if digest(origins[name]) != expected:
            raise ValueError(f'{name} differs from accepted input')
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    (RUN / 'job.slurm').write_text(job)
    subprocess.run(['bash', '-n', str(RUN / 'job.slurm')], check=True)
    hashes = {'job.slurm': digest(RUN / 'job.slurm')}
    for name, source in origins.items():
        shutil.copy2(source, RUN / name)
        hashes[name] = digest(RUN / name)
        if hashes[name] != digest(source):
            raise ValueError(f'Copy changed: {name}')
    with (RUN / 'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    manifest = {'source_job': 437456, 'source_checkpoint': str(CHECKPOINT),
                'checkpoint_sha256': CHECKPOINT_SHA, 'parent_time': START,
                'target_time': .05, 'binary_sha256': BINARY_SHA,
                'initial_sha256': INITIAL_SHA, 'mask_sha256': MASK_SHA,
                'jacobi_iterations': 64, 'inputs': hashes}
    with (RUN / 'stage_manifest.json').open('x') as out:
        json.dump(manifest, out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(RUN), 'start': START, 'target': .05}, indent=2))


if __name__ == '__main__':
    main()
