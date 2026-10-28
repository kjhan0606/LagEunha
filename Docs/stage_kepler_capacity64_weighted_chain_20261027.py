"""Stage two afterok-gated v5 continuations from the short recovery."""

import json
from pathlib import Path
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import (
    BINARY_SHA, INITIAL_SHA, MASK_SHA, ROOT, RUN as RECOVERY, digest,
)

RUNS = [
    (ROOT / 'gpu_energy_neutral_capacity64_v5_recover_to0825_a10_v1',
     RECOVERY, 'audit_t005.json', .05, .825, 'audit_t0825.json'),
    (ROOT / 'gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1',
     ROOT / 'gpu_energy_neutral_capacity64_v5_recover_to0825_a10_v1',
     'audit_t0825.json', .825, 6.5, 'audit_t65.json'),
]
WALLTIME = '03:00:00'


def main():
    manifest = json.loads((RECOVERY / 'stage_manifest.json').read_text())
    if (manifest['binary_sha256'] != BINARY_SHA or manifest['initial_sha256'] != INITIAL_SHA
            or manifest['mask_sha256'] != MASK_SHA or manifest['target_time'] != .05):
        raise ValueError('Recovery manifest changed')
    template = (RECOVERY / 'job.slurm').read_text()
    if digest(RECOVERY / 'job.slurm') != manifest['inputs']['job.slurm']:
        raise ValueError('Recovery runner changed')
    sources = [p for p in manifest['inputs'] if p not in ('job.slurm', 'parent_audit.json')]
    for name in sources:
        if digest(RECOVERY / name) != manifest['inputs'][name]:
            raise ValueError(f'Recovery input changed: {name}')
    for index, (run, parent, receipt, start, stop, output) in enumerate(RUNS):
        if run.exists():
            raise FileExistsError(run)
        begin = template.index("/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'")
        end = template.index('\nexport GFS_GPU_VOLW2=', begin)
        restart = '''/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

run = Path.cwd()
parent = Path(%r)
receipt_path = parent / %r
receipt = json.loads(receipt_path.read_text())
if (len(receipt['checks']) != 21 or not all(receipt['checks'].values())
        or receipt['binary_sha256'] != %r or receipt['initial_sha256'] != %r
        or abs(receipt['time'] - %r) > 1e-12):
    raise ValueError('Parent weighted physics gate failed')
source = Path(receipt['endpoint'])
if (source.parent != parent or digest(source) != receipt['endpoint_sha256']
        or abs(snapshot_time(source, 55351) - %r) > 1e-10):
    raise ValueError('Parent checkpoint changed or is incomplete')
target = run / 'khout.000001.dat'
if target.exists():
    raise FileExistsError(target)
shutil.copy2(source, target)
if digest(target) != receipt['endpoint_sha256']:
    raise ValueError('Copied checkpoint differs from accepted parent')
with (run / 'restart_origin.json').open('x') as output:
    json.dump({'receipt': str(receipt_path), 'source': str(source),
               'sha256': receipt['endpoint_sha256'], 'time': %r}, output, indent=2)
    output.write('\\n')
PY_RESTART''' % (str(parent), receipt, BINARY_SHA, INITIAL_SHA, start, start, start)
        job = template[:begin] + restart + template[end:]
        replacements = {
            'LagQ2W64Recover005v5': f'LagQ2W64v5Chain{index+1}',
            '--time=00:20:00': f'--time={WALLTIME}',
            'HYDRO_TSTOP=.05 EUNHA_DUMP_DT=.05': f'HYDRO_TSTOP={stop} EUNHA_DUMP_DT=0.825',
            '--parent-time 0.045616533734626756 --time .05': f'--parent-time {start} --time {stop}',
            '--output audit_t005.json': f'--output {output}',
        }
        for old, new in replacements.items():
            if job.count(old) != 1:
                raise ValueError(f'Unexpected recovery job template: {old!r}')
            job = job.replace(old, new)
        run.mkdir()
        (run / 'job.slurm').write_text(job)
        subprocess.run(['bash', '-n', str(run / 'job.slurm')], check=True)
        hashes = {'job.slurm': digest(run / 'job.slurm')}
        for name in sources:
            shutil.copy2(RECOVERY / name, run / name)
            hashes[name] = digest(run / name)
            if hashes[name] != manifest['inputs'][name]:
                raise ValueError(f'Staged input changed: {name}')
        with (run / 'inputs.sha256').open('x') as out:
            for name, value in sorted(hashes.items()):
                out.write(f'{value}  {name}\n')
        with (run / 'stage_manifest.json').open('x') as out:
            json.dump({'parent': str(parent), 'parent_receipt': receipt,
                       'start': start, 'stop': stop, 'binary_sha256': BINARY_SHA,
                       'inputs': hashes}, out, indent=2)
            out.write('\n')
        print(f'{run}: {start} -> {stop}')


if __name__ == '__main__':
    main()
