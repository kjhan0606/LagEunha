"""Compare v5 and isolated cache candidate at a common accepted checkpoint."""

import json
from pathlib import Path
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import (
    BINARY_SHA, CHECKPOINT, CHECKPOINT_SHA, DOCS, INITIAL_SHA, MASK_SHA,
    RECEIPT, ROOT, RUN as RECOVERY, START, digest,
)

BUILD = ROOT/'gpu_energy_neutral_accel_cache_capacity64_build_v7'
RUN = ROOT/'gpu_energy_neutral_cache_capacity64_pair_t006_v1'
BINARY_NAME = 'eunha2_cache_capacity64_v7.exe'


def main():
    parent = json.loads(RECEIPT.read_text())
    if (len(parent['checks']) != 21 or not all(parent['checks'].values())
            or parent['endpoint_sha256'] != CHECKPOINT_SHA
            or parent['initial_sha256'] != INITIAL_SHA
            or digest(CHECKPOINT) != CHECKPOINT_SHA):
        raise ValueError('Accepted parent failed')
    build_sha = (BUILD/'binary.sha256').read_text().split()[0]
    binary = BUILD/BINARY_NAME
    if digest(binary) != build_sha or build_sha == BINARY_SHA:
        raise ValueError('Candidate binary missing, changed or identical to v5')
    manifest = json.loads((RECOVERY/'stage_manifest.json').read_text())
    if (manifest['binary_sha256'] != BINARY_SHA or manifest['initial_sha256'] != INITIAL_SHA
            or manifest['mask_sha256'] != MASK_SHA):
        raise ValueError('v5 manifest changed')
    recovery_job = (RECOVERY/'job.slurm').read_text()
    if digest(RECOVERY/'job.slurm') != manifest['inputs']['job.slurm']:
        raise ValueError('v5 job changed')
    header = recovery_job[:recovery_job.index('sha256sum -c inputs.sha256')]
    header = header.replace('LagQ2W64Recover005v5', 'LagQ2W64CachePair')
    env = recovery_job[recovery_job.index('export GFS_GPU_VOLW2='):recovery_job.index('nvidia-smi -L')]
    env = env.replace('HYDRO_TSTOP=.05 EUNHA_DUMP_DT=.05',
                      'HYDRO_TSTOP=.06 EUNHA_DUMP_DT=.06')
    # The mask path must follow the per-variant working directory.
    job = header + env + '''root=$PWD
nvidia-smi -L
for variant in control candidate; do
    cd "$root/$variant"
    sha256sum -c inputs.sha256
    export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
    if [[ $variant == control ]]; then
        expected_binary=''' + BINARY_SHA + '''
    else
        expected_binary=''' + build_sha + '''
    fi
    TIMEFORMAT='[RUN_TIME] wall=%R'
    time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
    /home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \\
        initial.dat khout.000001.dat . --parent-time ''' + repr(START) + ''' --time .06 \\
        --expected-binary-sha256 "$expected_binary" --output audit_t006.json > audit.log 2>&1
done
cd "$root"
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_fast_pair_20261027.py \\
    control/audit_t006.json candidate/audit_t006.json \\
    --control-binary-sha256 ''' + BINARY_SHA + ''' --candidate-binary-sha256 ''' + build_sha + ''' \\
    --parent-time ''' + repr(START) + ''' --time .06 --output pair_audit.json > pair_audit.log
'''
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    (RUN/'job.slurm').write_text(job)
    subprocess.run(['bash', '-n', str(RUN/'job.slurm')], check=True)
    root_sources = {'audit_kepler_weighted_fast_pair_20261027.py':
                    DOCS/'audit_kepler_weighted_fast_pair_20261027.py',
                    'audit_kepler_sparse_area_t0825_20261027.py':
                    RECOVERY/'audit_kepler_sparse_area_t0825_20261027.py',
                    'prepare_kepler_sparse_area_restart_20261027.py':
                    RECOVERY/'prepare_kepler_sparse_area_restart_20261027.py'}
    root_hashes = {'job.slurm': digest(RUN/'job.slurm')}
    for name, source in root_sources.items():
        shutil.copy2(source, RUN/name)
        root_hashes[name] = digest(RUN/name)
        if root_hashes[name] != digest(source):
            raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(root_hashes.items()):
            out.write(f'{value}  {name}\n')
    names = [name for name in manifest['inputs'] if name not in
             ('job.slurm', 'parent_audit.json', 'eunha2.exe')]
    for variant, exe in (('control', RECOVERY/'eunha2.exe'), ('candidate', binary)):
        directory = RUN/variant
        directory.mkdir()
        sources = {name: RECOVERY/name for name in names}
        sources.update({'eunha2.exe': exe, 'khout.000001.dat': CHECKPOINT})
        hashes = {}
        for name, source in sources.items():
            shutil.copy2(source, directory/name)
            hashes[name] = digest(directory/name)
            if hashes[name] != digest(source):
                raise ValueError(f'Copy changed: {variant}/{name}')
        if (hashes['khout.000001.dat'] != CHECKPOINT_SHA or
                hashes['initial.dat'] != INITIAL_SHA or hashes['disk_active_mask.bin'] != MASK_SHA):
            raise ValueError('Variant lineage changed')
        with (directory/'inputs.sha256').open('x') as out:
            for name, value in sorted(hashes.items()):
                out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent_receipt': str(RECEIPT), 'parent_sha256': CHECKPOINT_SHA,
                   'parent_time': START, 'target_time': .06,
                   'control_binary_sha256': BINARY_SHA,
                   'candidate_binary_sha256': build_sha}, out, indent=2)
        out.write('\n')
    print(RUN)


if __name__ == '__main__':
    main()
