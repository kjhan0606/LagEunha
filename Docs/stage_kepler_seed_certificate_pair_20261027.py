"""Stage same-GPU W0 and weighted certification and fallback replays."""

import json
from pathlib import Path
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import (
    BINARY_SHA, DOCS, ROOT, RUN as RECOVERY, digest,
)

BUILD = ROOT/'gpu_energy_neutral_seed_certificate_capacity64_build_v9'
RUN = ROOT/'gpu_energy_neutral_seedcert_w0_weighted_pair_v1'
W0_RECEIPT = DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json'
W_RECEIPT = ROOT/'gpu_energy_neutral_capacity64_v5_recover_to0825_a10_v1/audit_t0825.json'


def main():
    v5 = json.loads((RECOVERY/'stage_manifest.json').read_text())
    parent = {'w0': json.loads(W0_RECEIPT.read_text()),
              'weighted': json.loads(W_RECEIPT.read_text())}
    for label, receipt in parent.items():
        if not all(receipt['checks'].values()):
            raise ValueError(f'{label}: parent physics gate failed')
        for name in ('initial', 'endpoint'):
            if digest(Path(receipt[name])) != receipt[name+'_sha256']:
                raise ValueError(f'{label}: parent {name} changed')
    if abs(parent['w0']['time']-6.5)>1e-12 or abs(parent['weighted']['time']-.825)>1e-12:
        raise ValueError('Unexpected parent times')
    candidate = BUILD/'eunha2_seedcert_capacity64_v9.exe'
    candidate_sha = (BUILD/'binary.sha256').read_text().split()[0]
    if digest(candidate) != candidate_sha or candidate_sha == BINARY_SHA:
        raise ValueError('Candidate build missing or changed')
    template = (RECOVERY/'job.slurm').read_text()
    if digest(RECOVERY/'job.slurm') != v5['inputs']['job.slurm']:
        raise ValueError('Recovery runner changed')
    header = template[:template.index('sha256sum -c inputs.sha256')]
    header = header.replace('LagQ2W64Recover005v5', 'LagQ2SeedCertPair')
    env = template[template.index('export GFS_GPU_VOLW2='):template.index('nvidia-smi -L')]
    env = env.replace('HYDRO_TSTOP=.05 EUNHA_DUMP_DT=.05',
                      'HYDRO_TSTOP=.826 EUNHA_DUMP_DT=.826')
    job = header + env + '''root=$PWD
nvidia-smi -L
for case_name in w0 weighted; do
  for mode in baseline enabled fallback; do
    cd "$root/$case_name/$mode"
    sha256sum -c inputs.sha256
    unset GFS_GPU_SEED_BIN_CERT GFS_GPU_SEED_BIN_CERT_DIAG GFS_GPU_SEED_BIN_CERT_TEST_FALLBACK
    if [[ $mode != baseline ]]; then
      export GFS_GPU_SEED_BIN_CERT=1 GFS_GPU_SEED_BIN_CERT_DIAG=1
      if [[ $mode == fallback ]]; then export GFS_GPU_SEED_BIN_CERT_TEST_FALLBACK=1; fi
      expected_binary=''' + candidate_sha + '''
    else
      expected_binary=''' + BINARY_SHA + '''
    fi
    if [[ $case_name == w0 ]]; then
      unset GFS_GPU_VOLW2 GFS_GPU_VOLW2_ITERS GFS_GPU_VOLW2_MASK
      export HYDRO_TSTOP=6.501 EUNHA_DUMP_DT=6.501
    else
      export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
      export HYDRO_TSTOP=.826 EUNHA_DUMP_DT=.826
    fi
    TIMEFORMAT='[RUN_TIME] wall=%R'
    time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
    if [[ $case_name == w0 ]]; then
      /home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_short_replay_20261027.py . \\
        --time 6.501 --output audit.json > audit.log 2>&1
    else
      /home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \\
        initial.dat khout.000001.dat . --parent-time .825 --time .826 \\
        --expected-binary-sha256 "$expected_binary" --output audit.json > audit.log 2>&1
    fi
  done
done
cd "$root"
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 audit_kepler_seed_certificate_pair_20261027.py . --output pair_audit.json > pair_audit.log
'''
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    (RUN/'job.slurm').write_text(job)
    subprocess.run(['bash', '-n', str(RUN/'job.slurm')], check=True)
    common = [name for name in v5['inputs'] if name not in
              ('job.slurm', 'parent_audit.json', 'eunha2.exe', 'initial.dat')]
    common += ['audit_kepler_bounded_short_replay_20261027.py']
    root_names = ('audit_kepler_seed_certificate_pair_20261027.py',
                  'audit_kepler_sparse_area_t0825_20261027.py',
                  'prepare_kepler_sparse_area_restart_20261027.py')
    root_hashes = {'job.slurm': digest(RUN/'job.slurm')}
    for name in root_names:
        shutil.copy2(DOCS/name, RUN/name)
        root_hashes[name] = digest(RUN/name)
        if root_hashes[name] != digest(DOCS/name):
            raise ValueError(f'Root input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(root_hashes.items()):
            out.write(f'{value}  {name}\n')
    for label, receipt in parent.items():
        (RUN/label).mkdir()
        for mode in ('baseline','enabled','fallback'):
            directory=RUN/label/mode
            directory.mkdir()
            sources={name: DOCS/name if name == 'audit_kepler_bounded_short_replay_20261027.py'
                     else RECOVERY/name for name in common}
            sources.update({'eunha2.exe': RECOVERY/'eunha2.exe' if mode == 'baseline' else candidate,
                            'initial.dat': Path(receipt['initial']),
                            'khout.000001.dat': Path(receipt['endpoint'])})
            hashes={}
            for name, source in sources.items():
                shutil.copy2(source,directory/name)
                hashes[name]=digest(directory/name)
                if hashes[name] != digest(source):
                    raise ValueError(f'Copy changed: {label}/{mode}/{name}')
            if (hashes['khout.000001.dat'] != receipt['endpoint_sha256'] or
                    hashes['initial.dat'] != receipt['initial_sha256'] or
                    hashes['eunha2.exe'] != (BINARY_SHA if mode == 'baseline' else candidate_sha)):
                raise ValueError(f'Variant lineage changed: {label}/{mode}')
            with (directory/'inputs.sha256').open('x') as out:
                for name, value in sorted(hashes.items()):
                    out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'w0_receipt': str(W0_RECEIPT), 'weighted_receipt': str(W_RECEIPT),
                   'baseline_binary_sha256': BINARY_SHA, 'candidate_binary_sha256': candidate_sha,
                   'parents': {k:v['endpoint_sha256'] for k,v in parent.items()}}, out,indent=2)
        out.write('\n')
    print(RUN)


if __name__ == '__main__':
    main()
