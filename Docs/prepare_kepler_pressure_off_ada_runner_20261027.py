"""Make a hash-recorded single-host Ada runner from the staged A10 runner."""

import hashlib
import json
from pathlib import Path


RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/pressure_force_off_jitter_from0_to0825_a10_v1')
SOURCE=RUN/'job.slurm'
TARGET=RUN/'job_ada.sh'
RECEIPT=RUN/'ada_runner_receipt.json'
OLD='export I_MPI_HYDRA_BOOTSTRAP=slurm\n'
NEW='unset I_MPI_HYDRA_BOOTSTRAP\n'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    original=SOURCE.read_text()
    if original.count(OLD)!=1 or TARGET.exists() or RECEIPT.exists():
        raise ValueError('Unexpected runner source or prior Ada runner')
    staged=json.loads((RUN/'stage_manifest.json').read_text())
    if digest(SOURCE)!=staged['inputs']['job.slurm']['sha256']:
        raise ValueError('Staged Slurm runner changed')
    TARGET.write_text(original.replace(OLD,NEW,1))
    receipt={'slurm_runner':str(SOURCE),'slurm_runner_sha256':digest(SOURCE),
             'ada_runner':str(TARGET),'ada_runner_sha256':digest(TARGET),
             'sole_change':'Remove Slurm-only MPI bootstrap for one-rank manual Ada execution',
             'diagnostic_binary_sha256':staged['diagnostic_binary_sha256'],
             'initial_sha256':staged['initial_sha256']}
    with RECEIPT.open('x') as stream:
        json.dump(receipt,stream,indent=2)
        stream.write('\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
