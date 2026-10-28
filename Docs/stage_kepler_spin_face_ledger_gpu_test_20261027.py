"""Stage the isolated device execution test of the spin face ledger."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/spin_face_ledger_gpu_test_v1')


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if RUN.exists():raise FileExistsError(RUN)
    sources={'job.slurm':DOCS/'kepler_spin_face_ledger_gpu_test_20261027.slurm',
             'gfs_spin_face_exchange_20261027.h':DOCS/'gfs_spin_face_exchange_20261027.h',
             'test_gfs_spin_face_exchange_cuda_20261027.cu':DOCS/'test_gfs_spin_face_exchange_cuda_20261027.cu'}
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    for name,source in sources.items():
        shutil.copy2(source,RUN/name)
        if digest(RUN/name)!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as f:
        for name in sorted(sources):f.write(f'{digest(RUN/name)}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump({'run':str(RUN),'inputs':{name:digest(RUN/name) for name in sources},
                   'scope':'One-thread CUDA execution of an isolated conservative spin face ledger; '
                           'not a native hydrodynamic run.'},f,indent=2);f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(sources)}))


if __name__=='__main__':main()
