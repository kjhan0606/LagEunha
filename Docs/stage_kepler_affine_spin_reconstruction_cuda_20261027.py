"""Stage a small device test of the affine-spin algebra."""

from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/affine_spin_reconstruction_gpu_v1')


def main():
    sources={'job.slurm':DOCS/'kepler_affine_spin_reconstruction_cuda_20261027.slurm',
             'test.cu':DOCS/'test_gfs_affine_spin_reconstruction_cuda_20261027.cu',
             'gfs_affine_spin_reconstruction_20261027.h':DOCS/'gfs_affine_spin_reconstruction_20261027.h'}
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir();hashes={}
    for name,source in sources.items():
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
