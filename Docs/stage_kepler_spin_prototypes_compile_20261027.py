"""Compile host and GPU spin prototypes without waiting for a GPU slot."""

from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/spin_prototypes_compile_v1')


def main():
    sources={'job.slurm':DOCS/'kepler_spin_prototypes_compile_20261027.slurm',
             'host_reconstruct.cpp':DOCS/'test_gfs_affine_spin_reconstruction_20261027.cpp',
             'host_exchange.cpp':DOCS/'test_gfs_spin_affine_face_exchange_20261027.cpp',
             'gpu_reconstruct.cu':DOCS/'test_gfs_affine_spin_reconstruction_cuda_20261027.cu',
             'gpu_exchange.cu':DOCS/'test_gfs_spin_affine_face_exchange_cuda_20261027.cu',
             'gfs_affine_spin_reconstruction_20261027.h':DOCS/'gfs_affine_spin_reconstruction_20261027.h',
             'gfs_spin_affine_face_exchange_20261027.h':DOCS/'gfs_spin_affine_face_exchange_20261027.h',
             'gfs_spin_face_exchange_20261027.h':DOCS/'gfs_spin_face_exchange_20261027.h'}
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
