"""Stage the superset device test including the excess-spin reduction."""

import stage_kepler_affine_spin_reconstruction_cuda_20261027 as stage

stage.RUN=stage.RUN.parent/'affine_spin_reconstruction_gpu_v2'

if __name__=='__main__':stage.main()
