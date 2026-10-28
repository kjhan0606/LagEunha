"""Stage v2 CUDA test including both terms in moving power-face speed."""
import stage_kepler_weight_motion_remap_cuda_20261027 as stage

stage.RUN=stage.RUN.parent/'weight_motion_remap_gpu_v2'

if __name__=='__main__':stage.main()
