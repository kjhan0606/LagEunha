"""Stage corrected GPU-test timing for the two-term power-face speed."""
import stage_kepler_weight_motion_remap_cuda_20261027 as stage

stage.RUN=stage.RUN.parent/'weight_motion_remap_gpu_v3'

if __name__=='__main__':stage.main()
