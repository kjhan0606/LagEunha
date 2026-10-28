"""Record the requested native flags and allocation; never evolves state."""
import json
import os
from pathlib import Path

KEYS=['SLURM_JOB_ID','SLURM_NODELIST','CUDA_VISIBLE_DEVICES',
      'GFS_GPU_VERTEX_HALO_CERT','GFS_KEPLER_HALO_CELLS',
      'GFS_REUSE_TESS_PADDING','GFS_GPU_LOCAL_PERIODIC_SYNC']

if __name__=='__main__':
    with Path('launch_context.json').open('x') as stream:
        json.dump({key:os.environ.get(key,'') for key in KEYS},stream,indent=2)
        stream.write('\n')
