"""Representation comparison at the failed exact-pressure patch face."""
import sys,json,numpy as np
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
root=Path('/gpfs/kjhan/LagForce/Codex_checkpoint_patch_rhs_exact_failure_r0.5_20261018')
f=json.loads((root/'report.json').read_text())['results'][0]['failure_frames'];a=f[0]['values'];b=f[1]['values'];n=np.array(a['n']);i,j=a['i'],a['j']
s=np.load(root/'source.npz');states=s['state'][:-4].reshape(-1,20);ids=s['ids']
ic,t=read(Path('/gpfs/kjhan/LagForce/Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'))
original=np.column_stack([ic['vx'][ids],ic['vy'][ids]]);mean=states[:,9:11]/s['mass'][:,None]
print(json.dumps(dict(cell_ids=[int(ids[k]) for k in (i,j)],cell_radii=[float(np.linalg.norm(states[k,:2])) for k in (i,j)],frozen=[bool(s['frozen'][k]) for k in (i,j)],face_radius=float(np.linalg.norm(a['point'])),original_point_velocity_jump=float((original[j]-original[i])@n),remapped_mean_velocity_jump=float((mean[j]-mean[i])@n),affine_face_velocity_jump=b['ur']-b['ul'],face_vacuum_threshold=2*(b['cl']+b['cr'])/(b['gamma']-1),scope='Original point velocities are not exact face traces. No physical vacuum conclusion for production solution.'),indent=2))
