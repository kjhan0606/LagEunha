"""Uniform affine-field checks on saved t0 common faces; no evolution."""
import json
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000')
i,j=f[:,:2].astype(int).T;dx=f[:,2:4];area=f[:,4]
d=np.linalg.norm(dx,axis=1);e=dx/d[:,None];off=f[:,5:7]-.5*dx
x=a[:,1:3]-3;r=np.linalg.norm(x,axis=1);n=len(a)
masks={'inner_background':r<.5,'inner_disk_edge':(r>=.5)&(r<.6),'disk_core':(r>=.6)&(r<1.9)}
out={'scope':__doc__,'boundary_limit':'Only r<1.9 assessed, no global boundary or evolution claim.',
     'cases':{}}
for name,G in [('rotation',[[0,-1],[1,0]]),('dilation',[[1,0],[0,1]]),
               ('simple_shear',[[0,1],[0,0]]),('diagonal_strain',[[1,0],[0,-1]])]:
    G=np.array(G,dtype=float);vi=x[i]@G.T
    # Use neighbor image displacement, not wrapped absolute coordinates.
    vj=vi+dx@G.T
    face=area*(.5*np.sum((vi+vj)*e,axis=1)+np.sum((vi-vj)*off,axis=1)/d)
    rate=np.bincount(i,weights=face,minlength=n)/a[:,6];error=rate-np.trace(G)
    out['cases'][name]={'gradient':G.tolist(),'expected_material_volume_rate':float(np.trace(G)),'regions':{}}
    for label,mask in masks.items():
        out['cases'][name]['regions'][label]={'error_mass_RMS':float(np.sqrt(np.average(error[mask]**2,weights=a[mask,7]))),
          'error_max_abs':float(abs(error[mask]).max()),'error_mass_mean':float(np.average(error[mask],weights=a[mask,7]))}
    if name in ('rotation','dilation'):assert abs(error[r<1.9]).max()<1e-7,(name,abs(error[r<1.9]).max())
with (root/'global_affine_volume.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
