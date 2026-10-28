"""Affine reconstruction from manufactured cell averages on actual t0 geometry.
Uses uniform density so geometric centroid equals mass centroid.
"""
import json
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000');i,j=f[:,:2].astype(int).T
x=a[:,1:3]-3;r=np.linalg.norm(x,axis=1);n=len(a);active=r<1.9
dx=f[:,2:4];area=f[:,4];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None];cf=f[:,5:7]
def acc(z):
    return np.column_stack([np.bincount(i,weights=z[:,k],minlength=n) for k in range(z.shape[1])])
volume=.5*np.bincount(i,weights=area*np.sum(cf*e,axis=1),minlength=n)
volume_error=float(np.max(abs(volume[active]/a[active,6]-1)))
# Common faces are intersected and can leave tiny closure gaps. Use the
# original closed polygon edges for the cell volume and centroid instead.
raw=np.loadtxt(root/'raw_faces.dat.rank0000');owner=raw[:,0].astype(int)
va=raw[:,6:8];vb=raw[:,8:10];cross=va[:,0]*vb[:,1]-va[:,1]*vb[:,0]
signed_volume=.5*np.bincount(owner,weights=cross,minlength=n)
polygon_error=float(np.max(abs(abs(signed_volume[active])/a[active,6]-1)))
assert polygon_error<1e-10,polygon_error
q=np.column_stack([np.bincount(owner,weights=cross*(va[:,k]+vb[:,k])/6,minlength=n) for k in range(2)])
offset=np.divide(q,signed_volume[:,None],out=np.zeros_like(q),where=signed_volume[:,None]!=0)
centroid=x+offset
dc=dx+offset[j]-offset[i]

from run_canonical_av_branches_20261018 import read
ic,_=read(root/'khout.000001.dat')
du=np.column_stack((ic['vx'][j]-ic['vx'][i],ic['vy'][j]-ic['vy'][i]))
weight=area/np.maximum(np.sum(dc*dc,axis=1),1e-30)
normal=acc((weight[:,None,None]*dc[:,:,None]*dc[:,None,:]).reshape(-1,4)).reshape(n,2,2)
rhs=acc((weight[:,None,None]*dc[:,:,None]*du[:,None,:]).reshape(-1,4)).reshape(n,2,2)
grad=np.zeros((n,2,2));grad[active]=np.linalg.solve(normal[active],rhs[active]).transpose(0,2,1)
# Positive three-node triangle integration of centered geometric covariance.
cov=np.zeros((n,2,2))
triangle_area=.5*abs(cross)
for u,v in [(1/6,1/6),(2/3,1/6),(1/6,2/3)]:
    pos=u*va+v*vb-offset[owner]
    for k in range(2):
        for l in range(2):
            cov[:,k,l]+=np.bincount(owner,weights=triangle_area/3*pos[:,k]*pos[:,l],minlength=n)
cov*=np.divide(ic['mass'],abs(signed_volume),out=np.zeros(n),where=signed_volume!=0)[:,None,None]
kin=.5*np.einsum('nij,njk,nik->n',grad,cov,grad)
out={'scope':'Actual t0 masses, velocities and U. Hypothetical replacement of per-cell constant velocity by centroid-affine field at same mean momentum; no conversion performed.',
     'cross_term_caveat':'For old constant-velocity lift C=T=0. Treating a neighbor gradient as an actual new affine field adds K_affine. Keeping it merely as a basis requires cross/residual terms; negative T-GSGT alone is not an impossibility proof.',
     'regions':{}}
for name,mask in [('inner',(r<.5)&(ic['den']<.5)),('edge',(r<.6)&(ic['den']>.5)),('core',(r>=.6)&(r<1.9)&(ic['den']>.5))]:
    ratio=kin[mask]/ic['ie'][mask]
    out['regions'][name]={'cells':int(mask.sum()),'added_affine_K_over_original_U_total':float(kin[mask].sum()/ic['ie'][mask].sum()),
       'per_cell_K_over_U_quantiles':np.quantile(ratio,[.5,.9,.99,1]).tolist(),
       'cells_K_greater_than_U':int((ratio>1).sum()),
       'mass_fraction_K_greater_than_U':float(ic['mass'][mask][ratio>1].sum()/ic['mass'][mask].sum())}
with (root/'restart_affine_energy_audit.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
