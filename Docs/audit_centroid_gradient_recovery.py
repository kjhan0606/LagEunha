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
out={'scope':__doc__,'common_face_area_relative_error_max':volume_error,
     'closed_polygon_area_relative_error_max':polygon_error,'cases':{}}
for name,G in [('shear',[[0,1],[0,0]]),('rotation',[[0,-1],[1,0]]),('general',[[.3,-.7],[1.2,-.3]])]:
    G=np.asarray(G,dtype=float);ubar=centroid@G.T
    # Same cell-average differences in both fits, with periodic image correction.
    du=ubar[j]-ubar[i]+(dx-(x[j]-x[i]))@G.T
    out['cases'][name]={}
    for location,z in [('mass_centroid',dc),('generator',dx)]:
        weight=area/np.maximum(np.sum(z*z,axis=1),1e-30)
        normal=acc((weight[:,None,None]*z[:,:,None]*z[:,None,:]).reshape(-1,4)).reshape(n,2,2)
        rhs=acc((weight[:,None,None]*z[:,:,None]*du[:,None,:]).reshape(-1,4)).reshape(n,2,2)
        fitted=np.linalg.solve(normal[active],rhs[active]).transpose(0,2,1)
        err=np.linalg.norm(fitted-G,axis=(1,2))
        recovered=ubar[active]-np.einsum('nij,nj->ni',fitted,offset[active])
        velocity_error=np.linalg.norm(recovered-x[active]@G.T,axis=1)
        rec={'max_gradient_error':float(err.max()),'gradient_error_volume_RMS':float(np.sqrt(np.average(err**2,weights=a[active,6]))),
             'max_generator_velocity_error':float(velocity_error.max()),
             'normal_matrix_max_condition':float(np.linalg.cond(normal[active]).max())}
        if location=='mass_centroid':assert err.max()<1e-10 and velocity_error.max()<1e-10
        out['cases'][name][location]=rec
out['limits']=['Prescribed affine averages, not canonical point momenta reinterpreted as averages.',
               'No nonuniform-density centroid, limiter, evolved gradient, shock or long-time test.']
with (root/'centroid_gradient_recovery.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
