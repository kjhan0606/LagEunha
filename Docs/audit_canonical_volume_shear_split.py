"""Apply saved geometric volume operator to local circular Taylor fields.

Each row uses its own Taylor expansion: these are diagnostics, not a global
velocity field or an alternative conservative evolution operator.
"""
import json
from pathlib import Path
import numpy as np
from audit_canonical_pair_history import gravity_coefficient

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000')
i,j=f[:,:2].astype(int).T;dx=f[:,2:4];area=f[:,4];dist=np.linalg.norm(dx,axis=1)
e=dx/dist[:,None];off=f[:,5:7]-.5*dx;n=len(a)
x=a[:,1:3]-3;v=a[:,3:5];r=np.linalg.norm(x,axis=1)
omega=np.sqrt(gravity_coefficient(r));eps=1e-5
domega=(np.sqrt(gravity_coefficient(r*np.exp(eps)))-np.sqrt(gravity_coefficient(r*np.exp(-eps))))/(2*eps)
jx=np.column_stack((-x[:,1],x[:,0]))
def row_rate(vi,vj):
    z=area*(.5*np.sum((vi+vj)*e,axis=1)+np.sum((vi-vj)*off,axis=1)/dist)
    return np.bincount(i,weights=z,minlength=n)/a[:,6]
actual=row_rate(v[i],v[j]);active=r<1.9
err=float(np.max(abs(actual[active]-a[active,12]/a[active,6])))
assert err<1e-9,err
translation=row_rate(v[i],v[i])
rot_delta=omega[i,None]*np.column_stack((-dx[:,1],dx[:,0]))
shear_delta=jx[i]*(domega[i]*np.sum(x[i]*dx,axis=1)/np.maximum(r[i]**2,1e-30))[:,None]
rigid=row_rate(v[i],v[i]+rot_delta)
affine=row_rate(v[i],v[i]+rot_delta+shear_delta)
curvature=actual-affine
out={'scope':__doc__,'reconstruction_max_abs_error_r_lt_1p9':err,
     'translation_rate_max_abs_r_lt_1p9':float(abs(translation[active]).max()),
     'rigid_rotation_rate_max_abs_r_lt_1p9':float(abs(rigid[active]).max()),'bins':[]}
edges=[0,.1,.2,.3,.4,.5,.6,.8,1,1.5,1.9]
for lo,hi in zip(edges[:-1],edges[1:]):
    mask=(r>=lo)&(r<hi);w=a[mask,7]
    if not mask.any():continue
    avg=lambda z:float(np.average(z[mask],weights=w))
    rms=lambda z:np.sqrt(avg(z*z))
    out['bins'].append({'r_interval':[lo,hi],'actual_RMS':rms(actual),
        'local_affine_RMS':rms(affine),'curvature_remainder_RMS':rms(curvature),
        'affine_curvature_cross_mean':avg(affine*curvature),
        'variance_identity_residual':avg(actual**2-affine**2-curvature**2-2*affine*curvature)})
with (root/'volume_shear_split.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
