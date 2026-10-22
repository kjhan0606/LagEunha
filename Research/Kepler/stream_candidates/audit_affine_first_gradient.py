"""Evaluate the native LSF and TVD formulas on recorded initial native faces."""
from pathlib import Path
import json
import numpy as np
from audit_snapshot import read_snapshot

root=Path(__file__).resolve().parent
d=root/'irregular_affine_probe/rigid_irregular_lsf'
p,_=read_snapshot(d/'khout.000001.dat')
f=np.loadtxt(d/'face_dump_r0_call1.tsv',skiprows=1)
ids=(p['u4if']&np.uint64((1<<56)-1)).astype(int)
lookup=np.empty(len(ids),dtype=int);lookup[ids]=np.arange(len(ids))
i,j=lookup[f[:,4].astype(int)],lookup[f[:,5].astype(int)]
dx=p['x'][j]-p['x'][i];dy=p['y'][j]-p['y'][i]
dx-=24*np.rint(dx/24);dy-=24*np.rint(dy/24)
cf=.5*(f[:,12:14]+f[:,14:16])-np.column_stack((p['x'][i],p['y'][i]))
cf-=24*np.rint(cf/24)
area=np.linalg.norm(f[:,10:12],axis=1)
w=area/(dx*dx+dy*dy)
add=lambda a:np.bincount(i,weights=a,minlength=len(p))
Mxx,Mxy,Myy=add(w*dx*dx),add(w*dx*dy),add(w*dy*dy)
det=Mxx*Myy-Mxy*Mxy
grad=np.empty((len(p),2,2));alpha=np.ones((len(p),2))
for comp,name in enumerate(('vx','vy')):
    v=p[name]
    dv=v[j]-v[i]
    bx,by=add(w*dx*dv),add(w*dy*dv)
    grad[:,comp,0]=(Myy*bx-Mxy*by)/det
    grad[:,comp,1]=(-Mxy*bx+Mxx*by)/det
    vmin,vmax=v.copy(),v.copy()
    np.minimum.at(vmin,i,v[j]);np.maximum.at(vmax,i,v[j])
    delta=grad[i,comp,0]*cf[:,0]+grad[i,comp,1]*cf[:,1]
    ratio=np.ones(len(i));pos=delta>1e-30;neg=delta< -1e-30
    ratio[pos]=np.maximum(0,(vmax[i[pos]]-v[i[pos]])/delta[pos])
    ratio[neg]=np.maximum(0,(vmin[i[neg]]-v[i[neg]])/delta[neg])
    np.minimum.at(alpha[:,comp],i,ratio)
limited=grad*alpha[:,:,None]
B=np.array([[0,-.1],[.1,0]])
normal=f[:,10:12]/area[:,None]
vel=np.column_stack((p['vx'],p['vy']))
relj=cf-np.column_stack((dx,dy))
def jumps(g):
    left=vel[i]+np.einsum('ijk,ik->ij',g[i],cf)
    right=vel[j]+np.einsum('ijk,ik->ij',g[j],relj)
    return np.sum((right-left)*normal,axis=1)
raw_jump,limited_jump=jumps(grad),jumps(limited)
r=np.hypot(p['x']-12,p['y']-12)
out={'scope':'Independent evaluation of existing native formula on first native face geometry; no temporal evolution. Interior only avoids nonperiodic affine boundary seam.', 'radial':[]}
for lo,hi in ((0,.5),(.5,1),(1,2),(2,4),(4,8)):
    sel=(r>=lo)&(r<hi);fs=sel[i]
    out['radial'].append({'Rlo':lo,'Rhi':hi,'n':int(sel.sum()),
                         'raw_gradient_max_relative_error':float(np.max(np.linalg.norm(grad[sel]-B,axis=(1,2))/np.linalg.norm(B))),
                         'limited_gradient_max_relative_error':float(np.max(np.linalg.norm(limited[sel]-B,axis=(1,2))/np.linalg.norm(B))),
                         'cells_clipped':int(np.sum(np.any(alpha[sel]<1-1e-10,axis=1))),
                         'alpha_min':float(alpha[sel].min()),
                         'normal_jump_raw_max':float(abs(raw_jump[fs]).max()),
                         'normal_jump_limited_max':float(abs(limited_jump[fs]).max()),
                         'limited_false_strain_max':float(np.max(np.linalg.norm(.5*(limited[sel]+limited[sel].transpose(0,2,1)),axis=(1,2))))})
(root/'affine_limiter_intervention/first_gradient_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
