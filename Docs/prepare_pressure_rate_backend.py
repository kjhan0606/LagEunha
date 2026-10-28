"""Manufactured affine pressure and velocity on original closed t0 cells."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');src=base/'Codex_pair_hg_modes_t0_20261018'
root=base/'Codex_pressure_rate_backend_v2_20261018';root.mkdir(exist_ok=False)
a=np.loadtxt(src/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])];n=len(a);x=a[:,1:3]-3
active=np.linalg.norm(x,axis=1)<1.9
f=np.loadtxt(src/'raw_faces.dat.rank0000');f=f[active[f[:,0].astype(int)]];i=f[:,0].astype(int)
va=f[:,6:8];vb=f[:,8:10];line=vb-va;A=np.linalg.norm(line,axis=1);t=line/A[:,None]
normal=f[:,4:6]/np.linalg.norm(f[:,4:6],axis=1)[:,None];c=x[i]+.5*(va+vb)
cross=va[:,0]*vb[:,1]-va[:,1]*vb[:,0]
V=.5*np.bincount(i,weights=cross,minlength=n)
q=np.column_stack([np.bincount(i,weights=cross*(va[:,k]+vb[:,k])/6,minlength=n) for k in range(2)])+V[:,None]*x
sgn=np.sign(V);V=abs(V);q*=sgn[:,None]
polygon_normal=sgn[i,None]*np.column_stack((t[:,1],-t[:,0]))
normal_difference=float(np.max(np.linalg.norm(polygon_normal-normal,axis=1)))
normal=polygon_normal
assert np.max(abs(V[active]/a[active,6]-1))<1e-10
gp=np.array([.001,.002]);G=np.array([[.3,-.7],[1.2,-.1]])
p0=.01+c@gp;pt=t@gp;v=c@G.T;gt=t@G.T
Pint=.01*V+q@gp;nf=np.bincount(i,minlength=n)
inp=np.column_stack((A,p0,pt,c,t,normal,v,gt,Pint[i]/nf[i],np.tile(G.ravel(),(len(i),1))))
ref=np.zeros((len(i),17));nodes,weights=np.polynomial.legendre.leggauss(3)
for z0,w in zip(nodes,weights):
    z=.5*A*z0;pos=c+z[:,None]*t;vel=v+z[:,None]*gt;P=p0+pt*z;force=-(.5*A*w*P)[:,None]*normal
    ref[:,6:8]+=force
    ref[:,8:12]+=(pos[:,:,None]*force[:,None,:]).reshape(-1,4)
    ref[:,12]+=2*vel[:,0]*force[:,0];ref[:,13]+=vel[:,0]*force[:,1]+vel[:,1]*force[:,0];ref[:,14]+=2*vel[:,1]*force[:,1]
    ref[:,15]+=np.sum(vel*force,axis=1)
share=Pint[i]/nf[i];ref[:,8]+=share;ref[:,11]+=share
ref[:,12]+=2*share*G[0,0];ref[:,13]+=share*(G[0,1]+G[1,0]);ref[:,14]+=2*share*G[1,1]
expected=np.zeros((n,17));expected[:,6:8]=-V[:,None]*gp
expected[:,8:12]=(-q[:,:,None]*gp[None,None,:]).reshape(n,4)
p=q@G.T;expected[:,12]=-2*p[:,0]*gp[0];expected[:,13]=-p[:,0]*gp[1]-p[:,1]*gp[0];expected[:,14]=-2*p[:,1]*gp[1]
expected[:,15]=-p@gp-Pint*np.trace(G)
summed=np.column_stack([np.bincount(i,weights=ref[:,k],minlength=n) for k in range(17)])
err=float(np.max(abs((summed-expected)[active])/V[active,None]));assert err<1e-10,err
inp.astype('<f8').tofile(root/'inputs.bin');ref.astype('<f8').tofile(root/'expected.bin')
np.savez(root/'assembly.npz',owner=i,active=active,volume=V,expected=expected,Pint=Pint,G=G)
(root/'manifest.json').write_text(json.dumps({'cases':len(i),'cells':int(active.sum()),'generator_vs_polygon_normal_max_difference':normal_difference,'reference':'Independent3pointface quadrature plus volume terms; compared with volume pressure-gradient integral','assembled_reference_max_volume_scaled_error':err,'scope':'Smooth affine manufactured fields, closed polygon normals, no Riemann solver or time evolution.'},indent=2)+'\n')
print((root/'manifest.json').read_text())
