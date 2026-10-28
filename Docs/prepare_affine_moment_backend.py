"""Reuse saved face inputs; independent 3-point quadrature reference."""
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_affine_moment_backend_20261018');root.mkdir(exist_ok=False)
old=np.fromfile(root.parent/'Codex_affine_flux_backend_20261018/inputs.bin',dtype='<f8').reshape(-1,11)
a=np.column_stack((old,np.full(len(old),.001),np.full(len(old),.0002)))
a.astype('<f8').tofile(root/'inputs.bin');ref=np.zeros((len(a),17))
nodes,weights=np.polynomial.legendre.leggauss(3)
for node,weight in zip(nodes,weights):
    z=.5*a[:,0]*node;dm=.5*a[:,0]*weight*(a[:,1]+a[:,2]*z)
    x=a[:,3]+a[:,5]*z;y=a[:,4]+a[:,6]*z
    vx=a[:,7]+a[:,9]*z;vy=a[:,8]+a[:,10]*z;eth=a[:,11]+a[:,12]*z
    fields=np.column_stack((np.ones(len(a)),x,y,x*x,x*y,y*y,vx,vy,x*vx,x*vy,y*vx,y*vy,
                           vx*vx,vx*vy,vy*vy,eth+.5*(vx*vx+vy*vy),np.zeros(len(a))))
    ref+=dm[:,None]*fields
ref.astype('<f8').tofile(root/'expected.bin')
(root/'manifest.json').write_text(json.dumps({'cases':len(a),'reference':'Independent 3-point Gauss','scope':'17 transport moment rates only; thermal trace prescribed, potential flux zero.'},indent=2)+'\n')
