"""Check artificial face slip from a rigid-spin-only fluid reconstruction.

Exact smooth affine fields; internal faces only. Unit damping per unit face
length is a diagnostic, not a proposed physical viscosity coefficient.
"""
import json,sys,hashlib
from pathlib import Path
import numpy as np
source=Path('/gpfs/kjhan/LagForce/Codex_moving_cell_spin_20261018')
sys.path.insert(0,str(source))
from affine_moment_patch import bounded_mesh,integrate
J=np.array([[0.,-1.],[1.,0.]])
results=[]
for n in [4,8]:
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack((xx.ravel(),yy.ravel()))+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    polygons,edges=bounded_mesh(x,np.zeros(n*n))
    for label,G in [('rigid',.7*J),('shear',np.array([[0.,1.1],[0.,0.]]))]:
        velocity=lambda p:p@G.T+np.array([.2,-.1])
        mom=[integrate(p,velocity) for p in polygons]
        center=np.array([s['center'] for s in mom]);vbar=np.array([s['mean'] for s in mom])
        grad=np.array([np.linalg.solve(s['S'],s['C']).T for s in mom])
        omega=np.array([(s['C'][0,1]-s['C'][1,0])/np.trace(s['S']) for s in mom])
        heat_rigid=heat_affine=max_rigid=max_affine=0.;count=0
        nodes,weights=np.polynomial.legendre.leggauss(3)
        for i,j,a,b,normal,length in edges:
            if j<0:continue
            count+=1;tangent=(b-a)/length
            for z,w in zip(nodes,weights):
                point=.5*(a+b)+.5*z*(b-a)
                ri=point-center[i];rj=point-center[j]
                rigid=vbar[i]+omega[i]*(J@ri)-vbar[j]-omega[j]*(J@rj)
                affine=vbar[i]+grad[i]@ri-vbar[j]-grad[j]@rj
                sr=float(rigid@tangent);sa=float(affine@tangent)
                heat_rigid+=.5*length*w*sr**2;heat_affine+=.5*length*w*sa**2
                max_rigid=max(max_rigid,abs(sr));max_affine=max(max_affine,abs(sa))
        assert max_affine<1e-12
        if label=='rigid':assert max_rigid<1e-12
        else:assert heat_rigid>1e-6
        results.append(dict(n=n,field=label,internal_faces=count,rigid_trace_heat_rate=heat_rigid,
                            full_affine_trace_heat_rate=heat_affine,
                            max_rigid_trace_slip=max_rigid,max_full_affine_trace_slip=max_affine))
print(json.dumps(dict(scope=__doc__,cases=results,
    source_sha256=hashlib.sha256((source/'affine_moment_patch.py').read_bytes()).hexdigest(),
    decision='Do not connect rigid-spin-only slip damping to the Kepler solver. '
    'A full affine trace passes this null test, but force work must then include '
    'affine deformation degrees of freedom; pair torque alone does not close it.',
    limits=['No time evolution or physical viscosity calibrated.',
    'Affine null test is necessary, not sufficient for curved shear, residual transport or shocks.']),indent=2))
