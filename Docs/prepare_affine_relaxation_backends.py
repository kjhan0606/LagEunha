"""Actual shared-face quadrature inputs for affine relaxation backend checks."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_moving_cell_spin_20261018')
from affine_moment_patch import bounded_mesh,integrate
parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=Path('/gpfs/kjhan/LagForce/Codex_affine_relaxation_backends_20261018'))
root=parser.parse_args().root;root.mkdir(exist_ok=False)
n=8;xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
x=np.column_stack((xx.ravel(),yy.ravel()))+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
polys,edges=bounded_mesh(x,np.zeros(n*n));inputs=[];expected=[];labels=[]
for label in ['affine','curved']:
    velocity=lambda p:np.column_stack((.2+1.1*p[:,1] if label=='affine' else .5*np.sin(2*np.pi*p[:,1]),np.zeros(len(p))))
    states=[integrate(p,velocity) for p in polys]
    raw=[]
    for s in states:
        a=np.zeros(17);a[:16]=s['raw'][[0,1,2,3,4,6,7,8,9,10,11,12,13,14,16,17]];raw.append(a)
    nodes,weights=np.polynomial.legendre.leggauss(3)
    for i,j,a,b,normal,length in edges:
        if j<0:continue
        t=(b-a)/length
        for z,w in zip(nodes,weights):
            f=.5*(a+b)+.5*z*(b-a);gamma=.5*length*w
            pair=[states[i],states[j]];v=[];G=[];r=[]
            for s in pair:
                g=np.linalg.solve(s['S'],s['C']).T;arm=f-s['center']
                G.append(g);r.append(arm);v.append(s['mean']+g@arm)
            slip=float((v[0]-v[1])@t);F=-gamma*slip*t;heat=gamma*slip**2
            ref=np.zeros(35)
            for side,s in enumerate(pair):
                force=F if side==0 else -F
                dg=np.outer(force,np.linalg.solve(s['S'],r[side]))
                pdot=force
                ddot=np.outer(s['center'],force)+s['S']@dg.T
                pidot=np.outer(s['mean'],force)+np.outer(force,s['mean'])+dg@s['S']@G[side].T+G[side]@s['S']@dg.T
                d=ref[1+17*side:18+17*side];d[6:8]=pdot;d[8:12]=ddot.ravel()
                d[12:15]=[pidot[0,0],pidot[0,1],pidot[1,1]];d[15]=.5*np.trace(pidot)+.5*heat
            inputs.append(np.r_[raw[i],raw[j],f,t,gamma]);expected.append(ref);labels.append(label)
valid_template=inputs[-1].copy()
for kind in ['mass','inertia','tangent','nan']:
    a=valid_template.copy();code=1
    if kind=='mass':a[0]=0
    elif kind=='inertia':a[3]=a[5]=-1;code=2
    elif kind=='tangent':a[36:38]*=2
    else:a[7]=np.nan
    ref=np.full(35,12345.);ref[0]=code
    inputs.append(a);expected.append(ref);labels.append(kind)
np.asarray(inputs,dtype='<f8').tofile(root/'inputs.bin')
np.asarray(expected,dtype='<f8').tofile(root/'expected.bin')
(root/'manifest.json').write_text(json.dumps(dict(cases=len(inputs),labels=labels,scope=__doc__),indent=2))
print(len(inputs))
