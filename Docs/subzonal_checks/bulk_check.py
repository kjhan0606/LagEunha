import numpy as np
from sub2d import Mesh, cross
def vols(x):
    M=Mesh(x,1.0); out=[]
    for c in M.cells:
        a=np.array([v for v,_,_ in c]); b=np.roll(a,-1,0); out.append(0.5*cross(a,b).sum())
    return np.array(out)
rng=np.random.default_rng(4); n=12; h=1/n
i,j=np.meshgrid(np.arange(n),np.arange(n),indexing='ij'); lat=np.stack([(i.ravel()+.5)*h,(j.ravel()+.5)*h],1)
for lab,x in (('square lattice',lat),('glass +-0.25h',(lat+rng.uniform(-.25,.25,lat.shape)*h)%1.0)):
    V=vols(x); eps=1e-6
    # uniform shear v_x = (y-0.5), periodic-consistent for the interior; measure interior cells
    v=np.stack([x[:,1]-0.5,np.zeros(len(x))],1)
    Vp=vols((x+eps*v)%1.0); Vm=vols((x-eps*v)%1.0); D=(Vp-Vm)/(2*eps)/V
    m=np.max(np.abs(x-0.5),1)<0.3
    print('%-16s uniform shear: rms V_dot/V (interior) = %.3e  max = %.3e   (div v = 0)'%(lab,np.sqrt(np.mean(D[m]**2)),np.abs(D[m]).max()))
