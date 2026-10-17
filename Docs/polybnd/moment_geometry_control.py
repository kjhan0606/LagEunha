"""Kinematic RK4 control: exact affine generator paths, weight rate enforces A'=0.
No fluid moment evolution or pressure feedback. Runs on compute node.
"""
import json
import numpy as np
from moment_ale_evolve import initial,geometry
from periodic2d_hllc import quadrature


def run(dt,end=.02):
    mass,flat=initial(4);x0=flat[:-4].reshape(16,20)[:,:2]
    velocity=np.column_stack([.2+1.1*x0[:,1],np.zeros(16)]);weights=np.zeros(16)
    def rate(t,w):
        x=x0+t*velocity;_,edges=geometry(x,w,t);base_rate=np.zeros(16);L=np.zeros((16,16))
        for i,j,a,b,n,length in edges:
            mid=(a+b)/2
            if j>=0:
                distance=np.linalg.norm(x[j]-x[i]);base=((x[j]-x[i])@velocity[j]-(velocity[j]-velocity[i])@(mid-x[i]))/distance
                c=length/(2*distance);L[i,i]+=c;L[j,j]+=c;L[i,j]-=c;L[j,i]-=c
                base_rate[i]+=length*base;base_rate[j]-=length*base
            else:base_rate[i]+=length*np.dot([.2+1.1*mid[1],0.],n)
        z=np.r_[np.linalg.solve(L[:-1,:-1],-base_rate[:-1]),0.]
        assert np.max(abs(base_rate+L@z))<1e-12
        return z
    for step in range(round(end/dt)):
        t=step*dt;k1=rate(t,weights);k2=rate(t+dt/2,weights+dt*k1/2)
        k3=rate(t+dt/2,weights+dt*k2/2);k4=rate(t+dt,weights+dt*k3)
        weights+=dt*(k1+2*k2+2*k3+k4)/6
    polygons,_=geometry(x0+end*velocity,weights,end);area=np.array([quadrature(p)[1].sum() for p in polygons])
    return dict(dt=dt,end=end,max_relative_area_drift=float(np.max(abs(area/mass-1))),max_absolute_area_drift=float(np.max(abs(area-mass))))


if __name__=='__main__':
    for dt in (.005,.0025,.00125,.000625,.0003125):print(json.dumps(run(dt)),flush=True)
