"""Static oracle: prescribed exact shear traces, not a reconstruction/evolver.
Keep the same generator velocities; solve mass constraint with exact rho=1.
"""
import json
import numpy as np
from periodic2d_shear import initial,P0
from periodic2d_hllc import mesh


def check(n=8,amplitude=.5):
    mass,y,_=initial(n,amplitude,'resolved');x=y[:,:2];v=y[:,3:5]/mass[:,None]
    _,area,edges=mesh(x,y[:,2],audit=True);count=len(x);lap=np.zeros((count,count));r=np.zeros(count);data=[]
    gx,gw=np.polynomial.legendre.leggauss(12);gx=(gx+1)/2;gw=gw/2
    for i,j,a,b,delta,length in edges:
        distance=np.linalg.norm(delta);normal=delta/distance;points=a[None,:]+gx[:,None]*(b-a)
        base=(delta@v[j]-points@(v[j]-v[i]))/distance
        actual_y=points[:,1]+(x[i]%1)[1]
        vn=amplitude*np.sin(2*np.pi*actual_y)*normal[0]
        flux=length*(gw@(vn-base));r[i]+=flux;r[j]-=flux;c=length/(2*distance)
        lap[i,i]+=c;lap[j,j]+=c;lap[i,j]-=c;lap[j,i]-=c
        data.append((i,j,length,distance,normal,vn,base))
    z=np.r_[np.linalg.solve(lap[:-1,:-1],r[:-1]),0.];netmass=np.zeros(count);transverse=np.zeros(count)
    for i,j,length,distance,normal,vn,base in data:
        flux=length*(gw@(vn-base-(z[i]-z[j])/(2*distance)))
        netmass[i]-=flux;netmass[j]+=flux
        fy=length*P0*normal[1];transverse[i]-=fy;transverse[j]+=fy
    return dict(n=n,amplitude=amplitude,mass_residual=float(np.max(abs(netmass))),
        transverse_acceleration_rms=float(np.sqrt(np.sum(transverse**2/mass))),
        generator_velocity_change=0.,weight_rate_max=float(np.max(abs(z))))


if __name__=='__main__':print(json.dumps(check(),indent=2))
