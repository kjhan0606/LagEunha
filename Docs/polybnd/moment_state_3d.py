"""35-component raw moment layout and 3D kinetic face flux.

The layout retains full 3x3 tensors for compatibility with the 2D moment
construction: m,q,Q,p,D,Pi,E. Geometry and time integration are separate.
"""
import numpy as np
from moment_density_3d import Density3D

ALPHA=np.sqrt(3.)*np.vstack([np.eye(3),-np.eye(3)])


def pack(m,q,Q,p,D,Pi,E):
    return np.r_[m,np.asarray(q),np.asarray(Q).ravel(),np.asarray(p),
                 np.asarray(D).ravel(),np.asarray(Pi).ravel(),E]


def unpack(raw):
    r=np.asarray(raw,dtype=float)
    if r.shape!=(35,) or not np.all(np.isfinite(r)):raise ValueError('35 finite raw moments required')
    return r[0],r[1:4],r[4:13].reshape(3,3),r[13:16],r[16:25].reshape(3,3),r[25:34].reshape(3,3),r[34]


def decode(vertices,raw,gamma=5/3,order=6):
    m,q,Q,p,D,Pi,E=unpack(raw);density=Density3D(vertices,m,q,Q,order)
    c=q/m;mean=p/m;S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
    S=(S+S.T)/2;B=np.linalg.solve(S,C).T;R=Pi-np.outer(p,p)/m-B@S@B.T;R=(R+R.T)/2
    scale=max(float(np.linalg.norm(Pi)),np.finfo(float).tiny)
    if np.linalg.eigvalsh(R).min() < -1024*np.finfo(float).eps*scale:
        raise ValueError('non-realizable velocity moments; covariance not clipped')
    internal=E-.5*np.trace(Pi)
    if internal<=0 or not gamma>1:raise ValueError('positive internal energy and gamma>1 required')
    return dict(m=m,c=c,v=mean,S=S,B=B,R=R,P=(gamma-1)*internal/density.volume,density=density)


def face_flux(point,normal,left_velocity,right_velocity,left_rates,right_rates,
              left_internal_specific,right_internal_specific,star_pressure,star_speed):
    """Signed transport rates already include density and relative face speed.

    The same returned flux must be subtracted from the left cell and added to
    the right. Pressure/work is shared, with no independent side evaluation.
    """
    x=np.asarray(point);n=np.asarray(normal);vl=np.asarray(left_velocity);vr=np.asarray(right_velocity)
    a=np.asarray(left_rates);b=np.asarray(right_rates)
    if x.shape!=(3,) or n.shape!=(3,) or vl.shape!=(len(a),3) or vr.shape!=(len(b),3):
        raise ValueError('3D face samples required')
    if np.any(a<0) or np.any(b>0) or min(left_internal_specific,right_internal_specific,star_pressure)<=0:
        raise ValueError('positive upwind samples and thermodynamic states required')
    fm=a.sum()+b.sum();fp=a@vl+b@vr+star_pressure*n
    Pi=vl.T@(a[:,None]*vl)+vr.T@(b[:,None]*vr)
    energy=a@(.5*np.sum(vl*vl,axis=1)+left_internal_specific)+b@(.5*np.sum(vr*vr,axis=1)+right_internal_specific)+star_pressure*star_speed
    return pack(fm,x*fm,np.outer(x,x)*fm,fp,np.outer(x,fp),Pi,energy)
