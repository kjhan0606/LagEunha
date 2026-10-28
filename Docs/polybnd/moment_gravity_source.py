"""Acceleration source for raw18 moments, using positive mass quadrature.
Position, mass and potential energy have no acceleration-only source. Potential
work cancellation belongs to the continuity/transport equation, not this kick.
"""
import numpy as np

def source(points, mass_weights, mean_velocity, acceleration):
    x=np.asarray(points);dm=np.asarray(mass_weights);v=np.asarray(mean_velocity);a=np.asarray(acceleration)
    if x.ndim!=2 or x.shape[1]!=2 or dm.shape!=(len(x),) or v.shape!=x.shape or a.shape!=x.shape:
        raise ValueError('invalid gravity quadrature shape')
    if not all(np.all(np.isfinite(z)) for z in (x,dm,v,a)) or np.any(dm<=0):
        raise ValueError('gravity quadrature must be finite and positive')
    force=dm@a;D=x.T@(dm[:,None]*a);VA=v.T@(dm[:,None]*a)
    result=np.zeros(18)
    result[7:9]=force;result[9:13]=D.ravel();result[13:17]=(VA+VA.T).ravel()
    result[17]=np.sum(dm*np.sum(v*a,axis=1))
    return result


def covariance_source(raw, gravity_rate):
    """Fixed-position acceleration contribution to raw residual covariance.

    Includes resolved velocity curvature; zero for an exactly affine mean.
    """
    r=np.asarray(raw,dtype=float);d=np.asarray(gravity_rate,dtype=float)
    if r.shape!=(18,) or d.shape!=(18,) or not np.isfinite(r).all() or not np.isfinite(d).all():
        raise ValueError('finite raw18 state and source required')
    if r[0]<=0 or np.any(d[:7]!=0):
        raise ValueError('positive mass and acceleration-only source required')
    m=r[0];q=r[1:3];p=r[7:9]
    S=r[3:7].reshape(2,2)-np.outer(q,q)/m
    C=r[9:13].reshape(2,2)-np.outer(q,p)/m
    B=np.linalg.solve(S,C).T
    dC=d[9:13].reshape(2,2)-np.outer(q,d[7:9])/m
    dT=d[13:17].reshape(2,2)-(np.outer(d[7:9],p)+np.outer(p,d[7:9]))/m
    rate=dT-B@dC-dC.T@B.T
    return (rate+rate.T)/2


def covariance_generator(factor, rate):
    """Symmetric A with A R+R A.T=rate, R=L L.T strictly positive.

    For adding gravity to the existing exponential factor update. This is
    an instantaneous identity, not an assertion of finite-step energy accuracy.
    Rank-deficient covariance requires a separate closure and is rejected.
    """
    L=np.asarray(factor,dtype=float);G=np.asarray(rate,dtype=float)
    if L.shape!=(2,2) or G.shape!=(2,2) or not np.isfinite(L).all() or not np.isfinite(G).all():
        raise ValueError('finite 2x2 factor and symmetric source required')
    if not np.allclose(G,G.T,rtol=1e-13,atol=0):
        raise ValueError('symmetric covariance source required')
    U,s,_=np.linalg.svd(L)
    if np.any(s<=0):raise ValueError('strictly positive covariance required')
    denominator=s[:,None]**2+s[None,:]**2
    if np.any(denominator==0):raise ValueError('covariance scale underflow')
    A=U@((U.T@G@U)/denominator)@U.T
    if not np.isfinite(A).all():raise ValueError('nonfinite covariance generator')
    return (A+A.T)/2
