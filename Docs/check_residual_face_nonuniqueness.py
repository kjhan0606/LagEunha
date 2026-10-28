"""Construct two smooth positive-pressure fields with identical stored moments.

Their face fluxes differ: covariance alone cannot recover subcell face traces.
This is an information-loss audit, not a proof that moment closures cannot converge.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

def field(x,kind):
    # On [-1,1]^2: both modes are orthogonal to 1,x,y and have equal L2 norm.
    if kind==2:return .1*(3*x*x-1)/2
    return .1*np.sqrt(7/5)*(5*x*x*x-3*x)/2

def moments(kind):
    z,w=leggauss(12);x,y=np.meshgrid(z,z);weights=np.outer(w,w).ravel()
    points=np.column_stack([x.ravel(),y.ravel()]);v=np.column_stack([field(points[:,0],kind),np.zeros(len(points))])
    rho=1.;pressure=.002777777777777778;gamma=5/3
    m=weights.sum();q=weights@points;Q=points.T@(weights[:,None]*points)
    p=weights@v;D=points.T@(weights[:,None]*v);Pi=v.T@(weights[:,None]*v)
    energy=pressure*m/(gamma-1)+.5*np.trace(Pi)
    return np.r_[m,q,Q[0,0],Q[0,1],Q[1,1],p,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],energy,0.]

def face(kind,x):
    u=float(field(x,kind));pressure=.002777777777777778;gamma=5/3
    # Flux in +x direction per unit face length, stationary face, rho=1.
    return dict(velocity=u,mass=u,momentum_x=u*u+pressure,
        total_energy=(.5*u*u+gamma*pressure/(gamma-1))*u)

a,b=moments(2),moments(3)
assert a.shape==(17,) and np.max(abs(a-b))<1e-14
assert abs(a[6])<1e-14 and np.max(abs(a[8:12]))<1e-14
assert abs(a[12]-.008)<1e-14
faces={str(k):{str(x):face(k,x) for x in (-1.,1.)} for k in (2,3)}
assert abs(faces['2']['1.0']['momentum_x']-faces['3']['1.0']['momentum_x'])>.0039
print(json.dumps(dict(scope=__doc__,moments_P2=a.tolist(),moments_P3=b.tolist(),
    maximum_moment_difference=float(np.max(abs(a-b))),faces=faces,
    affine_reconstruction_velocity=0.,residual_covariance_per_mass_xx=.002,
    donor_predicted_right_momentum_flux=.002777777777777778,
    gaussian_predicted_right_momentum_flux=.004777777777777778,
    conclusion='Neither stored affine velocity nor residual covariance selects the true face flux. Additional spatial reconstruction/closure and its convergence evidence are needed. No production force has been changed.'),indent=2))
