"""Constitutive frame-change audit of the auxiliary vector, not hydro stability.

Rebuild the same smooth residual on a rotated cloud. Compare zdot=Bv with
the candidate zdot=Bv+Wz, W=antisym(local reconstructed velocity gradient).
Time-dependent frames also need inertial terms in the momentum equation;
this test concerns only the auxiliary constitutive equation.
"""
import json,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_memory_circular_work_20261018')
from smooth_owner_memory import build_row

rng=np.random.default_rng(20261018)
x=np.vstack((np.zeros((1,2)),rng.uniform(-.5,.5,(100,2))))
v=x@np.array([[.4,1.2],[-.7,-.1]]).T+.03*rng.normal(size=x.shape)
z=np.array([.37,-.2]);J=np.array([[0.,-1.],[1.,0.]])
angle=.71; R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
omega=2.3;Omega=omega*J

def evaluate(xx,vv,zz):
    ids,b,meta=build_row(xx+5,0,box=10)
    d=xx[ids]-xx[0]; radius=meta['radius']; q=d/radius
    weights=np.maximum(1-np.sum(q*q,axis=1),0)**4
    P=np.column_stack((np.ones(len(ids)),q,q[:,0]**2,q[:,0]*q[:,1],q[:,1]**2))
    coef=np.linalg.lstsq(np.sqrt(weights)[:,None]*P,np.sqrt(weights)[:,None]*vv[ids],rcond=None)[0]
    grad=coef[1:3].T/radius
    W=.5*(grad-grad.T)
    rate=b@vv[ids]
    force=np.zeros_like(xx);force[ids]=-b[:,None]*zz
    return {'rate':rate,'W':W,'force':force,'power':float(np.sum(force*vv)+zz@(rate+W@zz)),
            'radius':radius,'moments':meta['moment_max']}

base=evaluate(x,v,z)
xr=x@R.T; zr=R@z
constant=evaluate(xr,v@R.T,zr)
moving=evaluate(xr,v@R.T+xr@Omega.T,zr)
expected_old=R@base['rate']+Omega@zr
expected_new=R@(base['rate']+base['W']@z)+Omega@zr
old_defect=moving['rate']-expected_old
new_defect=moving['rate']+moving['W']@zr-expected_new
result={'scope':__doc__,
 'constant_rotation_rate_residual':float(np.linalg.norm(constant['rate']-R@base['rate'])),
 'force_covariance_residual':float(np.linalg.norm(moving['force']-base['force']@R.T)),
 'gradient_spin_transform_residual':float(np.linalg.norm(moving['W']-(R@base['W']@R.T+Omega))),
 'current_time_dependent_rotation_defect':float(np.linalg.norm(old_defect)),
 'predicted_missing_spin_norm':float(np.linalg.norm(Omega@zr)),
 'current_defect_plus_missing_spin_residual':float(np.linalg.norm(old_defect+Omega@zr)),
 'candidate_time_dependent_rotation_residual':float(np.linalg.norm(new_defect)),
 'candidate_power_residuals':[base['power'],constant['power'],moving['power']],
 'moment_max':max(base['moments'],constant['moments'],moving['moments']),
 'limitations':'One cloud algebra audit, scale=1. No hydro integration, no proof this causes Kepler growth, no candidate production implementation. Spin is reconstructed from actual velocities; no prescribed circular velocity.'}
assert result['current_time_dependent_rotation_defect']>.1
for key in ('constant_rotation_rate_residual','force_covariance_residual','gradient_spin_transform_residual',
            'current_defect_plus_missing_spin_residual','candidate_time_dependent_rotation_residual','moment_max'):
    assert result[key]<1e-12
assert max(abs(x) for x in result['candidate_power_residuals'])<1e-12
print(json.dumps(result,indent=2))
