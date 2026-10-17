"""Mass-only evaluation of prescribed paths for the affine 3D kinetic model.

Same positive velocity samples and upwind split as kinetic_moment_rhs_3d.
Pressure traction is absent from the mass equation; no constraint is solved.
"""
import numpy as np
from moment_state_3d import decode,ALPHA


def evaluate(geometry,sites,weights,raw,factors,geometry_rates,bounds=None,
             boundary_velocity=None,face_order=3):
    x=np.asarray(sites);raw=np.asarray(raw);mass=raw[:,0];N=len(mass)
    gen,z=map(np.asarray,geometry_rates);factors=np.asarray(factors)
    if gen.shape!=(N,3) or z.shape!=(N,) or not np.isfinite(gen).all() or not np.isfinite(z).all():
        raise ValueError('invalid prescribed geometry rates')
    if factors.shape!=(N,3,3):raise ValueError('invalid primary factors')
    cells,faces,_=geometry.build(x,weights,bounds,face_order)
    states=[decode(p,r) for p,r in zip(cells,raw)]
    for s,L,m in zip(states,factors,mass):
        if not np.allclose(s['R'],L@L.T,rtol=1e-9,atol=2e-12*max(m,1.)):
            raise ValueError('primary factor mismatch')
    boundary_velocity=(lambda p:np.zeros(3)) if boundary_velocity is None else boundary_velocity
    flux=np.zeros(N)
    for face in faces:
        i,j,n=face['i'],face['j'],face['normal'];p=face['points'];w=face['weights']
        left=states[i];vl=left['v']+(p-left['c'])@left['B'].T
        sl=vl[:,None,:]+(ALPHA@factors[i].T/np.sqrt(mass[i]))[None,:,:]
        rl=left['density'].value(p)
        if j>=0:
            right=states[j];vr=right['v']+(p-right['c'])@right['B'].T
            sr=vr[:,None,:]+(ALPHA@factors[j].T/np.sqrt(mass[j]))[None,:,:]
            rr=right['density'].value(p);delta=x[j]-x[i]
            wn=(delta@gen[j]-(p-x[i])@(gen[j]-gen[i])+.5*(z[i]-z[j]))/np.linalg.norm(delta)
        else:
            wn=np.array([boundary_velocity(point)@n for point in p]);rr=rl
            sr=sl-2*(sl@n-wn[:,None])[:,:,None]*n
        a=rl/6*np.maximum(sl@n-wn[:,None],0.).sum(axis=1)
        b=rr/6*np.minimum(sr@n-wn[:,None],0.).sum(axis=1)
        value=float(w@(a+b));flux[i]-=value
        if j>=0:flux[j]+=value
    material=np.array([s['v']+s['B']@(p-s['c']) for s,p in zip(states,x)])
    return flux,dict(generator_velocity_defect=float(np.max(abs(gen-material))),
                    mass_root_solved=False,face_order=face_order)
