"""Integrate affine upwind mass samples after splitting their zero-speed lines.

This is an independent quadrature diagnostic of the existing kinetic model.
It does not modify evolution, weights, masses, or the prescribed face path.
"""
import numpy as np
from power_cells_3d import face_quadrature
from moment_state_3d import decode,ALPHA


def clip_face(vertices,normal,level):
    p=np.asarray(vertices);center=p.mean(axis=0);n=np.asarray(normal)
    axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);v=np.cross(n,u)
    p=p[np.argsort(np.arctan2((p-center)@v,(p-center)@u))]
    values=np.asarray(level(p));result=[]
    for a,b,fa,fb in zip(p,np.roll(p,-1,axis=0),values,np.roll(values,-1)):
        if fa>=0:result.append(a)
        if (fa<0<fb) or (fb<0<fa):result.append(a+(b-a)*(fa/(fa-fb)))
    unique=[]
    for q in result:
        if not unique or not np.array_equal(q,unique[-1]):unique.append(q)
    if len(unique)>1 and np.array_equal(unique[0],unique[-1]):unique.pop()
    if len(unique)<3:return np.empty((0,3))
    return np.asarray(unique)


def integrate_sample(vertices,normal,speed,density,sign,order,return_capacity=False):
    if sign not in (-1,1):raise ValueError('sign must be -1 or +1')
    p=clip_face(vertices,normal,lambda x:sign*speed(x))
    if len(p)<3:return (0.,0.) if return_capacity else 0.
    area=np.linalg.norm(np.cross(p-p.mean(0),np.roll(p,-1,axis=0)-p.mean(0)),axis=1).sum()/2
    if area==0:return (0.,0.) if return_capacity else 0.
    x,w=face_quadrature(p,normal,order)
    rho=density(x);value=float(w@(rho*sign*np.maximum(sign*speed(x),0.)))
    if return_capacity:
        capacity=float(w@rho)
        if np.all(speed(np.asarray(vertices))==0):capacity*=.5
        return value,capacity
    return value


def evaluate(geometry,sites,weights,raw,factors,geometry_rates,order=4,return_jacobian=False,bounds=None,return_diagnostics=False):
    x=np.asarray(sites);raw=np.asarray(raw);mass=raw[:,0];gen,z=map(np.asarray,geometry_rates)
    if gen.shape!=x.shape or z.shape!=mass.shape or not np.isfinite(gen).all() or not np.isfinite(z).all():
        raise ValueError('invalid prescribed geometry rates')
    cells,faces,_=geometry.build(x,weights,bounds=bounds,face_order=2)
    states=[decode(p,r) for p,r in zip(cells,raw)];flux=np.zeros(len(x))
    factors=np.asarray(factors)
    if factors.shape!=(len(x),3,3):raise ValueError('invalid primary factors')
    for s,L,m in zip(states,factors,mass):
        if not np.allclose(s['R'],L@L.T,rtol=1e-9,atol=2e-12*max(m,1.)):
            raise ValueError('primary factor mismatch')
    jac=np.zeros((len(x),len(x))) if return_jacobian else None
    for f in faces:
        i,j,n=f['i'],f['j'],f['normal']
        # Reflecting wall: identical densities and opposite relative velocities
        # cancel each paired sample pointwise, including moving-wall reflection.
        if j<0:continue
        delta=x[j]-x[i];dist=np.linalg.norm(delta)
        def wn(p):return (delta@gen[j]-(p-x[i])@(gen[j]-gen[i])+.5*(z[i]-z[j]))/dist
        total=0.;capacity=0.
        for k,sign in ((i,1),(j,-1)):
            s=states[k]
            for offset in ALPHA@np.asarray(factors[k]).T/np.sqrt(mass[k]):
                def speed(p):return (s['v']+(p-s['c'])@s['B'].T+offset)@n-wn(p)
                if return_jacobian:
                    value,cap=integrate_sample(f['vertices'],n,speed,s['density'].value,sign,order,True)
                    total+=value/6;capacity+=cap/6
                else:total+=integrate_sample(f['vertices'],n,speed,s['density'].value,sign,order)/6
        flux[i]-=total;flux[j]+=total
        if return_jacobian:
            c=capacity/(2*dist);jac[i,i]+=c;jac[j,j]+=c;jac[i,j]-=c;jac[j,i]-=c
    if return_diagnostics:
        if return_jacobian:raise ValueError('select diagnostics or Jacobian')
        material=np.array([s['v']+s['B']@(p-s['c']) for s,p in zip(states,x)])
        return flux,dict(generator_defect=float(np.max(abs(gen-material))),mass_root_solved=False)
    return (flux,jac) if return_jacobian else flux


def solve(geometry,sites,weights,raw,factors,generator,initial=None,order=4,
          tolerance=2e-13,callback=None,bounds=None):
    """Newton solve of split face mass balance with an explicit gauge."""
    z=np.zeros(len(sites)) if initial is None else np.array(initial,dtype=float,copy=True)
    z-=z[-1];history=[]
    for iteration in range(40):
        flux,J=evaluate(geometry,sites,weights,raw,factors,(generator,z),order,True,bounds)
        norm=float(np.max(abs(flux)));history.append(norm)
        if callback is not None:callback(iteration,z.copy(),flux.copy())
        if norm<tolerance:return z,dict(residual=norm,iterations=iteration,history=history)
        direction=np.r_[np.linalg.lstsq(J[:-1,:-1],-flux[:-1],rcond=1e-13)[0],0.]
        for backtrack in range(25):
            trial=z+direction*2.**(-backtrack)
            candidate=evaluate(geometry,sites,weights,raw,factors,(generator,trial),order,bounds=bounds)
            if np.max(abs(candidate))<norm:z=trial;break
        else:raise RuntimeError('split mass root line search failed')
    raise RuntimeError('split mass root iteration limit')
