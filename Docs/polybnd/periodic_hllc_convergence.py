"""1D periodic power-cell reduction of fixed-mass ALE; not a 2D mesh test.

Unwrapped ordered generators, ring neighbors, circumference one. No steering.
Analytic pre-shock isentropic simple wave via characteristic inversion.
"""
import json
import numpy as np
from fixedmass_hllc import fan,solve_weights,GAMMA

GAUSS_X,GAUSS_W=np.polynomial.legendre.leggauss(12)


def geometry(x,w):
    distance=np.roll(x,-1)-x; distance[-1]+=1.
    if np.any(distance<=0): raise RuntimeError('generators crossed')
    right=x+distance/2+(w-np.roll(w,-1))/(2*distance)
    left=np.roll(right,1); left[0]-=1.
    width=right-left
    if np.any(width<=0): raise RuntimeError('power cell vanished')
    return left,right,width,distance


def exact(points,t,kind,boost):
    if kind=='wave':
        xi=points-(boost+1)*t
        for _ in range(12):
            r=1+.01*np.sin(2*np.pi*xi); c=r**((GAMMA-1)/2)
            dc=(GAMMA-1)/2*r**((GAMMA-3)/2)*.01*2*np.pi*np.cos(2*np.pi*xi)
            speed=boost+1+(GAMMA+1)/(GAMMA-1)*(c-1)
            derivative=1+t*(GAMMA+1)/(GAMMA-1)*dc
            if np.any(derivative<=0): raise RuntimeError('characteristics crossed')
            step=(xi+t*speed-points)/derivative; xi-=step
            if np.max(abs(step))<2e-15: break
        r=1+.01*np.sin(2*np.pi*xi); c=r**((GAMMA-1)/2)
        v=boost+2*(c-1)/(GAMMA-1); p=r**GAMMA/GAMMA
        assert np.max(abs(xi+t*(v+c)-points))<2e-13
    elif kind=='smooth_contact':
        r=1+.2*np.sin(2*np.pi*(points-boost*t)); v=np.full_like(r,boost); p=np.full_like(r,1/GAMMA)
    elif kind=='contact':
        phase=(points-boost*t)%1
        r=np.where((phase>=.25)&(phase<.75),2.,1.); v=np.full_like(r,boost); p=np.full_like(r,1/GAMMA)
    else: raise ValueError(kind)
    return np.stack([r,r*v,p/(GAMMA-1)+.5*r*v*v],axis=-1)


def average(left,right,t,kind,boost):
    points=(left[:,None]+right[:,None])/2+(right-left)[:,None]*GAUSS_X/2
    return np.einsum('ijk,j->ik',exact(points,t,kind,boost),GAUSS_W)/2


def initial(n,kind,boost):
    x=(np.arange(n)+.5)/n; w=np.zeros(n); left,right,a,_=geometry(x,w)
    u=average(left,right,0.,kind,boost)
    mass=u[:,0]*a
    return mass,np.column_stack([x,w,u[:,1]*a,u[:,2]*a])


def rhs(y,mass):
    x,w,p,e=y.T; left,right,a,d=geometry(x,w)
    v=p/mass; rho=mass/a; pressure=(GAMMA-1)*(e-.5*p*p/mass)/a
    if np.any(pressure<=0): raise RuntimeError('nonpositive pressure')
    u=np.column_stack([rho,p/a,np.zeros(len(x)),e/a]); records=[]
    base=(v+np.roll(v,-1))/2-(w-np.roll(w,-1))*(np.roll(v,-1)-v)/(2*d*d)
    normal=np.array([1.,0.])
    for i in range(len(x)):
        j=(i+1)%len(x)
        records.append((i,j,1.,base[i],d[i],fan(u[i],u[j],normal)))
    zdot,samples,iterations,residual=solve_weights(records)
    flux=np.array([s[0] for s in samples]); change=np.roll(flux,1,axis=0)-flux
    return np.column_stack([v,zdot,change[:,1],change[:,3]]),dict(
        mass_flux=float(np.max(abs(change[:,0]))),pressure_min=float(pressure.min()),
        width_min=float(a.min()),iterations=iterations,edge_flux=float(np.max(abs(flux[:,0]))))


def run(n,kind,boost,cfl=.15,end=.2):
    mass,y=initial(n,kind,boost); y0=y.copy(); momentum0=y[:,2].sum(); energy0=y[:,3].sum()
    steps=int(np.ceil(end*n/cfl)); dt=end/steps
    diagnostics=dict(mass_flux=0.,pressure_min=float('inf'),width_min=float('inf'),iterations=0,edge_flux=0.)
    for _ in range(steps):
        k1,d1=rhs(y,mass); k2,d2=rhs(y+dt*k1/2,mass)
        k3,d3=rhs(y+dt*k2/2,mass); k4,d4=rhs(y+dt*k3,mass)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        for d in (d1,d2,d3,d4):
            for key in diagnostics:
                diagnostics[key]=(min if key in ('pressure_min','width_min') else max)(diagnostics[key],d[key])
    _,df=rhs(y,mass)
    for key in diagnostics:
        diagnostics[key]=(min if key in ('pressure_min','width_min') else max)(diagnostics[key],df[key])
    left,right,a,_=geometry(y[:,0],y[:,1]); u=np.column_stack([mass,y[:,2],y[:,3]])/a[:,None]
    reference=average(left,right,end,kind,boost)
    error=np.sum(a[:,None]*abs(u-reference),axis=0)
    points=(left[:,None]+right[:,None])/2+a[:,None]*GAUSS_X/2
    real=exact(points,end,kind,boost)
    reconstructed_l1=float(np.sum(a*np.einsum('ij,j->i',abs(u[:,None,0]-real[:,:,0]),GAUSS_W)/2))
    # Exact integral of the Fourier mode of the piecewise-constant density.
    integral=(np.exp(-2j*np.pi*right)-np.exp(-2j*np.pi*left))/(-2j*np.pi)
    mode=np.sum(u[:,0]*integral)
    refmode=np.sum(a*np.einsum('ij,j->i',real[:,:,0]*np.exp(-2j*np.pi*points),GAUSS_W)/2)
    phase=float(np.angle(mode/refmode)) if abs(refmode)>1e-12 else None
    result=dict(n=n,kind=kind,boost=boost,cfl=cfl,steps=steps,dt=dt,end=end,**diagnostics,
        cell_average_l1=error.tolist(),reconstructed_density_l1=reconstructed_l1,
        mode_amplitude_ratio=float(abs(mode/refmode)) if abs(refmode)>1e-12 else None,
        mode_phase_error=phase,mode_speed_error=-phase/(2*np.pi*end) if phase is not None else None,
        energy_relative=float(y[:,3].sum()/energy0-1),momentum_error=float(y[:,2].sum()-momentum0),
        translated_generator_error=float(np.max(abs(y[:,0]-y0[:,0]-boost*end))),
        weight_change=float(np.max(abs(y[:,1]-y0[:,1]))),final_state=y.tolist())
    assert diagnostics['mass_flux']<3e-13
    assert abs(result['energy_relative'])<1e-11 and abs(result['momentum_error'])<1e-11
    if kind in ('contact','smooth_contact'):
        assert error[0]<1e-10 and result['translated_generator_error']<1e-10 and result['weight_change']<1e-10
    return result


if __name__=='__main__':
    for kind,boost in [('contact',.7),('smooth_contact',.7),('wave',0.),('wave',.7)]:
        for n in (16,32,64,128): print(json.dumps(run(n,kind,boost)),flush=True)
    for boost in (0.,.7): print(json.dumps(run(128,'wave',boost,cfl=.075)),flush=True)
