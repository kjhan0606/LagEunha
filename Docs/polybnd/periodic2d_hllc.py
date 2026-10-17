"""Small 2D periodic power-cell ALE reference. Exact piecewise-linear face flux.
No regularization, no mass projection, piecewise constant states, fixed masses.
"""
import json
import numpy as np
from fixedmass_hllc import fan,sample,GAMMA
from weighted_mass_patch import clip_halfplane

SHIFTS=np.array([(a,b) for a in (-1,0,1) for b in (-1,0,1)])
K=np.array([1.,1.]); KN=np.linalg.norm(K); NORMAL=K/KN
GX,GW=np.polynomial.legendre.leggauss(6); GX=(GX+1)/2; GW=GW/2


def mesh(x,w,audit=False):
    x=x%1; count=len(x); polygons=[]; areas=[]; edges={}
    for i in range(count):
        delta=(x[None,:,:]+SHIFTS[:,None,:]-x[i]).reshape(-1,2)
        ids=np.tile(np.arange(count),9); shifts=np.repeat(SHIFTS,count,axis=0)
        norm=np.linalg.norm(delta,axis=1); bound=np.sum(delta*delta,axis=1)+w[i]-w[ids]
        order=np.argsort(norm); order=order[norm[order]>1e-14]
        p=np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]])
        # Self images bound the cell in this square. Each subsequent exclusion
        # uses a circle containing the CURRENT polygon, so it cannot drop a cut.
        pending=order
        while len(pending):
            radius=np.max(np.linalg.norm(p,axis=1))
            pending=pending[bound[pending]<2*norm[pending]*radius+1e-14]
            if not len(pending): break
            j=pending[0]; pending=pending[1:]
            p=clip_halfplane(p,2*delta[j],bound[j])
            if len(p)<3: raise RuntimeError('empty power cell')
        cross=p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)
        area=.5*cross.sum()
        if area<=0: raise RuntimeError('nonpositive signed area')
        polygons.append(p+x[i]); areas.append(area)
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-14: continue
            mid=(a+b)/2
            residual=abs(2*delta@mid-bound); residual[norm<1e-14]=np.inf
            j=int(np.argmin(residual))
            if residual[j]>2e-11: raise RuntimeError('unmatched periodic face')
            other=int(ids[j]); shift=tuple(int(s) for s in shifts[j]); key=(i,other,*shift)
            if other==i: raise RuntimeError('self face: grid too small')
            edges[key]=(a,b,delta[j],length)
    if abs(sum(areas)-1)>2e-10: raise RuntimeError('periodic area closure')
    output=[]
    for (i,j,sx,sy),(a,b,d,length) in edges.items():
        if i>j: continue
        reverse=edges.get((j,i,-sx,-sy))
        if reverse is None: raise RuntimeError('missing reverse face')
        aa,bb,dd,ll=reverse
        if abs(ll-length)>2e-10 or max(np.linalg.norm(aa+d-b),np.linalg.norm(bb+d-a))>2e-10:
            raise RuntimeError('periodic faces disagree')
        output.append((i,j,a,b,d,length))
    if audit:
        closure=np.zeros((count,2))
        for i,j,a,b,d,length in output:
            n=d/np.linalg.norm(d); closure[i]+=length*n; closure[j]-=length*n
        if np.max(abs(closure))>2e-10: raise RuntimeError('cell normal closure')
    return polygons,np.array(areas),output


def integrate_face(f,b0,b1,offset,length):
    s0=b0+offset; s1=b1+offset; points=[0.,1.]
    if s1!=s0:
        points.extend((speed-s0)/(s1-s0) for speed in f[0] if 0<(speed-s0)/(s1-s0)<1)
    points=sorted(points); flux=np.zeros(4); density=0.
    for lo,hi in zip(points,points[1:]):
        value,rho,_=sample(f,s0+(s1-s0)*(lo+hi)/2)
        flux+=(hi-lo)*value; density+=(hi-lo)*rho
    return length*flux,length*density,len(points)-1


def solve(records,count):
    def evaluate(z):
        residual=np.zeros(count); lap=np.zeros((count,count)); values=[]; cuts=0
        for i,j,distance,length,b0,b1,f in records:
            flux,rho,pieces=integrate_face(f,b0,b1,(z[i]-z[j])/(2*distance),length)
            residual[i]+=flux[0]; residual[j]-=flux[0]; c=rho/(2*distance)
            lap[i,i]+=c; lap[j,j]+=c; lap[i,j]-=c; lap[j,i]-=c
            values.append(flux); cuts+=pieces>1
        return residual,lap,values,cuts
    z=np.zeros(count)
    for iteration in range(40):
        r,lap,values,cuts=evaluate(z); norm=np.max(abs(r))
        if norm<2e-13: return z,values,iteration+1,cuts
        dz=np.r_[np.linalg.solve(lap[:-1,:-1],r[:-1]),0.]
        for backtrack in range(25):
            trial=z+dz*2.**(-backtrack)
            if np.max(abs(evaluate(trial)[0]))<norm:
                z=trial; break
        else: raise RuntimeError('line search failed')
    raise RuntimeError('weight constraint failed')


def exact(points,t,kind,boost):
    shape=points.shape[:-1]; q=(points-t*boost)@K
    if kind=='wave':
        q0=q-KN*t
        for _ in range(15):
            rho=1+.01*np.sin(2*np.pi*q0); c=rho**((GAMMA-1)/2)
            dc=(GAMMA-1)/2*rho**((GAMMA-3)/2)*.01*2*np.pi*np.cos(2*np.pi*q0)
            derivative=1+t*KN*(GAMMA+1)/(GAMMA-1)*dc
            if np.any(derivative<=0): raise RuntimeError('characteristics crossed')
            step=(q0+t*KN*(1+(GAMMA+1)/(GAMMA-1)*(c-1))-q)/derivative
            q0-=step
            if np.max(abs(step))<2e-15: break
        rho=1+.01*np.sin(2*np.pi*q0); c=rho**((GAMMA-1)/2)
        v=boost+2*(c[...,None]-1)/(GAMMA-1)*NORMAL; p=rho**GAMMA/GAMMA
    else:
        rho=1+.2*np.sin(2*np.pi*q) if kind=='contact' else np.ones(shape)
        v=np.broadcast_to(boost,shape+(2,)); p=np.full(shape,1/GAMMA)
    return np.concatenate([rho[...,None],rho[...,None]*v,(p/(GAMMA-1)+.5*rho*np.sum(v*v,axis=-1))[...,None]],axis=-1)


def quadrature(p):
    center=p.mean(axis=0); points=[]; weights=[]
    for a,b in zip(p,np.roll(p,-1,axis=0)):
        u=a-center; v=b-center; determinant=abs(u[0]*v[1]-u[1]*v[0])
        for r,wr in zip(GX,GW):
            points.extend(center+r*u+((1-r)*GX)[:,None]*v)
            weights.extend(determinant*(1-r)*wr*GW)
    return np.array(points),np.array(weights)


def averages(polygons,t,kind,boost):
    result=[]
    for p in polygons:
        q,w=quadrature(p); result.append(w@exact(q,t,kind,boost)/w.sum())
    return np.array(result)


def initial(n,kind,boost):
    a,b=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([a.ravel(),b.ravel()]); x+=np.random.default_rng(613).uniform(-.18,.18,x.shape)/n
    weights=np.zeros(n*n); polygons,area,_=mesh(x,weights,audit=True)
    u=averages(polygons,0.,kind,boost); mass=area*u[:,0]
    return mass,np.column_stack([x,weights,area[:,None]*u[:,1:3],area*u[:,3]])


def rhs(y,mass):
    x=y[:,:2]; w=y[:,2]; momentum=y[:,3:5]; energy=y[:,5]
    polygons,area,edges=mesh(x,w); velocity=momentum/mass[:,None]
    pressure=(GAMMA-1)*(energy-.5*np.sum(momentum*momentum,axis=1)/mass)/area
    if not np.isfinite(pressure).all() or np.any(pressure<=0): raise RuntimeError('nonpositive pressure')
    u=np.column_stack([mass/area,momentum/area[:,None],energy/area]); records=[]
    for i,j,a,b,delta,length in edges:
        distance=np.linalg.norm(delta); normal=delta/distance; dv=velocity[j]-velocity[i]
        b0=(delta@velocity[j]-dv@a)/distance; b1=(delta@velocity[j]-dv@b)/distance
        records.append((i,j,distance,length,b0,b1,fan(u[i],u[j],normal)))
    z,values,iterations,cuts=solve(records,len(x)); change=np.zeros((len(x),4))
    for record,value in zip(records,values):
        i,j=record[:2]; change[i]-=value; change[j]+=value
    return np.column_stack([velocity,z,change[:,1:3],change[:,3]]),dict(
        mass_flux=float(np.max(abs(change[:,0]))),pressure_min=float(pressure.min()),area_min=float(area.min()),
        iterations=iterations,split_faces=cuts,edge_mass_flux=float(max(abs(v[0]) for v in values)),
        topology=tuple(sorted((e[0],e[1]) for e in edges)))


def run(n,kind,boost,refine=1):
    boost=np.array(boost); mass,y=initial(n,kind,boost); y0=y.copy(); end=.1
    steps=int(np.ceil(end*n/.12))*refine; dt=end/steps
    stats=dict(mass_flux=0.,pressure_min=float('inf'),area_min=float('inf'),iterations=0,split_faces=0,edge_mass_flux=0.)
    previous=None; topology_changes=0
    def observed(y):
        nonlocal previous,topology_changes
        dy,d=rhs(y,mass)
        if previous is not None and previous!=d['topology']: topology_changes+=1
        previous=d['topology']
        for k in stats: stats[k]=(min if k in ('pressure_min','area_min') else max)(stats[k],d[k])
        return dy
    for _ in range(steps):
        k1=observed(y); k2=observed(y+dt*k1/2); k3=observed(y+dt*k2/2); k4=observed(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
    observed(y)
    polygons,area,_=mesh(y[:,:2],y[:,2],audit=True); u=np.column_stack([mass,y[:,3:5],y[:,5]])/area[:,None]
    ref=averages(polygons,end,kind,boost); errors=np.sum(area[:,None]*abs(u-ref),axis=0)
    mode=0j; refmode=0j
    for i,p in enumerate(polygons):
        q,w=quadrature(p); factor=np.exp(-2j*np.pi*(q@K)); mode+=u[i,0]*(w@factor)
        refmode+=w@(factor*exact(q,end,kind,boost)[:,0])
    phase=float(np.angle(mode/refmode)) if abs(refmode)>1e-10 else None
    result=dict(n=n,kind=kind,boost=boost.tolist(),refine=refine,end=end,steps=steps,dt=dt,**stats,
        topology_changes=topology_changes,cell_average_l1=errors.tolist(),
        energy_relative=float(y[:,5].sum()/y0[:,5].sum()-1),momentum_error=float(np.max(abs(y[:,3:5].sum(0)-y0[:,3:5].sum(0)))),
        translated_generator_error=float(np.max(abs(y[:,:2]-y0[:,:2]-end*boost))),weight_change=float(np.max(abs(y[:,2]))),
        amplitude_ratio=float(abs(mode/refmode)) if phase is not None else None,phase_error=phase,
        final_density=u[:,0].tolist(),final_state=y.tolist())
    assert stats['mass_flux']<3e-13 and abs(result['energy_relative'])<1e-10 and result['momentum_error']<1e-10
    if kind!='wave': assert errors[0]<1e-9 and result['translated_generator_error']<1e-9 and result['weight_change']<1e-9
    return result


if __name__=='__main__':
    cases=[(8,'uniform',[.37,-.23],1)]
    cases += [(n,kind,[.37,-.23] if kind=='contact' else [0.,0.],1) for kind in ('contact','wave') for n in (4,8,16)]
    cases += [(8,'wave',[.37,-.23],1),(16,'wave',[0.,0.],2)]
    for args in cases: print(json.dumps(run(*args)),flush=True)
