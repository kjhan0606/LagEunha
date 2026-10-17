"""Smooth affine-manifold evolution gate, NOT a general moment hydro solver.
Positive exponential-quadratic density fits transported m,q,Q in each polygon.
Affine velocity from C,S; Pi transported independently; R defect is measured,
never clipped or thermalized. Flux closure is intended only near R=0.
Moving slip boundary follows the manufactured material shear box.
"""
import json,os
from pathlib import Path
import numpy as np
from affine_moment_patch import bounded_mesh,integrate,DOMAIN,P0
from periodic2d_hllc import quadrature
from fixedmass_hllc import fan,state,sample,GAMMA
from weighted_mass_patch import power_cells


def domain(t):
    p=DOMAIN.copy();p[:,0]+=t*(.2+1.1*p[:,1]);return p


def boundary_velocity(point):
    return np.array([.2+1.1*point[1],0.])


def geometry(x,w,t):
    builder=power_cells
    if os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')=='1':
        from active_power_cells import power_cells_active
        builder=power_cells_active
    if os.environ.get('LAGFORCE_SEEDED_POWER_CELLS','0')=='1':
        if os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')!='1':raise ValueError('seeded clipping requires all-constraint active clipping')
        cells=builder(x,w,domain(t),seed_count=12)
    else:cells=builder(x,w,domain(t))
    polys=[np.array(p.exterior.coords[:-1]) for p in cells];edges=[]
    for i,p in enumerate(polys):
        if cells[i].area<=0:raise RuntimeError('empty power cell')
        delta=x-x[i];bound=np.sum(delta*delta,axis=1)+w[i]-w
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-14:continue
            normal=np.array([b[1]-a[1],a[0]-b[0]])/length
            residual=abs(2*delta@((a+b)/2-x[i])-bound);residual[i]=np.inf
            j=int(np.argmin(residual))
            if residual[j]>1e-10:j=-1
            if j<0 or i<j:edges.append((i,j,a,b,normal,length))
    return polys,edges


def features(points,origin,scale):
    r=(points-origin)/scale
    return np.column_stack([r[:,0],r[:,1],r[:,0]**2,r[:,0]*r[:,1],r[:,1]**2])


class Density:
    def __init__(self,poly,m,q,Q):
        from density_quadrature import rule,validate_order
        self.quadrature_order=validate_order(os.environ.get('LAGFORCE_DENSITY_QUADRATURE_ORDER','6'))
        pts,w=quadrature(poly) if self.quadrature_order==6 else rule(poly,self.quadrature_order)
        self.area=w.sum();self.origin=w@pts/self.area;self.scale=np.sqrt(self.area)
        self.points=pts;self.weights=w;self.mass=m
        extended=np.asarray(q).dtype==np.dtype(np.longdouble)
        if extended:
            m=float(m);self.mass=m;o=np.asarray(self.origin,dtype=np.longdouble);c=q/m
            second=Q/m-np.outer(o,c)-np.outer(c,o)+np.outer(o,o)
            target=np.asarray(np.r_[(c-o)/self.scale,second[0,0]/self.area,second[0,1]/self.area,second[1,1]/self.area],dtype=float)
        else:
            c=q/m;second=Q/m-np.outer(self.origin,c)-np.outer(c,self.origin)+np.outer(self.origin,self.origin)
            target=np.r_[(c-self.origin)/self.scale,second[0,0]/self.area,second[0,1]/self.area,second[1,1]/self.area]
        F=features(pts,self.origin,self.scale);theta=np.zeros(5)
        def evaluate(theta):
            exponent=F@theta;shift=exponent.max();weighted=w*np.exp(exponent-shift);Z=weighted.sum()
            probability=weighted/Z;mean=probability@F;centered=F-mean
            H=(centered*probability[:,None]).T@centered
            return mean-target,H,shift+np.log(Z),probability
        for iteration in range(30):
            residual,H,logZ,probability=evaluate(theta);norm=np.max(abs(residual))
            if norm<(1e-15 if extended else 3e-13):break
            step=np.linalg.solve(H,residual)
            for bt in range(24):
                trial=theta-step*2.**(-bt)
                if np.max(abs(evaluate(trial)[0]))<norm:theta=trial;break
            else:raise RuntimeError('positive density moment fit stalled')
        else:raise RuntimeError('positive density moments not fitted')
        self.theta=theta;self.logZ=logZ;self.residual=float(norm);self.iterations=iteration+1
        self.weighted_mass=m*probability
    def value(self,points):
        return self.mass*np.exp(features(np.atleast_2d(points),self.origin,self.scale)@self.theta-self.logZ)


def decode(poly,raw,density=None):
    m=raw[0];q=raw[1:3];Q=raw[3:7].reshape(2,2);p=raw[7:9];D=raw[9:13].reshape(2,2);Pi=raw[13:17].reshape(2,2)
    extended=np.asarray(raw).dtype==np.dtype(np.longdouble)
    if extended:m=float(m)
    if density is None:density=Density(poly,m,q,Q)
    elif density.mass!=float(m):raise ValueError('provided density has different cell mass')
    c=q/m;v=p/m;S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
    S=(S+S.T)/2
    if extended:
        c=np.asarray(c,dtype=float);v=np.asarray(v,dtype=float);S=np.asarray(S,dtype=float);C=np.asarray(C,dtype=float)
    if np.linalg.eigvalsh(S).min()<=0:raise RuntimeError('nonpositive position covariance')
    B=np.linalg.solve(S,C).T;T=Pi-np.outer(p,p)/m;R=T-B@S@B.T;R=(R+R.T)/2
    pressure=(GAMMA-1)*(raw[17]-.5*np.trace(Pi))/density.area
    if pressure<=0 or not np.isfinite(pressure):raise RuntimeError('nonpositive thermal pressure')
    if extended:R=np.asarray(R,dtype=float);pressure=float(pressure)
    return dict(density=density,c=c,v=v,B=B,R=R,P=pressure,raw=raw)


def affine_velocity(s,point):return s['v']+s['B']@(point-s['c'])


def velocity(s,point):
    value=affine_velocity(s,point)
    if 'velocity_quadratic' in s:
        from quadratic_velocity_reconstruction import correction
        value=value+correction(s,point)[0]
    return value


def solve(records,count):
    def evaluate(z):
        residual=np.zeros(count);lap=np.zeros((count,count));values=[]
        for i,j,length,base,distance,f in records:
            wn=base+(z[i]-z[j])/(2*distance) if j>=0 else base
            flux,rho,branch=sample(f,wn);values.append((flux,branch,wn))
            residual[i]+=length*flux[0]
            if j>=0:
                residual[j]-=length*flux[0];c=length*rho/(2*distance)
                lap[i,i]+=c;lap[j,j]+=c;lap[i,j]-=c;lap[j,i]-=c
        return residual,lap,values
    z=np.zeros(count)
    for iteration in range(30):
        r,L,values=evaluate(z);norm=np.max(abs(r))
        if norm<2e-13:return z,values,norm
        dz=np.r_[np.linalg.solve(L[:-1,:-1],r[:-1]),0.]
        for bt in range(24):
            trial=z+dz*2.**(-bt)
            if np.max(abs(evaluate(trial)[0]))<norm:z=trial;break
        else:raise RuntimeError('weight constraint incompatible/stalled')
    raise RuntimeError('weight constraint failed')


def initial(n):
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([xx.ravel(),yy.ravel()])+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    polys,_=geometry(x,np.zeros(n*n),0.)
    exact=lambda p:np.column_stack([.2+1.1*p[:,1],np.zeros(len(p))])
    raw=np.array([integrate(poly,exact)['raw'] for poly in polys]);mass=raw[:,0].copy()
    y=np.column_stack([x,np.zeros(n*n),raw[:,1:]])
    return mass,np.r_[y.ravel(),np.zeros(4)]


def rhs(t,flat,mass,face_density_provider=None,face_segments=1,face_order=4):
    count=len(mass);y=flat[:-4].reshape(count,20);x=y[:,:2];w=y[:,2]
    polys,edges=geometry(x,w,t);raw=np.column_stack([mass,y[:,3:]])
    states=[decode(poly,r) for poly,r in zip(polys,raw)]
    trace=None;density_diag={}
    if face_density_provider is not None:trace,density_diag=face_density_provider(polys,edges,raw)
    gen=np.array([velocity(s,point) for s,point in zip(states,x)])
    if face_segments<1 or face_order<1:raise ValueError('invalid face quadrature')
    records=[];data=[];gx,gw=np.polynomial.legendre.leggauss(face_order);gx=(gx+1)/2;gw=gw/2
    nodes=[((j+u)/face_segments,weight/face_segments) for j in range(face_segments) for u,weight in zip(gx,gw)]
    for edge_id,(i,j,a,b,normal,length) in enumerate(edges):
        distance=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        for u,weight in nodes:
            point=a+u*(b-a);vi=velocity(states[i],point);rhoi=float(states[i]['density'].value(point)[0]) if trace is None else trace(i,edge_id,u);pi=states[i]['P']
            if j>=0:
                vj=velocity(states[j],point);rhoj=float(states[j]['density'].value(point)[0]) if trace is None else trace(j,edge_id,u);pj=states[j]['P']
                base=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/distance
            else:
                base=np.dot(boundary_velocity(point),normal)
                vj=vi-2*(vi@normal-base)*normal;rhoj=rhoi;pj=pi
            f=fan(state(rhoi,vi,pi),state(rhoj,vj,pj),normal)
            records.append((i,j,length*weight,base,distance,f));data.append((point,normal))
    z,values,mass_residual=solve(records,count);derivative=np.zeros((count,18));ledger=np.zeros(4);area_rates=np.zeros(count)
    for rec,value,info in zip(records,values,data):
        i,j,length,base,distance,f=rec;flux,branch,wn=value;point,normal=info
        area_rates[i]+=length*wn
        if j>=0:area_rates[j]-=length*wn
        uf=f[1][branch];vf=uf[1:3]/uf[0];fm=flux[0];fp=flux[1:3]
        pressure_face=np.dot(fp-vf*fm,normal)
        fpi=np.outer(vf,vf)*fm+pressure_face*(np.outer(vf,normal)+np.outer(normal,vf))
        ff=np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,np.outer(point,fp).ravel(),fpi.ravel(),flux[3]]*length
        derivative[i]-=ff
        if j>=0:derivative[j]+=ff
        else:ledger-=np.r_[fp,flux[3],point[0]*fp[1]-point[1]*fp[0]]*length
    for i,s in enumerate(states):
        r=raw[i];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2);A=s['density'].area
        derivative[i,1:3]+=r[7:9];derivative[i,3:7]+=(D+D.T).ravel()
        derivative[i,9:13]+=(Pi+s['P']*A*np.eye(2)).ravel()
        derivative[i,13:17]+=(s['P']*A*(s['B']+s['B'].T)).ravel()
    # Pi closure has no modeled finite R stress. Record its departure from the
    # affine manifold; do not clip, reset, or disguise the defect as heat.
    diag=dict(mass_flux=mass_residual,min_pressure=min(s['P'] for s in states),
        density_fit_residual=max(s['density'].residual for s in states),
        density_fit_iterations=max(s['density'].iterations for s in states),
        residual_covariance_min=min(float(np.linalg.eigvalsh(s['R']).min()/mass[i]) for i,s in enumerate(states)),
        residual_covariance_max=max(float(np.max(abs(s['R']))/mass[i]) for i,s in enumerate(states)))
    area=np.array([s['density'].area for s in states]);P=np.array([s['P'] for s in states])
    diag.update(density_diag,mass_rates=derivative[:,0].copy(),area_rates=area_rates,
        generator_velocity_defect=0.,
        pressure_rates=(GAMMA-1)*(derivative[:,17]-.5*(derivative[:,13]+derivative[:,16]))/area-P*area_rates/area)
    out=np.column_stack([gen,z,derivative[:,1:]])
    return np.r_[out.ravel(),ledger],diag


def diagnostics(t,flat,mass,initial_flat):
    y=flat[:-4].reshape(len(mass),20);y0=initial_flat[:-4].reshape(len(mass),20)
    raw=np.column_stack([mass,y[:,3:]]);r0=np.column_stack([mass,y0[:,3:]])
    polys,edges=geometry(y[:,:2],y[:,2],t);states=[decode(p,r) for p,r in zip(polys,raw)]
    true=lambda p:np.column_stack([.2+1.1*p[:,1],np.zeros(len(p))])
    reference=np.array([integrate(p,true)['raw'] for p in polys])
    expected_x=y0[:,:2].copy();expected_x[:,0]+=t*(.2+1.1*expected_x[:,1])
    ptotal=raw[:,7:9].sum(0);pstart=r0[:,7:9].sum(0)
    spin=np.sum(raw[:,10]-raw[:,11]);spin0=np.sum(r0[:,10]-r0[:,11])
    return dict(time=t,topology=sorted((e[0],e[1]) for e in edges if e[1]>=0),
        pressure_error=max(abs(s['P']-P0) for s in states),
        density_error=max(float(np.max(abs(s['density'].value(s['density'].points)-1))) for s in states),
        generator_error=float(np.max(abs(y[:,:2]-expected_x))),raw_moment_error=float(np.max(abs(raw-reference))),
        momentum_ledger_error=float(np.max(abs(ptotal-pstart-flat[-4:-2]))),
        energy_ledger_error=float(abs(raw[:,17].sum()-r0[:,17].sum()-flat[-2])),
        angular_ledger_error=float(abs(spin-spin0-flat[-1])),
        min_covariance=min(float(np.linalg.eigvalsh(s['R']).min()/mass[i]) for i,s in enumerate(states)),
        max_covariance_defect=max(float(np.max(abs(s['R']))/mass[i]) for i,s in enumerate(states)))


def run(n,dt,end=.02):
    mass,y=initial(n);y0=y.copy();steps=round(end/dt)
    stats=dict(mass_flux=0.,min_pressure=float('inf'),density_fit_residual=0.,density_fit_iterations=0,
               residual_covariance_min=0.,residual_covariance_max=0.)
    def observed(t,y):
        dy,d=rhs(t,y,mass)
        for k in stats:stats[k]=(min if k in ('min_pressure','residual_covariance_min') else max)(stats[k],d[k])
        return dy
    name=f'n{n}_dt{dt}_end{end}';history=[]
    previous=diagnostics(0.,y,mass,y0)['topology']
    for step in range(steps):
        t=step*dt;k1=observed(t,y);k2=observed(t+dt/2,y+dt*k1/2)
        k3=observed(t+dt/2,y+dt*k2/2);k4=observed(t+dt,y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        d=diagnostics((step+1)*dt,y,mass,y0)
        d['topology_changed']=d['topology']!=previous;previous=d['topology'];history.append(d)
        Path(name+'.checkpoint.json').write_text(json.dumps(dict(time=(step+1)*dt,mass=mass.tolist(),state=y.tolist())))
        print(name,'step',step+1,'/',steps,'Perr',d['pressure_error'],'Rmin',d['min_covariance'],flush=True)
    observed(end,y)
    Path(name+'.history.json').write_text(json.dumps(history,indent=2))
    return dict(n=n,dt=dt,end=end,steps=steps,**stats,final=history[-1],
        covariance_monitor_tolerance=2e-11,
        final_covariance_above_monitor_threshold=history[-1]['min_covariance']>=-2e-11,
        max_completed_covariance_defect=max(d['max_covariance_defect'] for d in history))


if __name__=='__main__':
    failures=[]
    with open('results.jsonl','x') as out:
        for n in (3,4):
            for dt in (.005,.0025,.00125):
                try:r=run(n,dt);r['status']='complete'
                except Exception as error:r=dict(n=n,dt=dt,status='failed',error=repr(error));failures.append(r)
                out.write(json.dumps(r)+'\n');out.flush()
    if failures:raise RuntimeError(f'{len(failures)} cases failed')
