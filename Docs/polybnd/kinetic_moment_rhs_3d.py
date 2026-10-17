"""3D affine kinetic moment RHS on actual weighted C cells.

First coupled spatial operator: cell-constant pressure, affine mean velocity,
positive covariance samples and instantaneous fixed-mass weight rates. It does
not yet provide the production quadratic reconstruction or a time driver.
"""
import numpy as np
from moment_state_3d import decode,unpack,face_flux,ALPHA
from euler_star_pressure import star


def rhs(geometry,sites,weights,raw,factors,bounds=None,boundary_velocity=None,face_order=4,pressure_reconstruction='constant',audit_geometry_rates=None,face_integration='point'):
    x=np.asarray(sites);raw=np.asarray(raw);N=len(x);mass=raw[:,0]
    polys,faces,gdiag=geometry.build(x,weights,bounds,face_order)
    states=[decode(poly,r) for poly,r in zip(polys,raw)]
    if pressure_reconstruction not in ('constant','linear','quadratic'):raise ValueError('unsupported 3D pressure reconstruction')
    pressure_polynomials=None;pressure_diag={}
    if pressure_reconstruction=='quadratic':
        from pressure_reconstruction_3d import prepare,value
        pressure_polynomials,pressure_diag=prepare(states,polys)
    centers=np.array([s['density'].origin for s in states]);gradients=np.zeros((N,3))
    if pressure_reconstruction=='linear':
        neighbors=[set() for _ in range(N)]
        for f in faces:
            i,j=f['i'],f['j']
            if j>=0:neighbors[i].add(j);neighbors[j].add(i)
        for i,ids in enumerate(neighbors):
            ids=list(ids);dr=centers[ids]-centers[i];length=np.linalg.norm(dr,axis=1)
            dp=np.array([states[j]['P']-states[i]['P'] for j in ids])
            gradient,_,rank,_=np.linalg.lstsq(dr/length[:,None],dp/length,rcond=None)
            if rank<3:raise RuntimeError('3D linear pressure stencil deficient')
            gradients[i]=gradient
    def pressure(i,point):
        if pressure_polynomials is not None:
            from pressure_reconstruction_3d import value
            result=float(value(pressure_polynomials[i],point)[0])
        else:result=states[i]['P']+gradients[i]@(point-centers[i])
        if result<=0:raise RuntimeError('nonpositive 3D reconstructed pressure')
        return result
    factors=np.asarray(factors)
    if factors.shape!=(N,3,3):raise ValueError('primary 3D covariance factors required')
    for s,L,m in zip(states,factors,mass):
        if not np.allclose(s['R'],L@L.T,rtol=1e-9,atol=2e-12*max(m,1.)):
            raise ValueError('3D raw covariance and primary factors disagree')
    def velocity(i,p):return states[i]['v']+states[i]['B']@(p-states[i]['c'])
    material_gen=np.array([velocity(i,p) for i,p in enumerate(x)]);gen=material_gen
    if audit_geometry_rates is not None:
        gen=np.asarray(audit_geometry_rates[0]);prescribed_z=np.asarray(audit_geometry_rates[1])
        if gen.shape!=(N,3) or prescribed_z.shape!=(N,) or not np.isfinite(gen).all() or not np.isfinite(prescribed_z).all():raise ValueError('invalid 3D diagnostic geometry rates')
    boundary_velocity=(lambda point:np.zeros(3)) if boundary_velocity is None else boundary_velocity
    if face_integration not in ('point','split'):raise ValueError('unknown face integration')
    if face_integration=='split':
        from split_face_mass_3d import solve as solve_split
        from split_face_rules_3d import rule
        if audit_geometry_rates is None:
            split_z,split_diag=solve_split(geometry,x,weights,raw,factors,gen,order=face_order,bounds=bounds)
        else:split_z=prescribed_z;split_diag=dict(iterations=0)
        split_faces=[]
        for face in faces:
            i,j,n=face['i'],face['j'],face['normal']
            if j<0:
                if any(float(boundary_velocity(p)@n)!=0. for p in np.vstack([face['vertices'],face['points']])):
                    raise ValueError('split face mode currently requires fixed reflecting walls')
                def wall_speeds(points):
                    mean=states[i]['v']+(points-states[i]['c'])@states[i]['B'].T
                    samples=mean[:,None,:]+(ALPHA@factors[i].T/np.sqrt(mass[i]))[None,:,:]
                    return samples@n
                points,fw=rule(face['vertices'],n,wall_speeds,face_order)
                split_faces.append(dict(face,points=points,weights=fw));continue
            delta=x[j]-x[i];distance=np.linalg.norm(delta)
            def speeds(points):
                wn=(delta@gen[j]-(points-x[i])@(gen[j]-gen[i])+.5*(split_z[i]-split_z[j]))/distance
                values=[]
                for k in (i,j):
                    mean=states[k]['v']+(points-states[k]['c'])@states[k]['B'].T
                    samples=mean[:,None,:]+(ALPHA@factors[k].T/np.sqrt(mass[k]))[None,:,:]
                    values.append(samples@n-wn[:,None])
                return np.concatenate(values,axis=1)
            points,fw=rule(face['vertices'],n,speeds,face_order)
            split_faces.append(dict(face,points=points,weights=fw))
        faces=split_faces
    records=[]
    for f in faces:
        i,j,n=f['i'],f['j'],f['normal'];d=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        for point,scale in zip(f['points'],f['weights']):
            vl=velocity(i,point);rl=float(states[i]['density'].value(point)[0]);pl=pressure(i,point);vls=vl+ALPHA@factors[i].T/np.sqrt(mass[i])
            if j>=0:
                vr=velocity(j,point);rr=float(states[j]['density'].value(point)[0]);pr=pressure(j,point);vrs=vr+ALPHA@factors[j].T/np.sqrt(mass[j])
                wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
            else:
                wn=float(boundary_velocity(point)@n);vr=vl-2*(vl@n-wn)*n;rr=rl;pr=pl;vrs=vls-2*(vls@n-wn)[:,None]*n
            ps,sm=star(rl,vl,pl,rr,vr,pr,n)
            records.append((i,j,scale,d,point,n,wn,rl,rr,vls,vrs,pl,pr,ps,sm))
    def evaluate(z):
        residual=np.zeros(N);jac=np.zeros((N,N));values=[]
        for i,j,scale,d,point,n,wn,rl,rr,vl,vr,pl,pr,ps,sm in records:
            wn=wn+(z[i]-z[j])/(2*d) if j>=0 else wn
            sl=vl@n-wn;sr=vr@n-wn;a=rl/6*np.maximum(sl,0.);b=rr/6*np.minimum(sr,0.)
            values.append((a,b,wn));residual[i]+=scale*(a.sum()+b.sum())
            if j>=0:
                residual[j]-=scale*(a.sum()+b.sum())
                deriv=rl/6*np.where(sl>0,1.,np.where(sl<0,0.,.5)).sum()+rr/6*np.where(sr<0,1.,np.where(sr>0,0.,.5)).sum()
                c=scale*deriv/(2*d);jac[i,i]+=c;jac[j,j]+=c;jac[i,j]-=c;jac[j,i]-=c
        return residual,jac,values
    if face_integration=='split':
        z=split_z;residual,_,values=evaluate(z);norm=max(abs(residual));iteration=split_diag['iterations']
        if audit_geometry_rates is None and norm>=2e-13:
            raise RuntimeError('split shared-face mass quadrature disagrees with root')
    elif audit_geometry_rates is not None:
        z=prescribed_z;residual,_,values=evaluate(z);norm=max(abs(residual));iteration=0
    else:
        z=np.zeros(N)
        for iteration in range(40):
            residual,J,values=evaluate(z);norm=max(abs(residual))
            if norm<2e-13:break
            step=np.r_[np.linalg.lstsq(J[:-1,:-1],residual[:-1],rcond=1e-13)[0],0.]
            for bt in range(30):
                trial=z+step*2.**(-bt)
                if max(abs(evaluate(trial)[0]))<norm:z=trial;break
            else:raise RuntimeError('3D mass weight solve stalled')
        else:raise RuntimeError('3D mass weight solve exhausted')
    dr=np.zeros_like(raw);ledger=np.zeros(7);F=np.zeros((N,3));H=np.zeros((N,3,3));K=np.zeros((N,3,3));incoming=[[] for _ in range(N)]
    for rec,(a,b,wn) in zip(records,values):
        i,j,scale,d,point,n,_,rl,rr,vl,vr,pl,pr,ps,sm=rec
        flux=scale*face_flux(point,n,vl,vr,a,b,pl/((5/3-1)*rl),pr/((5/3-1)*rr),ps,sm)
        dr[i]-=flux;traction=scale*ps*n;F[i]-=traction;H[i]-=np.outer(point-states[i]['c'],traction)
        if j>=0:dr[j]+=flux;F[j]+=traction;H[j]+=np.outer(point-states[j]['c'],traction)
        else:ledger-=np.r_[flux[13:16],flux[34],np.cross(point,flux[13:16])]
        K[i]+=scale/mass[i]*(ALPHA.T@(a[:,None]*ALPHA))
        if j>=0:K[j]+=scale/mass[j]*(ALPHA.T@((-b)[:,None]*ALPHA))
        for receiver,v,rates in ((i,vr,-scale*b),(j,vl,scale*a)):
            if receiver<0:continue
            delta=v-np.array([velocity(receiver,point)])
            incoming[receiver].append((np.sqrt(rates)[:,None]*delta).T)
    for i,s in enumerate(states):
        m,q,Q,p,D,Pi,E=unpack(raw[i]);PV=s['P']*s['density'].volume
        dr[i,1:4]+=p;dr[i,4:13]+=(D+D.T).ravel();dr[i,16:25]+=(Pi+PV*np.eye(3)).ravel()
        H[i]+=PV*np.eye(3);A=np.linalg.solve(s['S'],H[i]).T;a0=F[i]/m-A@s['c']
        dr[i,25:34]+=(np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T).ravel()
    return dr,gen,z,dict(**gdiag,**pressure_diag,mass_flux=float(norm),ledger=ledger,weight_iterations=iteration,cov_linearization_factors=factors,cov_A=[-s['B'] for s in states],cov_K=K,cov_Qfactor=[np.concatenate(v,axis=1) if v else np.zeros((3,0)) for v in incoming],mass_rates=dr[:,0].copy(),mass_root_solved=audit_geometry_rates is None,generator_velocity_defect=float(np.max(abs(gen-material_gen))),pressure_model=pressure_reconstruction+'_exact_euler_star',face_integration=face_integration)
