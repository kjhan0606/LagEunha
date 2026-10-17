"""Positive discrete-velocity face transport plus compatible pressure force.
Four symmetric velocity samples represent each residual covariance. Each sample
is upwinded by its own relative normal speed. HLLC supplies pressure/work only.
This is an experimental moment closure, not a proved general Euler solver.
"""
import numpy as np
import os
import moment_ale_evolve as base

ALPHA=np.array([[np.sqrt(2),0.],[-np.sqrt(2),0.],[0.,np.sqrt(2)],[0.,-np.sqrt(2)]])

def pressure_gradients(states,edges):
    centers=np.array([s['density'].origin for s in states]);grad=np.zeros((len(states),2))
    neighbors=[set() for _ in states]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].add(j);neighbors[j].add(i)
    for i,nn in enumerate(neighbors):
        ids=list(nn);dr=centers[ids]-centers[i];d=np.linalg.norm(dr,axis=1)
        dp=np.array([states[j]['P']-states[i]['P'] for j in ids])
        g,_,rank,_=np.linalg.lstsq(dr/d[:,None],dp/d,rcond=None)
        if rank<2:raise RuntimeError('pressure stencil deficient')
        grad[i]=g
    return centers,grad


def transport(velL,velR,rhoL,rhoR,n,wn,epsilon=0.):
    sl=velL@n-wn;sr=velR@n-wn
    if epsilon:
        from kinetic_speed_smoothing import split
        lp,_,ld,_=split(sl,epsilon);_,rm,_,rd=split(sr,epsilon)
        return rhoL*.25*lp,rhoR*.25*rm,float(rhoL*.25*ld.sum()+rhoR*.25*rd.sum())
    left=rhoL*.25*np.maximum(sl,0.);right=rhoR*.25*np.minimum(sr,0.)
    derivative=rhoL*.25*np.where(sl>0,1.,np.where(sl<0,0.,.5)).sum()+rhoR*.25*np.where(sr<0,1.,np.where(sr>0,0.,.5)).sum()
    return left,right,float(derivative)


def evaluate_records(records,N,z):
    r=np.zeros(N);L=np.zeros((N,N));values=[]
    for rec in records:
        i,j,scale,wn0,d,vl,vr,rl,rr,n,*_=rec
        wn=wn0+(z[i]-z[j])/(2*d) if j>=0 else wn0
        epsilon=rec[15] if len(rec)>15 and np.ndim(rec[15])==0 else 0.
        a,b,rho=transport(vl,vr,rl,rr,n,wn,epsilon);fm=sum(a)+sum(b)
        values.append((a,b,wn));r[i]+=scale*fm
        if j>=0:
            r[j]-=scale*fm;c=scale*rho/(2*d)
            L[i,i]+=c;L[j,j]+=c;L[i,j]-=c;L[j,i]-=c
    return r,L,values


def solve(records,N):
    if os.environ.get('LAGFORCE_VECTOR_KINETIC_ROOT','0')=='1':
        from vector_kinetic_root import solve_records
        return solve_records(records,N)
    def evaluate(z):return evaluate_records(records,N,z)
    z=np.zeros(N)
    for iteration in range(50):
        r,L,values=evaluate(z);norm=float(max(abs(r)))
        if norm<2e-13:return z,values,norm
        # Vacuum intervals may leave several independent weight gauges.
        step=np.r_[np.linalg.lstsq(L[:-1,:-1],r[:-1],rcond=1e-13)[0],0.]
        for k in range(30):
            trial=z+step*2.**(-k)
            if max(abs(evaluate(trial)[0]))<norm:z=trial;break
        else:raise RuntimeError('kinetic mass constraint line search failed')
    raise RuntimeError('kinetic mass constraint did not converge')


def rhs(t,flat,mass,linear_pressure=True,covariance_factors=None,return_covariance_terms=True,
        audit_geometry_rates=None,face_density_provider=None,face_segments=1,face_order=4,exact_bin_transport=False,frozen_mask=None):
    # Optional rates ONLY evaluate the residual of an existing dense path.
    # Evolution callers leave this unset, retaining unsteered material sites
    # and the instantaneous mass solve. Do not use the audit mode to advance.
    if covariance_factors is None:raise ValueError('positive covariance factors required')
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic and (face_density_provider is not None or exact_bin_transport):
        raise ValueError('quadratic velocity currently requires continuous moment density and point quadrature')
    from kinetic_speed_smoothing import parameters,face_speed
    smoothing,smoothing_length=parameters()
    if smoothing and (face_density_provider is not None or exact_bin_transport):raise ValueError('smooth splitting requires continuous density quadrature')
    N=len(mass);y=flat[:-4].reshape(N,20);x=y[:,:2];w=y[:,2]
    frozen=np.zeros(N,dtype=bool) if frozen_mask is None else np.asarray(frozen_mask,dtype=bool)
    if frozen.shape!=(N,):raise ValueError('invalid frozen-cell mask')
    has_frozen=bool(np.any(frozen))
    if has_frozen and (audit_geometry_rates is not None or exact_bin_transport):
        raise ValueError('frozen boundary requires its explicit point-quadrature evolution solve')
    polys,edges=base.geometry(x,w,t);raw=np.column_stack([mass,y[:,3:]])
    states=[base.decode(p,r) for p,r in zip(polys,raw)]
    trace=None;density_diag={}
    if face_density_provider is not None:trace,density_diag=face_density_provider(polys,edges,raw)
    for s,L in zip(states,covariance_factors):s['R']=L@L.T
    if quadratic:
        from quadratic_velocity_reconstruction import prepare,velocity_samples
        density_diag.update(prepare(states,polys,edges,mass,covariance_factors))
    material_gen=np.array([base.velocity(s,p) for s,p in zip(states,x)])
    gen=material_gen.copy() if has_frozen else material_gen
    if has_frozen:gen[frozen]=0.
    if audit_geometry_rates is not None:
        gen=np.asarray(audit_geometry_rates[0]);prescribed_z=np.asarray(audit_geometry_rates[1])
        if gen.shape!=(N,2) or prescribed_z.shape!=(N,) or not np.all(np.isfinite(gen)) or not np.all(np.isfinite(prescribed_z)):
            raise ValueError('invalid diagnostic geometry rates')
    centers=np.array([s['density'].origin for s in states]);grad=np.zeros((N,2))
    if linear_pressure:
        centers,grad=pressure_gradients(states,edges)
    from quadratic_pressure_reconstruction import enabled as quadratic_pressure_enabled
    pressure_polynomials=None
    if linear_pressure and quadratic_pressure_enabled():
        from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value
        pressure_polynomials,pressure_diag=prepare_pressure(states,polys);density_diag.update(pressure_diag)
    def pressure(i,r):
        p=states[i]['P']+grad[i]@(r-centers[i]) if pressure_polynomials is None else float(pressure_value(pressure_polynomials[i],r)[0])
        if p<=0:raise RuntimeError('nonpositive reconstructed pressure')
        return p
    if face_segments<1 or face_order<1:raise ValueError('invalid face quadrature')
    records=[];gx,gw=np.polynomial.legendre.leggauss(face_order);gx=(gx+1)/2;gw/=2
    nodes=[((j+g)/face_segments,q/face_segments) for j in range(face_segments) for g,q in zip(gx,gw)]
    bins=[];analytic_mass=None;mass_iterations=None
    if exact_bin_transport:
        if trace is None:raise ValueError('exact-bin integration requires constant-bin density provider')
        from affine_kinetic_mass import solve as solve_bins,cuts,evaluate as evaluate_bins,arrays as bin_arrays
        from affine_kinetic_mass import make_bins
        bins=make_bins(edges,x,states,covariance_factors,mass,gen,trace,face_segments,base.boundary_velocity,base.velocity,ALPHA)
        if audit_geometry_rates is None:z,analytic_mass,mass_iterations=solve_bins(bins,N)
        else:
            z=prescribed_z;analytic_mass=float(max(abs(evaluate_bins(bin_arrays(bins),N,z)[0])))
    for edge_id,(i,j,a,b,n,length) in enumerate(edges):
        epsilon=face_speed(states,i,j,smoothing,smoothing_length)
        d=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        edge_nodes=nodes
        if exact_bin_transport:
            edge_nodes=[]
            for k in range(face_segments):
                split=cuts(bins[edge_id*face_segments+k],z)
                for lo,hi in zip(split[:-1],split[1:]):
                    edge_nodes.extend(((k+lo+(hi-lo)*g)/face_segments,(hi-lo)*q/face_segments) for g,q in zip(gx,gw))
        for g,q in edge_nodes:
            point=a+g*(b-a);vl=base.velocity(states[i],point);rl=float(states[i]['density'].value(point)[0]) if trace is None else trace(i,edge_id,g);pl=pressure(i,point)
            vls=vl+ALPHA@covariance_factors[i].T/np.sqrt(mass[i])
            if quadratic:
                samples,coordinates=velocity_samples(states[i],covariance_factors[i],mass[i],point,ALPHA)
                vls=samples[0];xil=coordinates[0]
            if j>=0:
                vr=base.velocity(states[j],point);rr=float(states[j]['density'].value(point)[0]) if trace is None else trace(j,edge_id,g);pr=pressure(j,point)
                vrs=vr+ALPHA@covariance_factors[j].T/np.sqrt(mass[j])
                if quadratic:
                    samples,coordinates=velocity_samples(states[j],covariance_factors[j],mass[j],point,ALPHA)
                    vrs=samples[0];xir=coordinates[0]
                wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
            else:
                wn=base.boundary_velocity(point)@n;vr=vl-2*(vl@n-wn)*n;rr=rl;pr=pl
                vrs=vls-2*((vls@n-wn)[:,None]*n)
                if quadratic:xir=ALPHA
            pressure_solver=os.environ.get('LAGFORCE_PRESSURE_RIEMANN','hllc_davis')
            if has_frozen and frozen[i] and (j<0 or frozen[j]):
                # Neither cell evolves. No physical pressure traction is needed;
                # keep kinetic records for the existing reservoir ledger algebra.
                ps=0.;sm=0.
            elif pressure_solver=='exact_euler_star':
                from euler_star_pressure import star
                ps,sm=star(rl,vl,pl,rr,vr,pr,n,base.GAMMA)
            elif pressure_solver=='hllc_davis':
                fan=base.fan(base.state(rl,vl,pl),base.state(rr,vr,pr),n)
                sm=fan[0][1];ps=float((fan[2][1][1:3]-fan[1][1][1:3]*sm)@n)
            else:raise ValueError('unknown pressure Riemann solver')
            record=(i,j,length*q,wn,d,vls,vrs,rl,rr,n,point,pl,pr,ps,sm,epsilon)
            records.append(record+(xil,xir) if quadratic else record)
    if exact_bin_transport:
        residual,_,values=evaluate_records(records,N,z);mass_residual=float(np.max(abs(residual)))
        if abs(mass_residual-analytic_mass)>1e-11:raise RuntimeError('split-bin mass integral mismatch')
    elif has_frozen:
        from frozen_weight_root import solve as solve_frozen
        z,values,mass_residual=solve_frozen(records,N,frozen)
    elif audit_geometry_rates is None:
        z,values,mass_residual=solve(records,N)
    else:
        z=prescribed_z
        residual,_,values=evaluate_records(records,N,z)
        mass_residual=float(np.max(abs(residual)))
    dr=np.zeros((N,18));ledger=np.zeros(4);F=np.zeros((N,2));H=np.zeros((N,2,2));K=np.zeros((N,2,2));incoming=[[] for _ in range(N)];area_rates=np.zeros(N)
    from vector_moment_flux import backend_name
    if backend_name()=='batched':
        from vector_moment_flux import accumulate
        dr,ledger,F,H,K,incoming,area_rates=accumulate(records,values,states,mass,base.GAMMA)
    else:
        for rec,(wl,wr,wn) in zip(records,values):
            i,j,scale,_,_,vl,vr,rl,rr,n,point,pl,pr,ps,sm,epsilon,*coordinates=rec
            area_rates[i]+=scale*wn
            if j>=0:area_rates[j]-=scale*wn
            fm=sum(wl)+sum(wr);fp=wl@vl+wr@vr+ps*n
            Pi=vl.T@(wl[:,None]*vl)+vr.T@(wr[:,None]*vr)
            fe=wl@(.5*np.sum(vl*vl,axis=1)+pl/((base.GAMMA-1)*rl))+wr@(.5*np.sum(vr*vr,axis=1)+pr/((base.GAMMA-1)*rr))+ps*sm
            ff=scale*np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,np.outer(point,fp).ravel(),Pi.ravel(),fe]
            dr[i]-=ff;traction=scale*ps*n;F[i]-=traction;H[i]-=np.outer(point-states[i]['c'],traction)
            if j>=0:dr[j]+=ff;F[j]+=traction;H[j]+=np.outer(point-states[j]['c'],traction)
            else:ledger-=scale*np.r_[fp,fe,point[0]*fp[1]-point[1]*fp[0]]
            xil,xir=coordinates if coordinates else (ALPHA,ALPHA)
            K[i]+=scale/mass[i]*(xil.T@(wl[:,None]*xil))
            if j>=0:K[j]+=scale/mass[j]*(xir.T@((-wr)[:,None]*xir))
            # Incoming samples are positive contributions to the receiver residual.
            for receiver,vel,rates in [(i,vr,-scale*wr),(j,vl,scale*wl)]:
                if receiver<0:continue
                delta=vel-base.affine_velocity(states[receiver],point)
                incoming[receiver].append((np.sqrt(rates)[:,None]*delta).T)
    for i,s in enumerate(states):
        r=raw[i];p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2);PV=s['P']*s['density'].area
        dr[i,1:3]+=p;dr[i,3:7]+=(D+D.T).ravel();dr[i,9:13]+=(Pi+PV*np.eye(2)).ravel()
        S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/mass[i];H[i]+=PV*np.eye(2)
        A=np.linalg.solve(S,H[i]).T;a0=F[i]/mass[i]-A@s['c']
        dr[i,13:17]+=(np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T).ravel()
    frozen_exchange=np.zeros(4)
    if has_frozen:
        removed=dr[frozen]
        frozen_exchange=-np.r_[removed[:,7:9].sum(axis=0),removed[:,17].sum(),(removed[:,10]-removed[:,11]).sum()]
        ledger+=frozen_exchange;dr[frozen]=0.
    diag=dict(cov_linearization_factors=covariance_factors,mass_flux=mass_residual,cov_A=[-s['B'] for s in states],cov_K=list(K),
        cov_Qfactor=[np.concatenate(a,axis=1) if a else np.zeros((2,0)) for a in incoming],
        min_pressure=min(s['P'] for s in states),mass_rates=dr[:,0].copy(),
        generator_velocity_defect=float(np.max(abs((gen-material_gen)[~frozen]),initial=0.)))
    if has_frozen:
        for i in np.flatnonzero(frozen):
            diag['cov_A'][i]=np.zeros((2,2));diag['cov_K'][i]=np.zeros((2,2));diag['cov_Qfactor'][i]=np.zeros((2,0))
        diag.update(frozen_cells=int(frozen.sum()),frozen_exchange_ledger=frozen_exchange)
    area=np.array([s['density'].area for s in states]);P=np.array([s['P'] for s in states])
    diag.update(density_diag,area_rates=area_rates,
        pressure_rates=(base.GAMMA-1)*(dr[:,17]-.5*(dr[:,13]+dr[:,16]))/area-P*area_rates/area)
    if exact_bin_transport:diag.update(analytic_mass_flux=analytic_mass,mass_iterations=mass_iterations,face_quadrature_records=len(records))
    if not exact_bin_transport:diag['mass_root_backend']='vector_sparse' if os.environ.get('LAGFORCE_VECTOR_KINETIC_ROOT','0')=='1' else 'legacy_dense'
    if has_frozen:diag['mass_root_backend']='frozen_dirichlet_sparse'
    if quadratic:
        from quadratic_velocity_reconstruction import closure_name
        diag['velocity_reconstruction']=closure_name()
    else:diag['velocity_reconstruction']='affine'
    diag.update(kinetic_smoothing_coefficient=smoothing,kinetic_smoothing_length=smoothing_length,
                pressure_riemann_solver=os.environ.get('LAGFORCE_PRESSURE_RIEMANN','hllc_davis'))
    if 'pressure_reconstruction' not in diag:diag['pressure_reconstruction']='linear' if linear_pressure else 'constant'
    return np.r_[np.column_stack([gen,z,dr[:,1:]]).ravel(),ledger],diag
