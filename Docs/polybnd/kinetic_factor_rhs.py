"""Positive discrete-velocity face transport plus compatible pressure force.
Four symmetric velocity samples represent each residual covariance. Each sample
is upwinded by its own relative normal speed. HLLC supplies pressure/work only.
This is an experimental moment closure, not a proved general Euler solver.
"""
import numpy as np
import os
import moment_ale_evolve as base

ALPHA=np.array([[np.sqrt(2),0.],[-np.sqrt(2),0.],[0.,np.sqrt(2)],[0.,-np.sqrt(2)]])

def angular_rate_partition(raw, rates, mass, pressure_centered_moment, frozen):
    """Hydro-only internal angular rate, including moving-center corrections.

    pressure_centered_moment includes the isotropic volume pressure term,
    whose antisymmetric part vanishes. The remainder is transport; gravity
    and any later factor rebuild are outside this diagnostic.
    """
    q=raw[:,1:3];p=raw[:,7:9]
    dq=rates[:,1:3];dp=rates[:,7:9];dm=rates[:,0]
    cross=lambda a,b:a[:,0]*b[:,1]-a[:,1]*b[:,0]
    internal=rates[:,10]-rates[:,11]-(cross(dq,p)+cross(q,dp))/mass+cross(q,p)*dm/mass**2
    pressure=pressure_centered_moment[:,0,1]-pressure_centered_moment[:,1,0]
    pressure=pressure.copy();pressure[frozen]=0.
    return dict(internal_angular_rate=internal,pressure_internal_torque=pressure,
                transport_internal_angular_rate=internal-pressure)

def angular_face_partition(records, values, raw, mass, selected):
    """Resolve instantaneous transport into sample-mean and residual flux.

    Sample mean is the equally weighted velocity ensemble mean at the face,
    not the cell bulk velocity. This is a diagnostic, not a flux replacement.
    """
    rows={i:{} for i in selected}
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    for rec,(wl,wr,_) in zip(records,values):
        i,j,scale,_,_,vl,vr,_,_,_,point,*_=rec
        if i not in rows and j not in rows:continue
        fm=sum(wl)+sum(wr)
        adv=wl@vl+wr@vr
        mean=sum(wl)*np.mean(vl,axis=0)+sum(wr)*np.mean(vr,axis=0)
        for cell,neighbor,sign in ((i,j,-1.),(j,i,1.)):
            if cell not in rows:continue
            c=raw[cell,1:3]/mass[cell];v=raw[cell,7:9]/mass[cell]
            lever=point-c
            item=rows[cell].setdefault(neighbor,dict(mean=0.,residual=0.,total=0.,absolute_quadrature_sum=0.))
            total=sign*scale*cross(lever,adv-v*fm)
            item['mean']+=sign*scale*cross(lever,mean-v*fm)
            item['residual']+=sign*scale*cross(lever,adv-mean)
            item['total']+=total;item['absolute_quadrature_sum']+=abs(total)
    return rows

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


def prepare_face_state_batch(edges,nodes,face_node_provider,states,gen,factors,mass):
    """One outer evaluation's cell batches, scattered to original node order.

    Face splitting is performed once by the unchanged provider. The cache
    belongs to this call only. Exterior right states remain mirror states
    constructed by the existing RHS, not independent density evaluations.
    """
    from quadratic_velocity_reconstruction import velocity_samples_and_mean_batch
    edge_nodes=[];offsets=[0];points=[];normals=[];owners=[]
    entries=[[] for _ in states]
    for edge_id,edge in enumerate(edges):
        i,j,a,b,n,length=edge
        chosen=nodes if face_node_provider is None else face_node_provider(edge_id,edge,states,gen)
        if face_node_provider is not None:
            if not chosen or any(not (0<=g<=1 and q>0 and np.isfinite(q)) for g,q in chosen) or abs(sum(q for g,q in chosen)-1)>1e-12:
                raise ValueError('invalid custom face quadrature')
        edge_nodes.append(chosen)
        for g,q in chosen:
            row=len(points)
            points.append(a+g*(b-a));normals.append(n);owners.append((i,j))
            entries[i].append((row,0))
            if j>=0:entries[j].append((row,1))
        offsets.append(len(points))
    points=np.asarray(points,dtype=float).reshape(-1,2)
    normals=np.asarray(normals,dtype=float).reshape(-1,2)
    owners=np.asarray(owners,dtype=np.int64).reshape(-1,2)
    count=len(points)
    rho=np.full((count,2),np.nan)
    samples=np.full((count,2,len(ALPHA),2),np.nan)
    coordinates=np.full_like(samples,np.nan)
    means=np.full((count,2,2),np.nan)
    sizes=[];branches=np.zeros(3,dtype=np.int64)
    for cell,items in enumerate(entries):
        if not items:continue
        indices=np.asarray(items,dtype=np.int64);rows=indices[:,0];sides=indices[:,1]
        cell_points=points[rows]
        inward=(2*sides-1)[:,None]*normals[rows]
        density=states[cell]['density']
        if not hasattr(density,'value_trace_batch'):
            raise ValueError('batch face state requires an explicit batch density trace')
        values,branch=density.value_trace_batch(cell_points,inward,return_branch=True)
        v,xi,mean=velocity_samples_and_mean_batch(states[cell],factors[cell],mass[cell],cell_points,ALPHA)
        rho[rows,sides]=values;samples[rows,sides]=v;coordinates[rows,sides]=xi;means[rows,sides]=mean
        sizes.append(len(rows));branches+=np.bincount(branch+1,minlength=3)
    offsets=np.asarray(offsets,dtype=np.int64)
    arrays=(points,normals,owners,rho,samples,coordinates,means,offsets)
    stats=dict(nodes=count,cell_batches=len(sizes),point_sides=sum(sizes),
        maximum_cell_batch=max(sizes,default=0),minimum_cell_batch=min(sizes,default=0),
        branch_counts_low_flat_high=branches.tolist(),array_bytes=sum(a.nbytes for a in arrays),
        scalar_face_density_calls=0,scalar_face_sample_calls=0,
        cache_scope='this outer RHS evaluation only')
    return dict(points=points,normals=normals,owners=owners,rho=rho,samples=samples,
                coordinates=coordinates,means=means,edge_nodes=edge_nodes,offsets=offsets,stats=stats)


def rhs(t,flat,mass,linear_pressure=True,covariance_factors=None,return_covariance_terms=True,
        audit_geometry_rates=None,face_density_provider=None,face_segments=1,face_order=4,exact_bin_transport=False,frozen_mask=None,density_provider=None,face_node_provider=None,
        audit_face_density=None,audit_weight_rates=None,audit_capture_nodes=False):
    # Optional rates ONLY evaluate the residual of an existing dense path.
    # Evolution callers leave this unset, retaining unsteered material sites
    # and the instantaneous mass solve. Do not use the audit mode to advance.
    if covariance_factors is None:raise ValueError('positive covariance factors required')
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    shared_velocity=quadratic and os.environ.get('LAGFORCE_SHARED_FACE_VELOCITY','0')=='1'
    if shared_velocity:
        from quadratic_velocity_reconstruction import velocity_samples_and_mean
    if quadratic and (face_density_provider is not None or exact_bin_transport):
        raise ValueError('quadratic velocity currently requires continuous moment density and point quadrature')
    from kinetic_speed_smoothing import parameters,face_speed
    smoothing,smoothing_length=parameters()
    if face_node_provider is not None:
        if smoothing or exact_bin_transport:raise ValueError('custom face nodes require unsmoothed point transport')
        if quadratic and getattr(face_node_provider,'velocity_degree',1)!=2:
            raise ValueError('quadratic velocity requires a quadratic face-node provider')
    if smoothing and (face_density_provider is not None or exact_bin_transport):raise ValueError('smooth splitting requires continuous density quadrature')
    N=len(mass);y=flat[:-4].reshape(N,20);x=y[:,:2];w=y[:,2]
    frozen=np.zeros(N,dtype=bool) if frozen_mask is None else np.asarray(frozen_mask,dtype=bool)
    if frozen.shape!=(N,):raise ValueError('invalid frozen-cell mask')
    has_frozen=bool(np.any(frozen))
    if has_frozen and (audit_geometry_rates is not None or exact_bin_transport):
        raise ValueError('frozen boundary requires its explicit point-quadrature evolution solve')
    polys,edges=base.geometry(x,w,t);raw=np.column_stack([mass,y[:,3:]])
    if density_provider is not None and (face_density_provider is not None or exact_bin_transport):
        raise ValueError('provided volume density requires matching continuous face density')
    states=[base.decode(p,r,density=None if density_provider is None else density_provider(i,p,r)) for i,(p,r) in enumerate(zip(polys,raw))]
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
    batch_face=os.environ.get('LAGFORCE_BATCH_FACE_STATE','0')=='1'
    face_cache=None
    if batch_face:
        if not shared_velocity or trace is not None or exact_bin_transport or smoothing or audit_face_density is not None:
            raise ValueError('batch face states require unsmoothed shared quadratic point states without density overrides')
        face_cache=prepare_face_state_batch(edges,nodes,face_node_provider,states,gen,covariance_factors,mass)
    import time
    record_begin=time.perf_counter();record_cpu=time.process_time()
    record_batch_stats=None;record_batch_audit=None
    if os.environ.get('LAGFORCE_BATCH_FACE_RECORDS','0')=='1':
        if face_cache is None or pressure_polynomials is not None or exact_bin_transport or audit_face_density is not None or smoothing or not quadratic or os.environ.get('LAGFORCE_PRESSURE_RIEMANN','hllc_davis')!='exact_euler_star':
            raise ValueError('native records require approved shared quadratic/cached scalar-pressure/exact-star contract')
        from face_record_batch import build as build_records
        records,record_batch_stats,record_batch_audit=build_records(face_cache,edges,states,x,gen,frozen,centers,grad,base.boundary_velocity,base.GAMMA,ALPHA)
    else:
        for edge_id,(i,j,a,b,n,length) in enumerate(edges):
            epsilon=face_speed(states,i,j,smoothing,smoothing_length)
            d=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
            edge_nodes=nodes
            if face_cache is not None:
                edge_nodes=face_cache['edge_nodes'][edge_id]
            elif face_node_provider is not None:
                edge_nodes=face_node_provider(edge_id,(i,j,a,b,n,length),states,gen)
                if not edge_nodes or any(not (0<=g<=1 and q>0 and np.isfinite(q)) for g,q in edge_nodes) or abs(sum(q for g,q in edge_nodes)-1)>1e-12:
                    raise ValueError('invalid custom face quadrature')
            if exact_bin_transport:
                edge_nodes=[]
                for k in range(face_segments):
                    split=cuts(bins[edge_id*face_segments+k],z)
                    for lo,hi in zip(split[:-1],split[1:]):
                        edge_nodes.extend(((k+lo+(hi-lo)*g)/face_segments,(hi-lo)*q/face_segments) for g,q in zip(gx,gw))
            for node_id,(g,q) in enumerate(edge_nodes):
                cached_row=None if face_cache is None else face_cache['offsets'][edge_id]+node_id
                point=a+g*(b-a) if face_cache is None else face_cache['points'][cached_row]
                dl=states[i]['density']
                rl=float(face_cache['rho'][cached_row,0]) if face_cache is not None else (float(dl.value_trace(point,-n)[0] if hasattr(dl,'value_trace') else dl.value(point)[0]) if trace is None else trace(i,edge_id,g))
                pl=pressure(i,point)
                if face_cache is not None:
                    vls=face_cache['samples'][cached_row,0];xil=face_cache['coordinates'][cached_row,0];vl=face_cache['means'][cached_row,0]
                elif shared_velocity:
                    samples,coordinates,vl=velocity_samples_and_mean(states[i],covariance_factors[i],mass[i],point,ALPHA)
                    vls=samples[0];xil=coordinates[0]
                else:
                    vl=base.velocity(states[i],point)
                    vls=vl+ALPHA@covariance_factors[i].T/np.sqrt(mass[i])
                if quadratic and not shared_velocity:
                    samples,coordinates=velocity_samples(states[i],covariance_factors[i],mass[i],point,ALPHA)
                    vls=samples[0];xil=coordinates[0]
                if j>=0:
                    dright=states[j]['density']
                    rr=float(face_cache['rho'][cached_row,1]) if face_cache is not None else (float(dright.value_trace(point,n)[0] if hasattr(dright,'value_trace') else dright.value(point)[0]) if trace is None else trace(j,edge_id,g))
                    pr=pressure(j,point)
                    if face_cache is not None:
                        vrs=face_cache['samples'][cached_row,1];xir=face_cache['coordinates'][cached_row,1];vr=face_cache['means'][cached_row,1]
                    elif shared_velocity:
                        samples,coordinates,vr=velocity_samples_and_mean(states[j],covariance_factors[j],mass[j],point,ALPHA)
                        vrs=samples[0];xir=coordinates[0]
                    else:
                        vr=base.velocity(states[j],point)
                        vrs=vr+ALPHA@covariance_factors[j].T/np.sqrt(mass[j])
                    if quadratic and not shared_velocity:
                        samples,coordinates=velocity_samples(states[j],covariance_factors[j],mass[j],point,ALPHA)
                        vrs=samples[0];xir=coordinates[0]
                    wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
                else:
                    wn=base.boundary_velocity(point)@n;vr=vl-2*(vl@n-wn)*n;rr=rl;pr=pl
                    vrs=vls-2*((vls@n-wn)[:,None]*n)
                    if quadratic:xir=ALPHA
                # Counterfactual audit only: preserve decoded density, velocity
                # basis, generator and volume moments; replace face traces only.
                if audit_face_density is not None:
                    rl,rr=map(float,audit_face_density(point,i,j))
                    if not np.isfinite(rl+rr) or min(rl,rr)<=0:
                        raise ValueError('invalid audit face density')
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
    record_timing=dict(wall_seconds=time.perf_counter()-record_begin,process_cpu_seconds=time.process_time()-record_cpu)
    if audit_weight_rates is not None:
        if exact_bin_transport:raise ValueError('fixed audit weights incompatible with exact bins')
        z=np.asarray(audit_weight_rates,dtype=float)
        if z.shape!=(N,) or not np.isfinite(z).all() or np.any(z[frozen]!=0):
            raise ValueError('invalid fixed audit weight rates')
        residual,_,values=evaluate_records(records,N,z)
        mass_residual=float(np.max(abs(residual[~frozen])))
    elif exact_bin_transport:
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
    accumulation_begin=time.perf_counter();accumulation_cpu=time.process_time()
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
    accumulation_timing=dict(wall_seconds=time.perf_counter()-accumulation_begin,process_cpu_seconds=time.process_time()-accumulation_cpu)
    capture_pressure=os.environ.get('LAGFORCE_CAPTURE_PRESSURE_RATES','0')=='1'
    pressure_raw=np.zeros((N,18)) if capture_pressure else None
    for i,s in enumerate(states):
        r=raw[i];p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2);PV=s['P']*s['density'].area
        dr[i,1:3]+=p;dr[i,3:7]+=(D+D.T).ravel();dr[i,9:13]+=(Pi+PV*np.eye(2)).ravel()
        S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/mass[i];H[i]+=PV*np.eye(2)
        A=np.linalg.solve(S,H[i]).T;a0=F[i]/mass[i]-A@s['c']
        force_Pi=np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T
        dr[i,13:17]+=force_Pi.ravel()
        if capture_pressure:
            pressure_raw[i,7:9]=F[i]
            pressure_raw[i,9:13]=(H[i]+np.outer(s['c'],F[i])).ravel()
            pressure_raw[i,13:17]=force_Pi.ravel()
    if capture_pressure:
        for rec in records:
            i,j,scale=rec[:3];work=scale*rec[13]*rec[14]
            pressure_raw[i,17]-=work
            if j>=0:pressure_raw[j,17]+=work
        pressure_raw[frozen]=0.
    frozen_exchange=np.zeros(4)
    if has_frozen:
        removed=dr[frozen]
        frozen_exchange=-np.r_[removed[:,7:9].sum(axis=0),removed[:,17].sum(),(removed[:,10]-removed[:,11]).sum()]
        ledger+=frozen_exchange;dr[frozen]=0.
    diag=dict(cov_linearization_factors=covariance_factors,mass_flux=mass_residual,cov_A=[-s['B'] for s in states],cov_K=list(K),
        cov_Qfactor=[np.concatenate(a,axis=1) if a else np.zeros((2,0)) for a in incoming],
        min_pressure=min(s['P'] for s in states),mass_rates=dr[:,0].copy(),
        generator_velocity_defect=float(np.max(abs((gen-material_gen)[~frozen]),initial=0.)))
    diag.update(angular_rate_partition(raw,dr,mass,H,frozen))
    if face_cache is not None:diag['batch_face_state']=face_cache['stats']
    if audit_capture_nodes:
        if audit_capture_nodes=='records':
            diag['audit_face_records']=(records,values,states,gen.copy(),z.copy())
        else:
            from kinetic_face_capture import capture
            diag['captured_face_nodes']=capture(records,values,base.GAMMA)
    if capture_pressure:
        diag.update(pressure_raw_rates=pressure_raw,transport_kinematic_raw_rates=dr-pressure_raw)
    diag.update(audit_face_density_override=audit_face_density is not None,
                audit_fixed_weight_rates=audit_weight_rates is not None)
    selected=os.environ.get('LAGFORCE_ANGULAR_FACE_CELLS','')
    if selected:
        selected=[int(k) for k in selected.split(',')]
        if any(k<0 or k>=N or frozen[k] for k in selected):
            raise ValueError('angular face diagnostic requires active local cell indices')
        diag['angular_face_partition']=angular_face_partition(records,values,raw,mass,selected)
    if has_frozen:
        for i in np.flatnonzero(frozen):
            diag['cov_A'][i]=np.zeros((2,2));diag['cov_K'][i]=np.zeros((2,2));diag['cov_Qfactor'][i]=np.zeros((2,0))
        diag.update(frozen_cells=int(frozen.sum()),frozen_exchange_ledger=frozen_exchange)
    area=np.array([s['density'].area for s in states]);P=np.array([s['P'] for s in states])
    diag.update(face_record_total_timing=record_timing,moment_accumulation_timing=accumulation_timing,face_record_batch=record_batch_stats)
    if record_batch_audit is not None and os.environ.get('LAGFORCE_AUDIT_RECORD_BATCH','0')=='1':diag['face_record_batch_audit']=record_batch_audit
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
