"""Static affine reconstruction / moment-ALE consistency gates, not an evolver.
Uniform density within cells. Raw state: m,q,Q,p,D,Pi,E (18 entries).
"""
import json
import numpy as np
from weighted_mass_patch import power_cells
from periodic2d_hllc import quadrature,mesh as periodic_mesh,initial as periodic_initial
from fixedmass_hllc import state,fan,sample,solve_weights,GAMMA

DOMAIN=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
P0=1/GAMMA


def bounded_mesh(x,w):
    cells=power_cells(x,w,DOMAIN);polygons=[np.array(c.exterior.coords[:-1]) for c in cells];edges=[]
    for i,p in enumerate(polygons):
        if cells[i].area<=0:raise RuntimeError('empty cell')
        delta=x-x[i];bound=np.sum(delta*delta,axis=1)+w[i]-w
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-14:continue
            normal=np.array([b[1]-a[1],a[0]-b[0]])/length
            residual=abs(2*delta@((a+b)/2-x[i])-bound);residual[i]=np.inf
            j=int(np.argmin(residual))
            if residual[j]>1e-10:j=-1
            if j<0 or i<j:edges.append((i,j,a,b,normal,length))
    return polygons,edges


def integrate(poly,velocity):
    points,weights=quadrature(poly);area=weights.sum();c=weights@points/area
    v=velocity(points);mean=weights@v/area;dx=points-c;dv=v-mean
    S=(dx*weights[:,None]).T@dx;C=(dx*weights[:,None]).T@dv;T=(dv*weights[:,None]).T@dv
    Q=(points*weights[:,None]).T@points;D=(points*weights[:,None]).T@v;Pi=(v*weights[:,None]).T@v
    E=P0*area/(GAMMA-1)+.5*np.trace(Pi)
    raw=np.r_[area,area*c,Q.ravel(),area*mean,D.ravel(),Pi.ravel(),E]
    return dict(area=area,center=c,mean=mean,S=S,C=C,T=T,raw=raw)


def least_squares(centers,means,neighbors):
    gradients=[]
    for i,stencil in enumerate(neighbors):
        displacement=np.array([centers[j]+shift-centers[i] for j,shift in stencil])
        difference=np.array([means[j]-means[i] for j,shift in stencil])
        scale=np.linalg.norm(displacement,axis=1)
        fit,_,rank,_=np.linalg.lstsq(displacement/scale[:,None],difference/scale[:,None],rcond=None)
        if rank!=2:raise RuntimeError('rank-deficient velocity stencil')
        gradients.append(fit.T)
    return np.array(gradients)


def affine_case(n=4, capture=None):
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([xx.ravel(),yy.ravel()])+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    w=np.zeros(n*n);B=np.array([[0.,1.1],[0.,0.]]);b=np.array([.2,0.])
    velocity=lambda points:points@B.T+b
    polys,edges=bounded_mesh(x,w);mom=[integrate(poly,velocity) for poly in polys]
    centers=np.array([s['center'] for s in mom]);means=np.array([s['mean'] for s in mom])
    neighbors=[[] for _ in x]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].append((j,np.zeros(2)));neighbors[j].append((i,np.zeros(2)))
    gradients=least_squares(centers,means,neighbors)
    wrong_gradients=least_squares(x,means,neighbors)
    generator_velocity=means+np.einsum('nij,nj->ni',gradients,x-centers)
    pressures=[];naive=[];gradient_energy=[];face_mismatch=0.;records=[];face_data=[]
    for s,g in zip(mom,gradients):
        bulk=.5*s['area']*np.dot(s['mean'],s['mean']);kg=.5*np.trace(g@s['S']@g.T)
        gradient_energy.append(kg);pressures.append((GAMMA-1)*(s['raw'][-1]-bulk-kg)/s['area'])
        naive.append((GAMMA-1)*(s['raw'][-1]-bulk)/s['area'])
    gx,gw=np.polynomial.legendre.leggauss(4);gx=(gx+1)/2;gw=gw/2
    for i,j,a,bb,normal,length in edges:
        distance=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        for t,weight in zip(gx,gw):
            point=a+t*(bb-a);vi=means[i]+gradients[i]@(point-centers[i])
            vj=means[j]+gradients[j]@(point-centers[j]) if j>=0 else vi
            pj=pressures[j] if j>=0 else pressures[i]
            if j>=0:face_mismatch=max(face_mismatch,float(np.max(abs(vi-vj))))
            riemann=fan(state(1.,vi,pressures[i]),state(1.,vj,pj),normal)
            base=((x[j]-x[i])@generator_velocity[j]-(generator_velocity[j]-generator_velocity[i])@(point-x[i]))/distance if j>=0 else 0.
            records.append((i,j,length*weight,base,distance,riemann))
            face_data.append((point,normal,(vi+vj)/2,(pressures[i]+pj)/2))
    z,samples,iterations,residual=solve_weights(records)
    count=len(x);derivative=np.zeros((count,18));area_rate=np.zeros(count)
    for record,sampled,data in zip(records,samples,face_data):
        i,j,length,base,distance,_=record;flux,branch,wn=sampled;point,normal,v,pressure=data
        fm=flux[0];fp=flux[1:3];fpi=np.outer(v,v)*fm+pressure*(np.outer(v,normal)+np.outer(normal,v))
        f=np.r_[fm,point*fm,(np.outer(point,point)*fm).ravel(),fp,np.outer(point,fp).ravel(),fpi.ravel(),flux[3]]*length
        derivative[i]-=f;area_rate[i]+=length*wn
        if j>=0:derivative[j]+=f;area_rate[j]-=length*wn
    for i,s in enumerate(mom):
        raw=s['raw'];D=raw[9:13].reshape(2,2);Pi=raw[13:17].reshape(2,2)
        derivative[i,1:3]+=raw[7:9]
        derivative[i,3:7]+=(D+D.T).ravel()
        derivative[i,9:13]+=(Pi+pressures[i]*s['area']*np.eye(2)).ravel()
        derivative[i,13:17]+=(pressures[i]*s['area']*(gradients[i]+gradients[i].T)).ravel()
    errors=[]
    if capture is not None:
        capture.update(raw=np.array([s['raw'] for s in mom]),rate=derivative.copy(),finite_states=[])
    for eps in (1e-4,1e-5,1e-6):
        plus=bounded_mesh(x+eps*generator_velocity,w+eps*z)[0]
        minus=bounded_mesh(x-eps*generator_velocity,w-eps*z)[0]
        raw_plus=np.array([integrate(p,velocity)['raw'] for p in plus])
        raw_minus=np.array([integrate(p,velocity)['raw'] for p in minus])
        fd=(raw_plus-raw_minus)/(2*eps)
        if capture is not None:
            capture['finite_states'].append((eps,raw_plus,raw_minus))
        errors.append(dict(eps=eps,max_raw_moment_rhs_error=float(np.max(abs(fd-derivative)))))
    kinetic_rate=.5*(derivative[:,13]+derivative[:,16])
    pressure_rate=(GAMMA-1)*(derivative[:,-1]-kinetic_rate)/np.array([s['area'] for s in mom])-np.array(pressures)*area_rate/np.array([s['area'] for s in mom])
    return dict(n=n,gradient_error=float(np.max(abs(gradients-B))),wrong_origin_gradient_error=float(np.max(abs(wrong_gradients-B))),
        face_trace_mismatch=face_mismatch,generator_point_velocity_error=float(np.max(abs(generator_velocity-velocity(x)))),
        point_vs_mean_velocity_difference=float(np.max(abs(generator_velocity-means))),
        split_pressure_error=float(np.max(abs(np.array(pressures)-P0))),naive_pressure_bias=float(np.max(np.array(naive)-P0)),
        total_gradient_energy=float(sum(gradient_energy)),mass_flux_residual=residual,area_rate_max=float(np.max(abs(area_rate))),
        transverse_momentum_rate_max=float(np.max(abs(derivative[:,8]))),pressure_rate_max=float(np.max(abs(pressure_rate))),
        total_energy_rate=float(derivative[:,-1].sum()),total_momentum_rate=derivative[:,7:9].sum(0).tolist(),
        angular_momentum_rate=float(np.sum(derivative[:,10]-derivative[:,11])),finite_difference=errors)


def curved_case(n=8):
    mass,y=periodic_initial(n,'uniform',np.zeros(2));x=y[:,:2];polys,area,edges=periodic_mesh(x,y[:,2])
    velocity=lambda points:np.column_stack([.5*np.sin(2*np.pi*points[:,1]),np.zeros(len(points))])
    moments=[integrate(poly,velocity) for poly in polys];centers=np.array([s['center'] for s in moments]);means=np.array([s['mean'] for s in moments])
    neighbors=[[] for _ in x]
    for i,j,a,b,delta,length in edges:
        shift=delta-((x[j]%1)-(x[i]%1));neighbors[i].append((j,shift));neighbors[j].append((i,-shift))
    ls=least_squares(centers,means,neighbors);ls_eigs=[];moment_eigs=[];pressure_errors=[];critical=[];spin_errors=[]
    actual_residual_eigs=[];cross_identity_errors=[];cross_energies=[]
    max_gradient_difference=0.
    for s,g in zip(moments,ls):
        B=np.linalg.solve(s['S'],s['C']).T;R=s['T']-B@s['S']@B.T
        # This additive remainder is NOT the actual residual covariance for
        # arbitrary least-squares g: the velocity residual has a cross moment.
        Rls=s['T']-g@s['S']@g.T
        Cres=s['C']-s['S']@g.T
        Ractual=s['T']-g@s['C']-s['C'].T@g.T+g@s['S']@g.T
        actual_residual_eigs.append(float(np.linalg.eigvalsh((Ractual+Ractual.T)/2).min()))
        cross_identity_errors.append(float(np.max(abs(s['T']-(g@s['S']@g.T+g@Cres+Cres.T@g.T+Ractual)))))
        cross_energies.append(float(np.trace(g@Cres)))
        ls_eigs.append(float(np.linalg.eigvalsh((Rls+Rls.T)/2).min()));moment_eigs.append(float(np.linalg.eigvalsh((R+R.T)/2).min()))
        max_gradient_difference=max(max_gradient_difference,float(np.max(abs(g-B))))
        bulk=.5*s['area']*np.dot(s['mean'],s['mean']);kgrad=.5*np.trace(B@s['S']@B.T);kr=.5*np.trace(R)
        pressure_errors.append(abs((GAMMA-1)*(s['raw'][-1]-bulk-kgrad-kr)/s['area']-P0))
        critical.append((GAMMA-1)*max(.5*np.trace(g@s['S']@g.T-s['T']),0)/s['area'])
        spin_errors.append(abs((s['S']@B.T)[0,1]-(s['S']@B.T)[1,0]-(s['C'][0,1]-s['C'][1,0])))
    return dict(n=n,ls_negative_additive_remainder_cells=sum(e<-1e-14 for e in ls_eigs),min_ls_additive_remainder_eigenvalue=min(ls_eigs),
        min_actual_ls_residual_covariance_eigenvalue=min(actual_residual_eigs),cross_term_identity_error=max(cross_identity_errors),
        min_ls_cross_energy=min(cross_energies),max_ls_cross_energy=max(cross_energies),
        min_moment_covariance_eigenvalue=min(moment_eigs),split_pressure_error=max(pressure_errors),spin_error=max(spin_errors),
        max_cold_pressure_threshold=max(critical),max_ls_vs_moment_gradient_difference=max_gradient_difference)


if __name__=='__main__':
    for n in (4,8):print(json.dumps(dict(kind='affine',**affine_case(n))),flush=True)
    for n in (4,8,16):print(json.dumps(dict(kind='curved',**curved_case(n))),flush=True)
