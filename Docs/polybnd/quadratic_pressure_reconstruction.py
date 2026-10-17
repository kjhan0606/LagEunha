"""Conservative cell-mean quadratic pressure from numerical thermal energies.

The polynomial has exactly the existing geometric cell mean. A continuous
weighted local fit uses neighbor cell means and geometric second moments.
Bernstein bounds limit its variation only when needed for positive pressure;
the stored energy, mean pressure, mass and generator velocity are not changed.
"""
import os
import numpy as np


def enabled():return os.environ.get('LAGFORCE_QUADRATIC_PRESSURE','0')=='1'


def closure_name():
    if not enabled():return 'linear'
    degree=int(os.environ.get('LAGFORCE_PRESSURE_DEGREE','2'))
    if degree==2:return 'quadratic_cell_mean'
    if degree in (3,4):return 'polynomial_cell_mean_degree'+str(degree)
    raise ValueError('unsupported pressure polynomial degree')


def polygon_statistics(poly):
    poly=np.asarray(poly);origin=poly.mean(0);local=poly-origin
    a=local;b=np.roll(local,-1,axis=0);areas=.5*abs(a[:,0]*b[:,1]-a[:,1]*b[:,0]);area=areas.sum()
    if area<=0:raise ValueError('nonpositive pressure cell area')
    mean=np.einsum('n,ni->i',areas,(a+b)/3)/area
    second=np.einsum('n,nij->ij',areas,(np.einsum('ni,nj->nij',a,a)+np.einsum('ni,nj->nij',b,b))/6+(np.einsum('ni,nj->nij',a,b)+np.einsum('ni,nj->nij',b,a))/12)/area
    return area,origin+mean,second-np.outer(mean,mean)


def design(z,second):
    z=np.atleast_2d(z)
    return np.column_stack([z,z[:,0]**2-second[0,0],z[:,0]*z[:,1]-second[0,1],z[:,1]**2-second[1,1]])


def value(metadata,points):
    if metadata.get('degree',2)>2:
        from polynomial_pressure_reconstruction import value as polynomial_value
        return polynomial_value(metadata,points)
    z=(np.atleast_2d(points)-metadata['center'])/metadata['scale']
    return metadata['mean']+design(z,metadata['second'])@metadata['coefficient']


def lower_bernstein_bound(metadata,poly):
    poly=np.asarray(poly);center=poly.mean(0);lower=float('inf')
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        tri=np.array([center,a,b]);v=value(metadata,tri);lower=min(lower,float(v.min()))
        for i,j in [(0,1),(0,2),(1,2)]:
            mixed=2*value(metadata,[(tri[i]+tri[j])/2])[0]-.5*(v[i]+v[j])
            lower=min(lower,float(mixed))
    return lower


def bernstein_lower_certificate(metadata,poly):
    # For a vertex/midpoint basis envelope M, every mixed Bernstein
    # coefficient has magnitude <= 2*M+.5*(M+M)=3*M. If this cheap
    # enclosure is positive, the original detailed test cannot limit.
    z=(np.asarray(poly)-metadata['center'])/metadata['scale'];x,y=np.max(abs(z),axis=0);S=metadata['second']
    envelope=np.array([x,y,x*x+abs(S[0,0]),x*y+abs(S[0,1]),y*y+abs(S[1,1])])
    variation=3*float(abs(metadata['coefficient'])@envelope)
    return metadata['mean']-variation-64*np.finfo(float).eps*(abs(metadata['mean'])+variation)


def prepare(states,cells):
    extension=os.environ.get('LAGFORCE_PRESSURE_REFLECT_Y','')
    if not extension:return prepare_extended(states,cells)
    bounds=np.asarray([float(v) for v in extension.split(',')])
    if bounds.shape!=(2,) or not np.isfinite(bounds).all() or bounds[0]>=bounds[1]:raise ValueError('invalid pressure reflection bounds')
    points=np.vstack(cells)
    if points[:,1].min()<bounds[0]-1e-12 or points[:,1].max()>bounds[1]+1e-12:raise ValueError('pressure cells outside reflection domain')
    extended_cells=list(cells);extended_states=list(states)
    for wall in bounds:
        for s,p in zip(states,cells):
            ghost=np.array(p,copy=True);ghost[:,1]=2*wall-ghost[:,1]
            extended_cells.append(ghost);extended_states.append(s)
    result,diag=prepare_extended(extended_states,extended_cells,len(states))
    diag['pressure_reflection_y']=extension
    return result,diag


def prepare_extended(states,cells,target_count=None):
    from pressure_fit_stencil import mode,tensor_stencil,free_mean_fit
    fit=mode()
    if target_count is None:target_count=len(states)
    degree=int(os.environ.get('LAGFORCE_PRESSURE_DEGREE','2'))
    if degree!=2:
        from polynomial_pressure_reconstruction import prepare as prepare_polynomial
        return prepare_polynomial(states,cells,degree,target_count=target_count)
    statistics=[polygon_statistics(p) for p in cells]
    areas=np.array([v[0] for v in statistics]);centers=np.array([v[1] for v in statistics]);cov=np.array([v[2] for v in statistics]);P=np.array([s['P'] for s in states])
    if np.any(P<=0):raise ValueError('quadratic pressure needs positive cell means')
    result=[];limited=0;mean_error=0.;min_ratio=1.;ridge=1e-14
    for i,(s,poly) in enumerate(zip(states[:target_count],cells[:target_count])):
        scale=np.sqrt(areas[i]);second=cov[i]/areas[i];radius=np.linalg.norm(centers-centers[i],axis=1)/(3*scale)
        ids=np.flatnonzero((radius<1)&(np.arange(len(states))!=i));z=(centers[ids]-centers[i])/scale
        if fit=='tensor_free_mean':
            ids,tensor_weight=tensor_stencil(centers,centers[i],scale,3.)
            z=(centers[ids]-centers[i])/scale
        neighbor_second=cov[ids]/areas[i]+np.einsum('ni,nj->nij',z,z)
        A=np.column_stack([z,neighbor_second[:,0,0]-second[0,0],neighbor_second[:,0,1]-second[0,1],neighbor_second[:,1,1]-second[1,1]])
        r=radius[ids];weight=(1-r)**4*(1+4*r)
        if fit=='tensor_free_mean':coefficient=free_mean_fit(A,P[ids]-P[i],tensor_weight,ridge)
        else:coefficient=np.linalg.solve(A.T@(weight[:,None]*A)+ridge*np.eye(5),A.T@(weight*(P[ids]-P[i])))
        metadata=dict(center=centers[i],scale=scale,second=second,mean=P[i],coefficient=coefficient)
        lower=bernstein_lower_certificate(metadata,poly);theta=1.
        if lower<.01*P[i]:lower=lower_bernstein_bound(metadata,poly)
        if lower<.01*P[i]:
            theta=.99*P[i]/(P[i]-lower);metadata['coefficient']=theta*coefficient;limited+=1
        metadata['positivity_scale']=theta;result.append(metadata);min_ratio=min(min_ratio,theta)
        # Basis geometric means are analytically zero. This diagnostic checks
        # the independent area used by the numerical density representation.
        mean_error=max(mean_error,abs(areas[i]-s['density'].area)/areas[i])
    if mean_error>1e-10:raise RuntimeError('pressure geometric mean normalization mismatch')
    return result,dict(pressure_fit=fit,pressure_reconstruction='quadratic_cell_mean',pressure_quadratic_limited_cells=limited,
        pressure_quadratic_min_positivity_scale=min_ratio,pressure_quadratic_area_relative_error=mean_error)
