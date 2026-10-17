"""Resolve velocity curvature within the stored positive covariance budget.

Quadratic basis functions are orthogonal to1,x under the numerical positive
density. A continuous dominant covariance direction preserves m,p,D,Pi while
moving part of R from independent velocity samples to a resolved mean field.
No stored moment or covariance eigenvalue is projected or clipped.
"""
import numpy as np
import os


def fit_mode():
    name=os.environ.get('LAGFORCE_VELOCITY_FIT','radial_fixed_mean')
    if name not in ('radial_fixed_mean','tensor_free_mean','tensor_free_affine'):
        raise ValueError('unknown velocity fit')
    return name


def curvature_fit(design,target,weight,mode,zmean=None):
    if mode=='radial_fixed_mean':
        return np.linalg.solve(design.T@(weight[:,None]*design)+1e-10*np.eye(3),
            design.T@(weight[:,None]*target))
    if mode not in ('tensor_free_mean','tensor_free_affine'):raise ValueError('unknown velocity fit')
    nuisance=np.ones((len(design),1))
    if mode=='tensor_free_affine':
        if zmean is None or np.shape(zmean)!=(len(design),2):raise ValueError('free affine fit requires neighbor coordinates')
        nuisance=np.column_stack([nuisance,zmean])
    matrix=np.sqrt(weight)[:,None]*np.column_stack([nuisance,design]);columns=matrix.shape[1]
    augmented=np.vstack([matrix,np.sqrt(1e-10)*np.eye(columns)])
    rhs=np.vstack([np.sqrt(weight)[:,None]*target,np.zeros((columns,2))])
    # Discard fitted constant/linear nuisance terms. The reconstructed
    # affine part remains prescribed by stored p,D; only curvature uses
    # this fit through the existing density-orthogonal basis.
    return np.linalg.lstsq(augmented,rhs,rcond=None)[0][nuisance.shape[1]:]


def closure_name():
    mode=os.environ.get('LAGFORCE_QUADRATIC_BUDGET','dominant')
    if mode=='dominant':return 'quadratic_budget_continuous'
    if mode=='projected':return 'quadratic_budget_projected'
    if mode=='projected_regularized':return 'quadratic_budget_projected_regularized'
    raise ValueError('unknown quadratic covariance budget')


def budget_backend_name():return 'batched' if os.environ.get('LAGFORCE_VECTOR_QUADRATIC_BUDGET','0')=='1' else 'scalar'


def projected_budget(L,mass,candidate,gram_factor,iterations=64,regularized=False):
    # Fixed-count accelerated projected gradient in covariance coordinates:
    # min ||L C - sqrt(m) H0 Ghalf.T||_F^2, ||C||_2 <= cap.
    # Unlike whitening H0 by L^{-1}, the physical objective has a continuous
    # finite-iteration selection at rank loss. All iterates are feasible.
    scale=np.linalg.svd(L,compute_uv=False)[0]
    if scale==0:return np.zeros_like(candidate),np.zeros_like(candidate),np.eye(2)
    physical=np.sqrt(mass)*candidate@gram_factor.T
    if regularized:scale=np.hypot(scale,np.linalg.norm(physical))
    F=L/scale;target=physical/scale
    current=np.zeros_like(candidate);extrapolated=current.copy();acceleration=1.;cap=1-1e-12
    for _ in range(iterations):
        trial=extrapolated-F.T@(F@extrapolated-target)
        left,values,right=np.linalg.svd(trial,full_matrices=False)
        allowed=np.minimum(values,cap);updated=(left*allowed)@right
        next_acceleration=(1+np.sqrt(1+4*acceleration**2))/2
        extrapolated=updated+(acceleration-1)/next_acceleration*(updated-current)
        current=updated;acceleration=next_acceleration
    standard=np.linalg.solve(gram_factor,current.T).T
    T0=(left*np.sqrt(1-allowed**2))@left.T
    # Polar rotation puts the residual samples on the unique symmetric
    # physical square root without subtracting nearly equal covariances.
    u,_,vh=np.linalg.svd(L@T0)
    residual=T0@vh.T@u.T
    return L@standard/np.sqrt(mass),standard,residual


def projected_budget_batch(factors,mass,candidates,gram_factors,iterations=64,regularized=False):
    """Same fixed iteration algorithm, with the independent cell axis batched."""
    L=np.asarray(factors);mass=np.asarray(mass);candidate=np.asarray(candidates);G=np.asarray(gram_factors)
    if len(mass)==0:return np.empty((0,2,3)),np.empty((0,2,3)),np.empty((0,2,2))
    physical=np.sqrt(mass)[:,None,None]*(candidate@np.swapaxes(G,-1,-2))
    scale=np.linalg.svd(L,compute_uv=False)[:,0]
    if regularized:scale=np.hypot(scale,np.linalg.norm(physical,axis=(-2,-1)))
    scale=np.where(scale>0,scale,1.);F=L/scale[:,None,None];target=physical/scale[:,None,None]
    current=np.zeros_like(candidate);extrapolated=current.copy();acceleration=1.;cap=1-1e-12
    for _ in range(iterations):
        trial=extrapolated-np.swapaxes(F,-1,-2)@(F@extrapolated-target)
        left,values,right=np.linalg.svd(trial,full_matrices=False);allowed=np.minimum(values,cap)
        updated=(left*allowed[:,None,:])@right
        next_acceleration=(1+np.sqrt(1+4*acceleration**2))/2
        extrapolated=updated+(acceleration-1)/next_acceleration*(updated-current)
        current=updated;acceleration=next_acceleration
    standard=np.swapaxes(np.linalg.solve(G,np.swapaxes(current,-1,-2)),-1,-2)
    T0=(left*np.sqrt(1-allowed**2)[:,None,:])@np.swapaxes(left,-1,-2)
    u,_,vh=np.linalg.svd(L@T0)
    residual=T0@np.swapaxes(vh,-1,-2)@np.swapaxes(u,-1,-2)
    return L@standard/np.sqrt(mass)[:,None,None],standard,residual


def basis(metadata,points):
    z=(np.atleast_2d(points)-metadata['center'])/metadata['scale']
    phi=np.column_stack([z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])
    return phi-np.column_stack([np.ones(len(z)),z])@metadata['projection']


def correction(state,points):
    metadata=state.get('velocity_quadratic')
    if metadata is None:return np.zeros_like(np.atleast_2d(points))
    return basis(metadata,points)@metadata['coefficient'].T


def velocity_samples(state,L,mass,points,alpha):
    points=np.atleast_2d(points)
    affine=state['v']+(points-state['c'])@state['B'].T
    metadata=state.get('velocity_quadratic')
    if metadata is None:
        coordinates=np.broadcast_to(alpha,(len(points),len(alpha),2))
    else:
        resolved=basis(metadata,points)@metadata['standard_coefficient'].T
        coordinates=resolved[:,None,:]+(alpha@metadata['residual_map'].T)[None,:,:]
    velocities=affine[:,None,:]+coordinates@L.T/np.sqrt(mass)
    return velocities,coordinates


def covariance_budget(L,mass,candidate,gram_factor):
    name=closure_name()
    if name.startswith('quadratic_budget_projected'):return projected_budget(L,mass,candidate,gram_factor,regularized=name.endswith('_regularized'))
    # Resolve curvature in the dominant physical covariance direction. The
    # anisotropy weight vanishes at eigenvalue degeneracy, removing the
    # ambiguous eigenvector. Tiny transverse variance cannot suppress the
    # finite longitudinal mean, as it did with whitened SVD clipping.
    U,singular,Vh=np.linalg.svd(L)
    largest,smallest=singular
    if largest==0:return np.zeros_like(candidate),np.zeros_like(candidate),np.eye(2)
    weight=1-(smallest/largest)**2
    physical=np.sqrt(mass)*candidate@gram_factor.T
    row=U[:,0]@physical
    size=np.linalg.norm(row);cap=1-1e-12
    normalized=row/max(largest,size/cap)
    eigencoordinates=np.vstack([weight*normalized,np.zeros_like(normalized)])
    bounded=Vh.T@eigencoordinates
    standard=np.linalg.solve(gram_factor,bounded.T).T
    remaining=np.sqrt(max(0.,1-weight**2*np.dot(normalized,normalized)))
    # L @ residual is the symmetric positive physical square root. Thus the
    # discrete velocity directions do not depend on the stored factor gauge.
    residual=Vh.T@np.diag([remaining,1.])@U.T
    coefficient=L@standard/np.sqrt(mass)
    return coefficient,standard,residual


def prepare(states,cells,edges,mass,factors):
    N=len(states)
    fit_name=fit_mode()
    centers=np.asarray([s['c'] for s in states]);means=np.asarray([s['v'] for s in states])
    second=np.asarray([s['raw'][3:7].reshape(2,2)/m for s,m in zip(states,mass)])
    active=rank_deficient=0;split_error=orthogonal_error=0.;fractions=[];prepared=[]
    for i,(s,L,m) in enumerate(zip(states,factors,mass)):
        R=L@L.T;norm=float(np.linalg.norm(R))
        if norm==0:continue
        center=s['c'];scale=np.sqrt(s['density'].area)
        radius=np.linalg.norm(centers-center,axis=1)/(3*scale)
        ids=np.flatnonzero((radius<1)&(np.arange(N)!=i))
        if fit_name.startswith('tensor_free_'):
            from pressure_fit_stencil import tensor_stencil
            ids,tensor_weight=tensor_stencil(centers,center,scale,3.)
        if len(ids)<3:rank_deficient+=1
        points=s['density'].points;prob=s['density'].weighted_mass/m
        z=(points-center)/scale;X=np.column_stack([np.ones(len(z)),z])
        phi=np.column_stack([z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])
        projection=np.linalg.solve(X.T@(prob[:,None]*X),X.T@(prob[:,None]*phi))
        psi=phi-X@projection;_,gram_factor=np.linalg.qr(np.sqrt(prob)[:,None]*psi,mode='reduced')
        cj=centers[ids];zmean=(cj-center)/scale
        # Neighbor expectations of a quadratic use only existing m,q,Q.
        c=np.asarray(center,dtype=np.longdouble);cc=np.asarray(cj,dtype=np.longdouble)
        zz=(np.asarray(second[ids],dtype=np.longdouble)-c[None,:,None]*cc[:,None,:]-cc[:,:,None]*c[None,None,:]+np.outer(c,c))/scale**2
        qmean=np.column_stack([zz[:,0,0],zz[:,0,1],zz[:,1,1]]).astype(float)
        design=qmean-np.column_stack([np.ones(len(ids)),zmean])@projection
        target=means[ids]-s['v']-(cj-center)@s['B'].T
        # Wendland C2 weights vanish continuously at the patch boundary.
        # A fixed positive ridge also makes deficient boundary fits unique.
        r=radius[ids];weight=(1-r)**4*(1+4*r)
        if fit_name.startswith('tensor_free_'):weight=tensor_weight
        fit=curvature_fit(design,target,weight,fit_name,zmean=zmean)
        orthogonal_error=max(orthogonal_error,float(np.max(abs(X.T@(prob[:,None]*psi)))))
        prepared.append((i,L,m,R,norm,center,scale,projection,fit.T,gram_factor))
    batched=budget_backend_name()=='batched'
    if batched and not closure_name().startswith('quadratic_budget_projected'):raise ValueError('batched quadratic budget requires projected closure')
    if batched:
        output=projected_budget_batch([v[1] for v in prepared],[v[2] for v in prepared],[v[8] for v in prepared],[v[9] for v in prepared],regularized=closure_name().endswith('_regularized'))
        budgets=zip(*output)
    else:budgets=(covariance_budget(v[1],v[2],v[8],v[9]) for v in prepared)
    for item,(coefficient,standard,residual) in zip(prepared,budgets):
        i,L,m,R,norm,center,scale,projection,_,gram_factor=item;s=states[i]
        metadata=dict(center=center,scale=scale,projection=projection,coefficient=coefficient,
            standard_coefficient=standard,residual_map=residual)
        s['velocity_quadratic']=metadata;active+=1
        resolved=m*coefficient@(gram_factor.T@gram_factor)@coefficient.T
        split_error=max(split_error,float(np.linalg.norm(resolved+L@residual@residual.T@L.T-R))/norm)
        fractions.append(float(np.trace(resolved)/np.trace(R)))
    if split_error>1e-9 or orthogonal_error>1e-10:
        raise RuntimeError('quadratic velocity moment budget failed: '+str(dict(split_error=split_error,orthogonal_error=orthogonal_error)))
    return dict(velocity_fit=fit_name,velocity_quadratic_cells=active,velocity_quadratic_rank_deficient=rank_deficient,
        velocity_quadratic_covariance_split_error=split_error,velocity_quadratic_orthogonality_error=orthogonal_error,
        velocity_quadratic_max_resolved_fraction=max(fractions,default=0.),velocity_quadratic_budget_backend='batched' if batched else 'scalar')
