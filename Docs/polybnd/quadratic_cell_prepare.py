"""Independent cell of the existing radial/area quadratic preparation.

Algebra and fixed-count projected covariance budget follow quad.prepare.
No analytical velocity values or new smoothing/regularization are introduced.
"""
import numpy as np
import quadratic_velocity_reconstruction as quad

def prepare_cell(i,state,density,mass,L,centers,means,second,tree):
    if quad.fit_mode()!='radial_fixed_mean' or quad.stencil_mode()!='area':raise ValueError('approved radial/area reconstruction required')
    R=L@L.T;norm=float(np.linalg.norm(R))
    if norm==0:raise ValueError('strictly positive covariance factor required')
    center=state['c'];scale=np.sqrt(density.area);ids,r=quad.radial_neighbors(centers,center,scale,i,tree)
    points=density.points;prob=density.weighted_mass/mass
    z=(points-center)/scale;X=np.column_stack([np.ones(len(z)),z]);phi=np.column_stack([z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])
    projection=np.linalg.solve(X.T@(prob[:,None]*X),X.T@(prob[:,None]*phi));psi=phi-X@projection
    _,gram_factor=np.linalg.qr(np.sqrt(prob)[:,None]*psi,mode='reduced')
    cj=centers[ids];zmean=(cj-center)/scale;c=np.asarray(center,dtype=np.longdouble);cc=np.asarray(cj,dtype=np.longdouble)
    zz=(np.asarray(second[ids],dtype=np.longdouble)-c[None,:,None]*cc[:,None,:]-cc[:,:,None]*c[None,None,:]+np.outer(c,c))/scale**2
    qmean=np.column_stack([zz[:,0,0],zz[:,0,1],zz[:,1,1]]).astype(float)
    design=qmean-np.column_stack([np.ones(len(ids)),zmean])@projection
    target=means[ids]-state['v']-(cj-center)@state['B'].T
    weight=(1-r)**4*(1+4*r);fit=quad.curvature_fit(design,target,weight,quad.fit_mode(),zmean=zmean)
    coefficient,standard,residual=quad.covariance_budget(L,mass,fit.T,gram_factor)
    metadata=dict(center=center,scale=scale,projection=projection,coefficient=coefficient,standard_coefficient=standard,residual_map=residual)
    resolved=mass*coefficient@(gram_factor.T@gram_factor)@coefficient.T
    split_error=float(np.linalg.norm(resolved+L@residual@residual.T@L.T-R))/norm
    orthogonality=float(np.max(abs(X.T@(prob[:,None]*psi))))
    if split_error>1e-9 or orthogonality>1e-10:raise RuntimeError(dict(status='quadratic_budget_failed',split_error=split_error,orthogonality=orthogonality))
    return metadata,dict(covariance_split_error=split_error,orthogonality=orthogonality,resolved_fraction=float(np.trace(resolved)/np.trace(R)),neighbors=len(ids))
