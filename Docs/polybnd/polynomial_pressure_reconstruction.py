"""Positive degree3/4 pressure with the existing thermal-energy cell mean.

Only numerical cell energies and exact polygon geometric moments are used.
Centered moments avoid global-coordinate fourth-moment cancellation. This is
an explicit alternative spatial closure, not a change to time integration.
"""
from functools import lru_cache
from math import comb,factorial
import numpy as np
import os


def evaluation_backend_name(degree=None):
    if degree is None:degree=int(os.environ.get('LAGFORCE_PRESSURE_DEGREE','2')) if os.environ.get('LAGFORCE_QUADRATIC_PRESSURE','0')=='1' else 0
    return 'horner' if degree>2 and os.environ.get('LAGFORCE_FAST_PRESSURE_POLYNOMIAL','0')=='1' else 'monomial'


def horner_coefficients(metadata):
    coefficient=np.asarray(metadata['coefficient']);return tuple(float(v) for v in np.r_[metadata['mean']-metadata['basis_means']@coefficient,coefficient,np.zeros(14-len(coefficient))])


def horner(x,y,c):
    return c[0]+y*(c[2]+y*(c[5]+y*(c[9]+y*c[14])))+x*(c[1]+y*(c[4]+y*(c[8]+y*c[13]))+x*(c[3]+y*(c[7]+y*c[12])+x*(c[6]+y*c[11]+x*c[10])))


def powers(degree):return [(total-k,k) for total in range(degree+1) for k in range(total+1)]


def monomials(z,exponents):
    z=np.atleast_2d(z)
    return np.column_stack([z[:,0]**a*z[:,1]**b for a,b in exponents])


def polygon_rule(poly):
    # Three positive Duffy nodes in each direction exactly integrate all
    # degree4 polynomials (the radial Jacobian raises degree to at most5).
    poly=np.asarray(poly);center=poly.mean(0);g,w=np.polynomial.legendre.leggauss(3);g=(g+1)/2;w=w/2;points=[];weights=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        u=a-center;v=b-center;det=abs(u[0]*v[1]-u[1]*v[0])
        points.append(center+g[:,None,None]*u+(1-g[:,None,None])*g[None,:,None]*v)
        weights.append(det*(1-g[:,None])*w[:,None]*w[None,:])
    return np.concatenate([q.reshape(-1,2) for q in points]),np.concatenate([q.ravel() for q in weights])


@lru_cache(maxsize=2)
def bernstein_transform(degree):
    alpha=np.array([(i,j,degree-i-j) for i in range(degree+1) for j in range(degree-i+1)])
    nodes=alpha/degree;factor=np.array([factorial(degree)/np.prod([factorial(int(v)) for v in a]) for a in alpha])
    V=factor[None,:]*np.prod(nodes[:,None,:]**alpha[None,:,:],axis=2)
    inverse=np.linalg.inv(V)
    return nodes,inverse,float(np.max(np.sum(abs(inverse),axis=1)))


def value(metadata,points):
    if 'horner_coefficients' in metadata:
        p=np.asarray(points)
        if p.ndim==1:
            x=float((p[0]-metadata['center'][0])/metadata['scale']);y=float((p[1]-metadata['center'][1])/metadata['scale'])
            return np.array([horner(x,y,metadata['horner_coefficients'])])
        z=(p-metadata['center'])/metadata['scale'];return horner(z[:,0],z[:,1],metadata['horner_coefficients'])
    z=(np.atleast_2d(points)-metadata['center'])/metadata['scale']
    return metadata['mean']+(monomials(z,metadata['powers'])-metadata['basis_means'])@metadata['coefficient']


def lower_bound(metadata,poly):
    nodes,inverse,_=bernstein_transform(metadata['degree']);poly=np.asarray(poly);center=poly.mean(0);lower=float('inf')
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        triangle=np.array([center,a,b]);v=value(metadata,nodes@triangle);coefficients=inverse@v
        lower=min(lower,float(coefficients.min()))
    return lower


def cheap_bound(metadata,poly):
    _,_,operator=bernstein_transform(metadata['degree']);z=(np.asarray(poly)-metadata['center'])/metadata['scale'];x,y=np.max(abs(z),axis=0)
    envelope=np.array([x**a*y**b for a,b in metadata['powers']])+abs(metadata['basis_means'])
    variation=operator*float(abs(metadata['coefficient'])@envelope)
    return metadata['mean']-variation-128*np.finfo(float).eps*(abs(metadata['mean'])+variation)


def shifted_expectations(centers,scales,moments,exponents,ids,center,scale):
    delta=(centers[ids]-center)/scale;ratio=scales[ids]/scale;indices={p:i for i,p in enumerate(exponents)};result=[]
    for a,b in exponents[1:]:
        q=np.zeros(len(ids))
        for k in range(a+1):
            for l in range(b+1):q+=comb(a,k)*comb(b,l)*delta[:,0]**(a-k)*delta[:,1]**(b-l)*ratio**(k+l)*moments[ids,indices[(k,l)]]
        result.append(q)
    return np.column_stack(result)


def prepare(states,cells,degree,target_count=None):
    from pressure_fit_stencil import mode,tensor_stencil,free_mean_fit
    fit=mode()
    if target_count is None:target_count=len(states)
    if degree not in (3,4):raise ValueError('polynomial pressure degree must be3 or4')
    from quadratic_pressure_reconstruction import polygon_statistics
    exponents=powers(degree);statistics=[polygon_statistics(p) for p in cells];areas=np.array([v[0] for v in statistics]);centers=np.array([v[1] for v in statistics]);scales=np.sqrt(areas);moments=[]
    for poly,c,h in zip(cells,centers,scales):
        points,w=polygon_rule(poly);moments.append((w@monomials((points-c)/h,exponents))/w.sum())
    moments=np.asarray(moments)
    # Exactly retain the stored thermal energy using the geometric area.
    means=np.array([s['P']*s['density'].area/A for s,A in zip(states,areas)])
    if np.any(means<=0):raise ValueError('polynomial pressure requires positive mean energy')
    ncoeff=len(exponents)-1;result=[];limited=deficient=0;minimum_theta=1.;ridge=1e-14;support=degree+1.
    for i,(s,poly) in enumerate(zip(states[:target_count],cells[:target_count])):
        radius=np.linalg.norm(centers-centers[i],axis=1)/(support*scales[i]);ids=np.flatnonzero((radius<1)&(np.arange(len(states))!=i))
        if fit=='tensor_free_mean':ids,tensor_weight=tensor_stencil(centers,centers[i],scales[i],support)
        A=shifted_expectations(centers,scales,moments,exponents,ids,centers[i],scales[i])-moments[i,1:]
        r=radius[ids];weight=(1-r)**4*(1+4*r);matrix=np.sqrt(weight)[:,None]*A;target=np.sqrt(weight)*(means[ids]-means[i])
        # Positive ridge makes the selection unique even on deficient stencils.
        augmented=np.vstack([matrix,np.sqrt(ridge)*np.eye(ncoeff)]);rhs=np.r_[target,np.zeros(ncoeff)]
        coefficient,_,rank,_=np.linalg.lstsq(augmented,rhs,rcond=None)
        if fit=='tensor_free_mean':coefficient=free_mean_fit(A,means[ids]-means[i],tensor_weight,ridge)
        if len(ids)<ncoeff:deficient+=1
        metadata=dict(degree=degree,center=centers[i],scale=scales[i],powers=exponents[1:],basis_means=moments[i,1:],mean=means[i],coefficient=coefficient)
        bound=cheap_bound(metadata,poly);theta=1.
        if bound<.01*means[i]:
            bound=lower_bound(metadata,poly)
            if bound<.01*means[i]:theta=.99*means[i]/(means[i]-bound);metadata['coefficient']=theta*coefficient;limited+=1
        metadata['positivity_scale']=theta
        if evaluation_backend_name(degree)=='horner':metadata['horner_coefficients']=horner_coefficients(metadata)
        result.append(metadata);minimum_theta=min(minimum_theta,theta)
    return result,dict(pressure_fit=fit,pressure_reconstruction='polynomial_cell_mean_degree'+str(degree),pressure_polynomial_degree=degree,pressure_polynomial_support_radius=support,
        pressure_polynomial_limited_cells=limited,pressure_polynomial_min_positivity_scale=minimum_theta,pressure_polynomial_neighbor_count_deficient_cells=deficient,pressure_polynomial_evaluation_backend=evaluation_backend_name(degree))
