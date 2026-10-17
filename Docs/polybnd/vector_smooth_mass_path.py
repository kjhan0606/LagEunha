"""Batched evaluation of the existing continuous-density path mass flux.

Density/velocity reconstruction, positive samples, face quadrature and soft
splitting are identical to smooth_mass_path. Grouping evaluations by cell
avoids constructing many tiny polynomial arrays for individual face points.
"""
import numpy as np
from kinetic_factor_rhs import ALPHA
from kinetic_speed_smoothing import face_speed,parameters,split
from quadratic_velocity_reconstruction import basis


def field_samples(states,factors,mass,owner,points):
    owner=np.asarray(owner);points=np.asarray(points);rho=np.empty(len(owner));mean=np.empty((len(owner),2));samples=np.empty((len(owner),4,2))
    order=np.argsort(owner,kind='stable');owners=owner[order];starts=np.r_[0,1+np.flatnonzero(np.diff(owners)),len(owner)]
    for lo,hi in zip(starts[:-1],starts[1:]):
        ids=order[lo:hi];i=int(owners[lo]);s=states[i];p=points[ids];affine=s['v']+(p-s['c'])@s['B'].T;metadata=s.get('velocity_quadratic')
        if metadata is None:
            velocity=affine;coordinates=np.broadcast_to(ALPHA,(len(ids),4,2))
        else:
            psi=basis(metadata,p);velocity=affine+psi@metadata['coefficient'].T
            coordinates=(psi@metadata['standard_coefficient'].T)[:,None,:]+(ALPHA@metadata['residual_map'].T)[None,:,:]
        rho[ids]=s['density'].value(p);mean[ids]=velocity;samples[ids]=affine[:,None,:]+coordinates@np.asarray(factors[i]).T/np.sqrt(mass[i])
    return rho,mean,samples


def evaluate(states,edges,sites,factors,mass,gen,z,order,boundary_velocity):
    N=len(states);sites=np.asarray(sites,float);gen=np.asarray(gen);z=np.asarray(z);mass=np.asarray(mass)
    i=np.array([e[0] for e in edges],int);j=np.array([e[1] for e in edges],int);internal=j>=0
    start=np.array([e[2] for e in edges]);end=np.array([e[3] for e in edges]);normal=np.array([e[4] for e in edges]);length=np.array([e[5] for e in edges])
    g,w=np.polynomial.legendre.leggauss(order);g=(g+1)/2;w=w/2;points=(start[:,None,:]+g[None,:,None]*(end-start)[:,None,:]).reshape(-1,2)
    left=np.repeat(i,order);right=np.repeat(j,order);inside=right>=0;count=len(points);ids=np.r_[np.arange(count),np.flatnonzero(inside)]
    rho,velocity,samples=field_samples(states,factors,mass,np.r_[left,right[inside]],points[ids])
    rl=rho[:count];rr=rl.copy();sl=samples[:count];sr=sl.copy();rr[inside]=rho[count:];sr[inside]=samples[count:]
    normals=np.repeat(normal,order,axis=0);wn=np.empty(count,dtype=np.result_type(gen,z,float));delta=sites[right[inside]]-sites[left[inside]];distance=np.linalg.norm(delta,axis=1)
    wn[inside]=(np.sum(delta*gen[right[inside]],axis=1)-np.sum((points[inside]-sites[left[inside]])*(gen[right[inside]]-gen[left[inside]]),axis=1)+.5*(z[left[inside]]-z[right[inside]]))/distance
    outside=~inside
    if np.any(outside):
        wn[outside]=np.array([boundary_velocity(p)@n for p,n in zip(points[outside],normals[outside])]);v=sl[outside];n=normals[outside]
        sr[outside]=v-2*(np.einsum('nki,ni->nk',v,n)-wn[outside,None])[:,:,None]*n[:,None,:]
    left_speed=np.einsum('nki,ni->nk',sl,normals)-wn[:,None];right_speed=np.einsum('nki,ni->nk',sr,normals)-wn[:,None]
    coefficient,unit_length=parameters();epsilon=np.repeat([face_speed(states,a,b,coefficient,unit_length) for a,b in zip(i,j)],order)
    if coefficient:
        plus,_,_,_=split(left_speed,epsilon[:,None]);_,minus,_,_=split(right_speed,epsilon[:,None])
    else:plus=np.maximum(left_speed,0.);minus=np.minimum(right_speed,0.)
    flux=.25*(rl*np.sum(plus,axis=1)+rr*np.sum(minus,axis=1));edge_flux=length*np.sum(w[None,:]*flux.reshape(-1,order),axis=1)
    # Retain the scalar edge accumulation order: left then right per edge.
    owners=np.column_stack([i,j]).ravel();values=np.column_stack([edge_flux,-edge_flux]).ravel();valid=owners>=0;residual=np.zeros(N);np.add.at(residual,owners[valid],values[valid])
    return -residual
