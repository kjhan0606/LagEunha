"""Shared positive face quadrature split at every affine sample-speed zero."""
import numpy as np
from split_face_mass_3d import clip_face


def piece_rule(vertices,normal,order):
    """Drop only exactly zero-area fan triangles created by coincident cuts."""
    p=clip_face(vertices,normal,lambda x:np.ones(len(x)))
    if len(p)<3:return np.empty((0,3)),np.empty(0)
    center=p.mean(0);z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2
    a,b=np.meshgrid(z,z,indexing='ij');wa,wb=np.meshgrid(w,w,indexing='ij')
    a=a.ravel();b=b.ravel();weights=(wa*wb).ravel()*a
    bary=np.column_stack([1-a,a*(1-b),a*b]);points=[];out=[]
    for left,right in zip(p,np.roll(p,-1,axis=0)):
        jac=np.linalg.norm(np.cross(left-center,right-center))
        if jac==0:continue
        points.append(bary@np.array([center,left,right]));out.append(weights*jac)
    if not points:return np.empty((0,3)),np.empty(0)
    return np.vstack(points),np.concatenate(out)


def rule(vertices,normal,speeds,order=4):
    pieces=[np.asarray(vertices)]
    count=np.asarray(speeds(pieces[0])).shape[1]
    for column in range(count):
        next_pieces=[]
        for p in pieces:
            values=speeds(p)[:,column]
            if values.min()>=0 or values.max()<=0:
                next_pieces.append(p);continue
            for sign in (1,-1):
                q=clip_face(p,normal,lambda x:sign*speeds(x)[:,column])
                if len(q)>=3:next_pieces.append(q)
        pieces=next_pieces
    points=[];weights=[]
    for p in pieces:
        x,w=piece_rule(p,normal,order)
        if len(w):points.append(x);weights.append(w)
    if not points:raise ValueError('split face has zero total area')
    return np.vstack(points),np.concatenate(weights)
