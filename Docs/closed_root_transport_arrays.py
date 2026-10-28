"""Feed existing closed massroot directly from its actual transport arrays.

No kinetic formula, gauge, tolerance, sparse solve or summation order changes.
Removes conversion of millions of array rows to tuples and back to arrays.
"""
import numpy as np
from closed_weight_root import solve_arrays

def root_arrays(geometry,transport):
    normal=geometry['normals'];owners=geometry['owners']
    vl,vr=transport['vl'],transport['vr']
    scale=geometry['record_scale']
    defect=abs(np.sum(normal*normal,axis=1)-1)
    if np.any(defect>64*np.finfo(float).eps):raise ValueError('unit face normals required')
    projectionL=(abs(vl)@abs(normal)[:,:,None])[:,:,0]
    projectionR=(abs(vr)@abs(normal)[:,:,None])[:,:,0]
    return dict(i=owners[:,0],j=owners[:,1],scale=scale,base=transport['base'],
                distance=geometry['record_distance'],
                vl=(vl@normal[:,:,None])[:,:,0],vr=(vr@normal[:,:,None])[:,:,0],
                rl=transport['rl'],rr=transport['rr'],epsilon=0.,
                wall_projection_absolute_sum=scale*.25*(transport['rl']*projectionL.sum(1)+transport['rr']*projectionR.sum(1)),
                normal_norm_defect=defect)

def solve_from_transport(geometry,transport,count):
    return solve_arrays(root_arrays(geometry,transport),count)
