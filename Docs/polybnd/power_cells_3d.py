"""Actual C Laguerre cells and shared face quadrature in an axis-aligned box.

Load an explicitly supplied library built from test_volume_wrapper.c and the
validated CPU laguerre.c. Current C topology limits are propagated as errors.
"""
import ctypes
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull


def face_quadrature(vertices,normal,order=4):
    face=np.asarray(vertices,dtype=float);n=np.asarray(normal,dtype=float)
    center=face.mean(axis=0);axis=np.eye(3)[np.argmin(abs(n))]
    u=np.cross(n,axis);u/=np.linalg.norm(u);v=np.cross(n,u)
    angle=np.arctan2((face-center)@v,(face-center)@u);face=face[np.argsort(angle)]
    z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2
    a,b=np.meshgrid(z,z,indexing='ij');wa,wb=np.meshgrid(w,w,indexing='ij')
    a=a.ravel();b=b.ravel();weights=(wa*wb).ravel()*a
    bary=np.column_stack([1-a,a*(1-b),a*b]);points=[];out=[]
    for left,right in zip(face,np.roll(face,-1,axis=0)):
        jac=np.linalg.norm(np.cross(left-center,right-center))
        if jac<=0:raise ValueError('degenerate face triangle')
        points.append(bary@np.array([center,left,right]));out.append(weights*jac)
    return np.vstack(points),np.concatenate(out)


class Geometry3D:
    def __init__(self,library):
        self.lib=ctypes.CDLL(str(Path(library).resolve()))
        array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.lib.cell.argtypes=[ctypes.c_int,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_double)]
        self.lib.cell.restype=ctypes.c_int

    def build(self,sites,weights,bounds=None,face_order=4):
        x=np.asarray(sites,dtype=float);w=np.asarray(weights,dtype=float)
        bounds=np.array([[0.,0.,0.],[1.,1.,1.]]) if bounds is None else np.asarray(bounds,dtype=float)
        if x.ndim!=2 or x.shape[1]!=3 or w.shape!=(len(x),) or bounds.shape!=(2,3):raise ValueError('invalid 3D geometry arrays')
        if not all(np.all(np.isfinite(a)) for a in (x,w,bounds)):raise ValueError('finite geometry required')
        extent=bounds[1]-bounds[0]
        if np.any(extent<=0) or np.any(extent>2) or np.any(x<=bounds[0]) or np.any(x>=bounds[1]):raise ValueError('interior sites and box extent in (0,2] required by wrapper')
        normals=np.vstack([np.eye(3),-np.eye(3)]);polys=[];volumes=[];face_candidates={}
        for i,site in enumerate(x):
            ids=np.flatnonzero(np.arange(len(x))!=i);delta=x[ids]-site
            distances=np.r_[bounds[1]-site,site-bounds[0]];ghosts=2*distances[:,None]*normals
            data=np.ascontiguousarray(np.vstack([np.column_stack([delta,w[ids]]),np.column_stack([ghosts,np.full(6,w[i])])]))
            out=np.zeros(1536);volume=ctypes.c_double();count=self.lib.cell(len(data),w[i],data,out,ctypes.byref(volume))
            if count<4:raise RuntimeError('C 3D cell empty or unsupported: '+str((i,count)))
            vertices=out[:3*count].reshape(-1,3)+site;polys.append(vertices);volumes.append(volume.value)
            for k,(dr,wj) in enumerate(zip(data[:,:3],data[:,3])):
                d=np.linalg.norm(dr);bound=.5*(dr@dr+w[i]-wj)
                face=vertices[abs((vertices-site)@dr-bound)<2e-11*d]
                if len(face)<3:continue
                n=dr/d;points,fw=face_quadrature(face,n,face_order)
                j=int(ids[k]) if k<len(ids) else -1
                key=(i,j) if j>=0 else (i,-1-k+len(ids))
                face_candidates[key]=dict(i=i,j=j,normal=n,points=points,weights=fw,area=float(fw.sum()),vertices=face)
        mismatch=0.;faces=[]
        for (i,j),f in face_candidates.items():
            if j>=0:
                other=face_candidates.get((j,i))
                if other is None:raise RuntimeError('missing reciprocal 3D face')
                mismatch=max(mismatch,abs(f['area']-other['area']))
                if i>j:continue
            faces.append(f)
        if mismatch>1e-10:raise RuntimeError('3D shared face area mismatch')
        volume_error=abs(sum(volumes)-np.prod(extent))
        if volume_error>1e-10*np.prod(extent):raise RuntimeError('3D box not partitioned')
        return polys,faces,dict(volumes=np.array(volumes),partition_error=volume_error,shared_face_area_mismatch=mismatch)
