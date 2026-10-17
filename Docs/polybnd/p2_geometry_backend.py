"""Optional signed triangle geometry kernels; no remap/force/physics changes."""
import ctypes,os
from pathlib import Path
import numpy as np
_lib=None

def all_geom(keys,X):
    ids=np.ascontiguousarray(keys,dtype=np.int64);X=np.ascontiguousarray(X,dtype=np.float64)
    n=len(ids);mode=os.environ.get('P2_GEOMETRY','numpy')
    if mode=='numpy':
        p=X[ids];a,b,c=p[:,0],p[:,1],p[:,2]
        area=.5*((b[:,0]-a[:,0])*(c[:,1]-a[:,1])-(c[:,0]-a[:,0])*(b[:,1]-a[:,1]))
        e=np.stack([c-b,a-c,b-a],axis=1);da=.5*np.stack([-e[:,:,1],e[:,:,0]],axis=2)
        return ids,area,da,p
    if mode!='openmp':raise ValueError('unknown geometry backend '+mode)
    global _lib
    if _lib is None:
        _lib=ctypes.CDLL(str(Path(__file__).with_name('p2_geometry_omp.so')))
        vec=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        ind=np.ctypeslib.ndpointer(dtype=np.int64,flags='C_CONTIGUOUS')
        _lib.p2_geometry.argtypes=[ctypes.c_int64,ind,vec,vec,vec,vec];_lib.p2_geometry.restype=None
        _lib.p2_threads.restype=ctypes.c_int
    area=np.empty(n);da=np.empty((n,3,2));p=np.empty((n,3,2))
    _lib.p2_geometry(n,ids,X,area,da,p)
    return ids,area,da,p

def threads_used():return _lib.p2_threads() if _lib is not None else None
