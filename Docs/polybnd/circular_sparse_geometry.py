"""Sparse native exact power supports; fixed box is an enclosure for the disk.

The exported polygons must still be intersected with the true circular wall.
They are NOT a chord approximation to that wall or a usable fluid mesh alone.
"""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

class SparsePowerGeometry:
    def __init__(self,library):
        self.library=Path(library).resolve(strict=True)
        self.lib=ctypes.CDLL(str(self.library))
        int64=np.ctypeslib.ndpointer(dtype=np.int64,ndim=1,flags='C_CONTIGUOUS')
        float64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        int32=np.ctypeslib.ndpointer(dtype=np.int32,ndim=1,flags='C_CONTIGUOUS')
        self.lib.circular_power_sparse_resolve.argtypes=[ctypes.c_int,ctypes.c_int,int64,float64,float64,int64,int64,ctypes.c_double,int32,ctypes.c_char_p]
        self.lib.circular_power_sparse_resolve.restype=ctypes.c_int

    def resolve(self,sites,weights,ids,radius,output,seed_count=12,threads=4):
        sites=np.ascontiguousarray(sites,dtype=float);weights=np.ascontiguousarray(weights,dtype=float);ids=np.ascontiguousarray(ids,dtype=np.int64)
        n=len(ids);output=Path(output)
        if sites.shape!=(n,2) or weights.shape!=(n,) or len(set(map(int,ids)))!=n or np.any(ids<0):raise ValueError('finite distinct stable-ID registry required')
        if not np.isfinite(sites).all() or not np.isfinite(weights).all() or radius<=0 or not np.isfinite(radius):raise ValueError('invalid finite registry/domain')
        if output.exists():raise FileExistsError(output)
        rows=[[] for _ in ids]
        if seed_count:
            _,seeds=cKDTree(sites).query(sites,k=min(n,seed_count+1))
            seeds=np.asarray(seeds).reshape(n,-1)
            rows=[sorted((int(j) for j in row if j!=i),key=lambda j:ids[j]) for i,row in enumerate(seeds)]
        offsets=np.ascontiguousarray(np.r_[0,np.cumsum([len(row) for row in rows])],dtype=np.int64)
        candidates=np.ascontiguousarray([j for row in rows for j in row],dtype=np.int64)
        status=np.full(n,-999,dtype=np.int32)
        np.savez(output.with_suffix('.inputs.npz'),sites=sites,weights=weights,ids=ids,radius=radius,offsets=offsets,candidates=candidates)
        code=self.lib.circular_power_sparse_resolve(threads,n,ids,sites,weights,offsets,candidates,radius,status,str(output).encode())
        np.save(output.with_suffix('.status.npy'),status)
        if code:
            failure=dict(code=code,status=status.tolist(),fallback=False)
            output.with_suffix('.failure.json').write_text(json.dumps(failure,indent=2));raise RuntimeError(failure)
        result=json.loads(output.read_text())
        result['library_sha256']=hashlib.sha256(self.library.read_bytes()).hexdigest()
        assert result['registry_count']==n and result['shared_exact_reversed'] and np.all(status==0)
        return result
