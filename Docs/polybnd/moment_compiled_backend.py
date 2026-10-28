"""Transactional C ABI wrapper for the validated planar moment operators."""
import ctypes as C
from pathlib import Path
import numpy as np

def f64(a):return np.ascontiguousarray(a,dtype=np.float64)
def i64(a):return np.ascontiguousarray(a,dtype=np.int64)
def ptr(a):return a.ctypes.data_as(C.c_void_p)

class Backend:
    def __init__(self,library,kind='cuda',threads=1):
        self.path=Path(library).resolve(strict=True);self.lib=C.CDLL(str(self.path))
        self.kind={'cpu':0,'cuda':1}[kind];self.threads=threads
        self.lib.lagmoment_stage.argtypes=[C.c_int]*4+[C.c_int64]*2+[C.c_void_p]*14
        self.lib.lagmoment_stage.restype=C.c_int
        self.lib.lagmoment_factor.argtypes=[C.c_int]*3+[C.c_int64]+[C.c_void_p]*6+[C.c_double]+[C.c_void_p]*3
        self.lib.lagmoment_factor.restype=C.c_int

    def stage(self,data):
        keys=('node_inputs','cells','owners','offsets','incidence','frozen','extra','factors','gravity_nodes','gravity_offsets')
        arrays={k:(i64(data[k]) if k in ('owners','offsets','incidence','frozen','gravity_offsets') else f64(data[k])) for k in keys}
        n=len(arrays['cells']);nn=len(arrays['node_inputs']);ni=len(arrays['incidence']);ng=len(arrays['gravity_nodes'])
        shapes={'node_inputs':(nn,33),'cells':(n,19),'owners':(nn,2),'offsets':(n+1,),'incidence':(ni,),'frozen':(n,),
                'extra':(nn,20),'factors':(n,2,2),'gravity_nodes':(ng,7),'gravity_offsets':(n+1,)}
        for key,shape in shapes.items():
            if arrays[key].shape!=shape:raise ValueError((key,arrays[key].shape,shape))
        out=dict(stage=np.empty((n,94)),Q=np.empty((ni,16)),hydro=np.empty((n,106)),ledgers=np.empty(12))
        status=self.lib.lagmoment_stage(self.kind,self.threads,n,nn,ni,ng,*[ptr(arrays[k]) for k in keys],*[ptr(out[k]) for k in ('stage','Q','hydro','ledgers')])
        if status:raise RuntimeError(f'compiled stage failed with status {status}')
        return out

    def factor(self,base,linear,A,K,Q,offsets,h,raw,frozen):
        base=f64(base);linear=f64(linear);A=f64(A);K=f64(K);Q=f64(Q);offsets=i64(offsets);raw=f64(raw);frozen=i64(frozen)
        n=len(base)
        for key,a,shape in (('base',base,(n,2,2)),('linear',linear,(n,2,2)),('A',A,(n,2,2)),('K',K,(n,2,2)),('raw',raw,(n,18)),('frozen',frozen,(n,)),('offsets',offsets,(n+1,))):
            if a.shape!=shape:raise ValueError((key,a.shape,shape))
        if Q.ndim!=2 or Q.shape[1]!=2:raise ValueError('Q columns must be rows [column][2]')
        out=np.empty((n,63))
        status=self.lib.lagmoment_factor(self.kind,self.threads,n,len(Q),*[ptr(a) for a in (base,linear,A,K,Q,offsets)],h,*[ptr(a) for a in (raw,frozen,out)])
        if status:raise RuntimeError(f'compiled factor failed with status {status}')
        return out
