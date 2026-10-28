"""C ABI for explicit host exact-resolution and CPU/CUDA sign filters."""
import ctypes as ct,json,time
from pathlib import Path
import numpy as np
D=ct.c_double;L=ct.c_int64;U=ct.c_uint64;I=ct.c_int;CAP=128
class Polygon(ct.Structure):
    _fields_=[('epoch',U),('count',I),('endpoint_rounding_difference',I),('endpoint_collapse',I),('endpoint_reversal',I),('max_relative_vector_error',D),('filtered_checks',U),('exact_checks',U),('added_candidates',U),('area',D),('centroid',D*2),('second',D*4)]+[(name,typ*CAP) for name,typ in [('x',D),('y',D),('support_a',L),('support_b',L),('neighbor',L),('neighbor_kind',I),('midpoint_x',D),('midpoint_y',D),('edge_dx',D),('edge_dy',D),('length',D),('length_lo',D),('length_hi',D),('nx',D),('ny',D),('geometry_wn',D)]]
def pointer(a,t):return a.ctypes.data_as(ct.POINTER(t))
class Backend:
    def __init__(self,library):
        self.lib=ct.CDLL(str(library));self.lib.moment_geometry_v2_size.restype=ct.c_size_t
        assert ct.sizeof(Polygon)==self.lib.moment_geometry_v2_size()
        self.lib.moment_geometry_v2_resolve.argtypes=[I,I,I,ct.POINTER(L),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(I),ct.POINTER(L),ct.POINTER(L),ct.POINTER(L),ct.POINTER(D),U,ct.POINTER(Polygon),ct.POINTER(I),ct.c_char_p]
        self.lib.moment_geometry_v2_filter.argtypes=[I,I,I,I,ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),I,ct.POINTER(L),ct.POINTER(D),ct.POINTER(ct.c_int8)]
        self.lib.moment_geometry_v2_exact_queries.argtypes=[I,ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),I,ct.POINTER(L),ct.POINTER(I)]
    def resolve(self,data,folder,threads):
        folder=Path(folder);folder.mkdir(exist_ok=False)
        arrays=[np.ascontiguousarray(data[k],dtype=d) for k,d in [('ids',np.int64),('xy',float),('W',float),('v',float),('z',float),('kind',np.int32),('selected',np.int64),('offsets',np.int64),('candidates',np.int64),('bounds',float)]]
        out=(Polygon*len(arrays[6]))();ct.memset(out,0xA5,ct.sizeof(out));before=bytes(out);status=np.full(len(out),-999,dtype=np.int32)
        pointers=[pointer(a,t) for a,t in zip(arrays,[L,D,D,D,D,I,L,L,L,D])];start=time.monotonic()
        code=self.lib.moment_geometry_v2_resolve(threads,len(out),len(arrays[0]),*pointers,int(data['epoch']),out,pointer(status,I),str(folder/'exact.json').encode());elapsed=time.monotonic()-start
        np.save(folder/'status.npy',status)
        if code:
            failure=dict(code=code,caller_unchanged=bytes(out)==before,seconds=elapsed);(folder/'failure.json').write_text(json.dumps(failure));raise RuntimeError(failure)
        result={name:np.array([getattr(p,name) for p in out]) for name,_ in Polygon._fields_}
        for i,n in enumerate(result['count']):
            for field,typ in Polygon._fields_:
                if getattr(typ,'_length_',0)==CAP:result[field][i,n:]=0
        assert np.all(status==0) and all(np.isfinite(a).all() for a in result.values())
        np.savez(folder/'export.npz',**result);return result,json.loads((folder/'exact.json').read_text()),elapsed
    def filter(self,data,queries,backend=0,threads=1,enclosures=None):
        xy=np.ascontiguousarray(data['xy'],dtype=float);W=np.ascontiguousarray(data['W'],dtype=float);bounds=np.ascontiguousarray(data['bounds'],dtype=float);queries=np.ascontiguousarray(queries,dtype=np.int64)
        mode=1 if enclosures is not None else 0
        if mode:
            assert queries.ndim==1;enclosures=np.ascontiguousarray(enclosures,dtype=float);assert enclosures.shape==(len(queries),4);out=np.empty((len(queries),len(xy)),dtype=np.int8)
        else:assert queries.ndim==2 and queries.shape[1]==4;out=np.empty(len(queries),dtype=np.int8)
        code=self.lib.moment_geometry_v2_filter(backend,threads,mode,len(xy),pointer(xy,D),pointer(W,D),pointer(bounds,D),len(queries),pointer(queries,L),pointer(enclosures,D) if mode else None,pointer(out,ct.c_int8))
        if code:raise RuntimeError(('filter',code))
        return out
    def exact(self,data,queries):
        xy=np.ascontiguousarray(data['xy'],dtype=float);W=np.ascontiguousarray(data['W'],dtype=float);bounds=np.ascontiguousarray(data['bounds'],dtype=float);queries=np.ascontiguousarray(queries,dtype=np.int64);assert queries.ndim==2 and queries.shape[1]==4;out=np.empty(len(queries),dtype=np.int32)
        code=self.lib.moment_geometry_v2_exact_queries(len(xy),pointer(xy,D),pointer(W,D),pointer(bounds,D),len(queries),pointer(queries,L),pointer(out,I))
        if code:raise RuntimeError(('exact signs',code))
        return out
