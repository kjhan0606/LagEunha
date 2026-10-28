"""Typed owned-cell ABI; global node identity and physical roles stay intact."""
import ctypes as C
import numpy as np
from moment_compiled_backend import Backend,f64,i64,ptr

class OwnedBackend(Backend):
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw)
        self.lib.lagmoment_canonical_nodes.argtypes=[C.c_int]*3+[C.c_void_p]*2
        self.lib.lagmoment_canonical_nodes.restype=C.c_int
        self.lib.lagmoment_owned_stage.argtypes=[C.c_int]*2+[C.c_uint64]*2+[C.c_int]*3+[C.c_int64]*2+[C.c_void_p]*16
        self.lib.lagmoment_owned_stage.restype=C.c_int

    def ordered_ledgers(self,data,output,node_flux,sentinel=None):
        function=self.lib.lagmoment_ordered_ledgers
        function.argtypes=[C.c_int]*4+[C.c_void_p]*6;function.restype=C.c_int
        nf=f64(node_flux);owners=i64(data['owners']);frozen=i64(data['frozen'])
        stage=f64(output['stage']);hydro=f64(output['hydro']);n=len(stage);nn=len(nf)
        assert nf.shape==(nn,19) and owners.shape==(nn,2) and frozen.shape==(n,)
        assert stage.shape==(n,94) and hydro.shape==(n,106)
        out=np.empty(12)
        if sentinel is not None:out.fill(sentinel)
        status=function(self.kind,self.threads,n,nn,*map(ptr,(nf,owners,frozen,stage,hydro,out)))
        if sentinel is not None:return status,out
        if status:raise RuntimeError(('ordered_ledger_failure',status))
        return out

    def nodes(self,node_inputs):
        a=f64(node_inputs);assert a.ndim==2 and a.shape[1]==33
        out=np.empty((len(a),19))
        status=self.lib.lagmoment_canonical_nodes(self.kind,self.threads,len(a),ptr(a),ptr(out))
        if status:raise RuntimeError(('canonical_node_failure',status))
        return out

    def owned_stage(self,registry_count,data,epoch=1,expected_epoch=1,sentinel=None):
        keys=('owned_registry_index','node_global_index','node_inputs','node_flux','cells','owners','offsets','incidence',
              'frozen','extra','factors','gravity_nodes','gravity_offsets')
        integers={'owned_registry_index','node_global_index','owners','offsets','incidence','frozen','gravity_offsets'}
        arrays={k:(i64(data[k]) if k in integers else f64(data[k])) for k in keys}
        n,nn,ni,ng=len(arrays['cells']),len(arrays['node_inputs']),len(arrays['incidence']),len(arrays['gravity_nodes'])
        shapes=dict(owned_registry_index=(n,),node_global_index=(nn,),node_inputs=(nn,33),node_flux=(nn,19),cells=(n,19),
            owners=(nn,2),offsets=(n+1,),incidence=(ni,),frozen=(n,),extra=(nn,20),factors=(n,2,2),gravity_nodes=(ng,7),gravity_offsets=(n+1,))
        for key,shape in shapes.items():assert arrays[key].shape==shape,(key,arrays[key].shape,shape)
        out=dict(stage=np.empty((n,94)),Q=np.empty((ni,16)),hydro=np.empty((n,106)))
        if sentinel is not None:
            for value in out.values():value.fill(sentinel)
        status=self.lib.lagmoment_owned_stage(self.kind,self.threads,epoch,expected_epoch,registry_count,n,nn,ni,ng,
            *[ptr(arrays[k]) for k in keys],*[ptr(out[k]) for k in ('stage','Q','hydro')])
        if sentinel is not None:return status,out
        if status:raise RuntimeError(('owned_stage_failure',status))
        return out

def partition(data,nf,owned):
    """Preserve each cell's canonical node order without constructing walls."""
    owned=np.asarray(owned,dtype=np.int64);off=data['offsets'];goff=data['gravity_offsets']
    incidence=np.concatenate([data['incidence'][off[i]:off[i+1]] for i in owned])
    nodes=np.unique(abs(incidence)-1)
    local_code=np.sign(incidence)*(np.searchsorted(nodes,abs(incidence)-1)+1)
    return dict(owned_registry_index=owned,node_global_index=nodes,node_inputs=data['node_inputs'][nodes],node_flux=nf[nodes],
        cells=data['cells'][owned],owners=data['owners'][nodes],offsets=np.r_[0,np.cumsum(np.diff(off)[owned])],incidence=local_code,
        frozen=data['frozen'][owned],extra=data['extra'][nodes],factors=data['factors'][owned],
        gravity_nodes=np.concatenate([data['gravity_nodes'][goff[i]:goff[i+1]] for i in owned]),
        gravity_offsets=np.r_[0,np.cumsum(np.diff(goff)[owned])])
