"""Unchanged conditioning/QR checks on native-owned tentative factor results.

No updates, state replacement, repair, GPU calls, MPI or commitment.
"""
import json,time
from pathlib import Path
import numpy as np
MAGIC=0x31544341464d464c

def evaluate(service,request):
    from moment_factor_error_bounds import make_bounds,check
    root=Path(bytes(request.output_dir).decode()).resolve()
    assert root.parent==service.root,'Foreign owned factor audit directory'
    header=np.fromfile(root/'header.bin',dtype='<u8');assert header.shape==(6,)
    magic,epoch,token,n,columns,rank=map(int,header)
    assert magic==MAGIC and epoch==request.epoch and token==request.run_token and rank==request.rank
    assert n==service.ctx['hi']-service.ctx['lo'] and columns>0
    def load(name,shape,dtype='<f8'):
        path=root/(name+'.bin');assert path.stat().st_size==int(np.prod(shape))*np.dtype(dtype).itemsize
        return np.memmap(path,mode='r',dtype=dtype,shape=shape)
    deck=load('deck',(n,36));out=load('factor',(n,63));Q=load('incoming_Q',(columns,2));offset=load('incoming_offsets',(n+1,),'<i8')
    assert offset[0]==0 and offset[-1]==columns and np.all(np.diff(offset)>=0)
    assert np.isfinite(deck).all() and np.isfinite(out).all() and np.isfinite(Q).all()
    assert (deck[:,34]>0).all() and (deck[:,34]<=.005).all() and np.all(deck[:,34]==deck[0,34]) and not deck[:,35].any()
    raw=deck[:,16:34]
    assert not out[:,62].any()
    for k in list(range(13))+[17]:assert np.array_equal(out[:,4+k].view(np.uint64),raw[:,k].view(np.uint64))
    assert np.array_equal(out[:,18].view(np.uint64),out[:,19].view(np.uint64))
    assert (out[:,0]>0).all() and (out[:,3]>0).all() and not out[:,1].any()
    assert (out[:,51]>0).all() and out[:,53].max()<128*np.finfo(float).eps
    started=time.perf_counter();tol,metrics=make_bounds(deck,out,Q,offset,n);maximum=0.
    for i,metric in enumerate(metrics):
        L=out[i,:4].reshape(2,2);U=metric['U'];thin=np.linalg.solve(L,U).T
        bound=64*metric['qr_gamma']*metric['condition_U']**2
        assert np.linalg.norm(thin.T@thin-np.eye(2),np.inf)<bound+64*metric['exp_bound']
        assert out[i,58]<bound
        check('UUT',L@L.T,U@U.T,2*tol[i,38:42].reshape(2,2))
        maximum=max(maximum,float(out[i,58]/bound))
    shadow=None
    if request.outer:
        from moment_bridge_audit import factor_reference,audit_factor
        from moment_compiled_backend import Backend
        inputs=(deck[:,:4].reshape(n,2,2),deck[:,4:8].reshape(n,2,2),deck[:,8:12].reshape(n,2,2),deck[:,12:16].reshape(n,2,2),Q,offset,float(deck[0,34]),raw,np.zeros(n,dtype=np.int64))
        reference=factor_reference(*inputs)
        cpu=Backend('/gpfs/kjhan/LagForce/Codex_full_disk_owned_ledger_build_v4_20261019/liblagmoment_owned.so','cpu',4).factor(*inputs)
        shadow=audit_factor(inputs,reference,dict(cpu4=cpu,cuda=out))
    result=dict(status='passed_owned_native_factor_audit',epoch=epoch,rank=rank,cells=n,
        seconds=time.perf_counter()-started,minimum_thermal=float(out[:,51].min()),
        maximum_QR_bound_fraction=maximum,shadow_audit=shadow,service_GPU_calls=0,service_MPI_calls=0,accepted_time_steps=0)
    path=root/'audit.json'
    with path.open('x') as f:json.dump(result,f,indent=2)
    return path
