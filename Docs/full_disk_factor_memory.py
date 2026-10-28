"""Existing midpoint factor map with conditioning/positivity checks per update.

Uses the validated CUDA ABI. Independent Python/CPU/CUDA shadows belong to
checkpoint audits, while the same QR conditioning bounds remain active here.
"""
import hashlib,time
import numpy as np
from moment_factor_error_bounds import make_bounds,check

def update(base,linear,saved,h,time_value,backend,shadow_audit=False):
    started=time.perf_counter();assert 0<h<=.005 and time_value>float(base['time'])
    mass,ids=base['mass'],base['ids'];n=len(ids)
    assert np.array_equal(linear['ids'],ids) and np.array_equal(linear['mass'],mass)
    assert np.array_equal(saved['ids'],ids) and np.array_equal(saved['mass'],mass)
    signature=hashlib.sha256(np.ascontiguousarray(linear['state']).view(np.uint8)).hexdigest()
    assert str(saved['initial_state_sha256'])==signature
    a=base['state'][:-4].reshape(n,20)
    raw=np.column_stack([mass,a[:,3:]])
    stage=saved['stage'];L=base['factors'];Ll=linear['factors']
    A=stage[:,56:60].reshape(n,2,2);K=stage[:,22:26].reshape(n,2,2)
    Q,offset=saved['incoming_Q'],saved['incoming_offsets']
    provisional=raw+h*stage[:,:18];provisional[:,0]=mass
    frozen=np.zeros(n,dtype=np.int64)
    out=backend.factor(L,Ll,A,K,Q,offset,h,provisional,frozen)
    assert np.isfinite(out).all() and not out[:,62].any()
    for k in list(range(13))+[17]:assert np.array_equal(out[:,4+k].view(np.uint64),provisional[:,k].view(np.uint64))
    assert np.array_equal(out[:,18].view(np.uint64),out[:,19].view(np.uint64))
    assert (out[:,0]>0).all() and (out[:,3]>0).all() and not out[:,1].any()
    assert (out[:,51]>0).all() and out[:,53].max()<128*np.finfo(float).eps
    deck=np.column_stack([L.reshape(n,4),Ll.reshape(n,4),A.reshape(n,4),K.reshape(n,4),provisional,np.full(n,h),frozen])
    tol,metrics=make_bounds(deck,out,Q,offset,n)
    maximum_QR_fraction=0.
    for i,metric in enumerate(metrics):
        newL=out[i,:4].reshape(2,2);U=metric['U'];thin=np.linalg.solve(newL,U).T
        bound=64*metric['qr_gamma']*metric['condition_U']**2
        assert np.linalg.norm(thin.T@thin-np.eye(2),np.inf)<bound+64*metric['exp_bound']
        assert out[i,58]<bound
        check('UUT',newL@newL.T,U@U.T,2*tol[i,38:42].reshape(2,2))
        maximum_QR_fraction=max(maximum_QR_fraction,float(out[i,58]/bound))
    state=base['state'].copy();b=state[:-4].reshape(n,20)
    b[:,:2]=a[:,:2]+h*saved['generator'];b[:,2]=a[:,2]+h*saved['weight_rates']
    b[:,3:]=out[:,5:22]
    state[-4:]=base['state'][-4:]+h*saved['ledgers'].reshape(3,4).sum(0)
    Psi=base['Psi'].astype(np.longdouble)+np.longdouble(h)*saved['Psi_rates'].astype(np.longdouble)
    assert np.isfinite(state).all() and np.isfinite(Psi).all()
    backend_audit=None
    if shadow_audit:
        from moment_bridge_audit import factor_reference,audit_factor
        from moment_compiled_backend import Backend
        inputs=(L,Ll,A,K,Q,offset,h,provisional,frozen)
        reference=factor_reference(*inputs)
        cpu=Backend(backend.path,'cpu',4).factor(*inputs)
        backend_audit=audit_factor(inputs,reference,dict(cpu4=cpu,cuda=out))
    result=dict(state=state,factors=out[:,:4].reshape(n,2,2).copy(),mass=mass,ids=ids,Psi=Psi,
        time=time_value,frozen=np.zeros(n,bool),psi_constrained_mask=np.ones(n,bool),weight_rate_seed=saved['weight_rates'].copy())
    report=dict(status='passed_production_factor_update',time=time_value,h=h,
        seconds=time.perf_counter()-started,minimum_thermal=float(out[:,51].min()),
        minimum_factor_diagonal=float(out[:,[0,3]].min()),maximum_QR_orthogonality=float(out[:,58].max()),
        maximum_QR_bound_fraction=maximum_QR_fraction,shadow_backends=2 if shadow_audit else 0,backend_audit=backend_audit)
    return result,report
