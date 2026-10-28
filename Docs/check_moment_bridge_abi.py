"""Bridge regression on completed captures and caller-output transaction tests."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_compiled_midpoint_modules_20261018';sys.path.insert(0,str(modules))
from moment_compiled_backend import Backend,ptr
from moment_bridge_audit import audit_stage,audit_factor
build=base/'Codex_compiled_midpoint_build_20261018';capture=base/'Codex_kinetic_face_capture_20261018';prep=base/'Codex_hydro_cell_backends_20261018';stage=base/'Codex_moment_stage_capture_20261018'
def binary(root,name,dtype,width=None):
    a=np.fromfile(root/name,dtype=dtype)
    return a.reshape(-1,width) if width else a
data=dict(node_inputs=binary(capture,'inputs.bin','<f8',33),cells=binary(prep,'cells.bin','<f8',19),owners=binary(prep,'owners.bin','<i8',2),offsets=binary(prep,'offsets.bin','<i8'),incidence=binary(prep,'incidence.bin','<i8'),frozen=binary(prep,'frozen.bin','<i8'),extra=binary(stage,'node_extra.bin','<f8',20),factors=binary(stage,'factors.bin','<f8',4).reshape(-1,2,2),gravity_nodes=binary(stage,'gravity_nodes.bin','<f8',7),gravity_offsets=binary(stage,'gravity_offsets.bin','<i8'))
st=np.load(stage/'reference.npz');expect=st['expected'];hy=np.load(prep/'reference.npz')['expected']
reference=dict(raw=expect[:,:18],unmasked_hydro=hy[:,:18],hydro=hy[:,18:36],pressure=hy[:,36:54],source=hy[:,78:96],F=hy[:,72:74],H=hy[:,74:78].reshape(-1,2,2),K=expect[:,22:26].reshape(-1,2,2),K_unmasked=expect[:,63:67].reshape(-1,2,2),Q=st['Q'],gravity=expect[:,30:48],B=expect[:,18:22].reshape(-1,2,2),G=expect[:,48:52].reshape(-1,2,2),Ag=expect[:,52:56].reshape(-1,2,2),cov_A=expect[:,56:60].reshape(-1,2,2))
backends={name:Backend(build/'liblagmoment.so',kind,threads) for name,kind,threads in [('cpu1','cpu',1),('omp4','cpu',4),('cuda','cuda',1)]}
outputs={name:b.stage(data) for name,b in backends.items()}
report=dict(stage=audit_stage(data,reference,outputs))
factor_root=base/'Codex_factor_update_backends_v2_20261018';n=len(data['cells'])
deck=binary(factor_root,'cpu1_input.bin','<f8',36)[:n];Q=binary(factor_root,'cpu1_Q.bin','<f8',2);off=binary(factor_root,'offsets.bin','<i8')
inputs=(deck[:,:4].reshape(-1,2,2),deck[:,4:8].reshape(-1,2,2),deck[:,8:12].reshape(-1,2,2),deck[:,12:16].reshape(-1,2,2),Q,off,float(deck[0,34]),deck[:,16:34],deck[:,35].astype(np.int64))
fout={name:b.factor(*inputs) for name,b in backends.items()};fref=np.load(factor_root/'cpu1_reference.npz')['expected'][0]
report['factor']=audit_factor(inputs,fref,fout)
# Raw C ABI checks: the caller buffers remain sentinel, including when cells
# before/after the failing cell have produced valid temporary outputs.
report['transactional']={}
for name,b in backends.items():
    bad={k:v.copy() for k,v in data.items()};i=int(np.flatnonzero(~data['frozen'].astype(bool))[0]);bad['cells'][i,3]=-1.
    out=[np.full((n,94),-947.),np.full((len(data['incidence']),16),-947.),np.full((n,106),-947.),np.full(12,-947.)]
    keys=('node_inputs','cells','owners','offsets','incidence','frozen','extra','factors','gravity_nodes','gravity_offsets')
    code=b.lib.lagmoment_stage(b.kind,b.threads,n,len(data['node_inputs']),len(data['incidence']),len(data['gravity_nodes']),*[ptr(bad[k]) for k in keys],*[ptr(a) for a in out])
    assert code>=100 and all(np.all(a==-947.) for a in out)
    baseL,linear,A,K,Q,off,h,raw,fr=[np.ascontiguousarray(x) if isinstance(x,np.ndarray) else x for x in inputs]
    linear=linear.copy();linear[i]=0.;outfactor=np.full((n,63),-947.)
    code2=b.lib.lagmoment_factor(b.kind,b.threads,n,len(Q),*[ptr(a) for a in (baseL,linear,A,K,Q,off)],h,*[ptr(a) for a in (raw,fr,outfactor)])
    assert code2>=100 and np.all(outfactor==-947.)
    report['transactional'][name]=dict(stage_status=code,factor_status=code2,outputs_unchanged=True)
report['library_sha256']=hashlib.sha256((build/'liblagmoment.so').read_bytes()).hexdigest()
(build/'abi_preflight.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
