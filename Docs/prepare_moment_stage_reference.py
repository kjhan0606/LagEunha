"""Independent fixed-stage reference from captured coordinates and quadrature."""
import json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_moment_stage_capture_20261018'
sys.path[:0]=[str(base/'Codex_kinetic_face_capture_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018')]
from moment_gravity_source import covariance_source,covariance_generator
inputs=np.load(root/'cell_inputs.npz');raw=inputs['raw'];L=inputs['factors'];N=len(raw)
a=np.fromfile(base/'Codex_kinetic_face_capture_20261018/inputs.bin',dtype='<f8').reshape(-1,33)
e=np.fromfile(root/'node_extra.bin',dtype='<f8').reshape(-1,20)
prep=base/'Codex_hydro_cell_backends_20261018'
off=np.fromfile(prep/'offsets.bin',dtype='<i8');inc=np.fromfile(prep/'incidence.bin',dtype='<i8')
K=np.zeros((N,2,2));Q=np.zeros((len(inc),4,2));incoming=np.zeros_like(K)
Kabs=np.zeros_like(K);Qabs=np.zeros_like(Q);Gabs=np.zeros_like(K)
for i in range(N):
    for t in range(off[i],off[i+1]):
        code=int(inc[t]);k=abs(code)-1;side=0 if code<0 else 1
        row=a[k];ext=e[k];xi=ext[8*side:8*side+8].reshape(4,2)
        loss=row[21:25] if side==0 else -row[25:29]
        rates=-row[4]*row[25:29] if side==0 else row[4]*row[21:25]
        velocity=row[13:21].reshape(4,2) if side==0 else row[5:13].reshape(4,2)
        affine=ext[16+2*side:18+2*side]
        K[i]+=row[4]/raw[i,0]*(xi.T@(loss[:,None]*xi))
        Kabs[i]+=row[4]/raw[i,0]*(abs(xi).T@(loss[:,None]*abs(xi)))
        Q[t]=np.sqrt(rates)[:,None]*(velocity-affine)
        Qabs[t]=np.sqrt(rates)[:,None]*(abs(velocity)+abs(affine))
        incoming[i]+=Q[t].T@Q[t]
Kmasked=K.copy();Qmasked=Q.copy();incmasked=incoming.copy();frozen=inputs['frozen']
Kmasked[frozen]=0.;incmasked[frozen]=0.
for i in np.flatnonzero(frozen):Qmasked[off[i]:off[i+1]]=0.
grav=inputs['gravity'];B=inputs['affine_B'];G=[];Ag=[]
for r,g,l in zip(raw,grav,L):
    rate=covariance_source(r,g);G.append(rate);Ag.append(covariance_generator(l,rate))
G=np.array(G);Ag=np.array(Ag);covA=-B+Ag;covAmasked=covA.copy();covAmasked[frozen]=0.
hydro=np.load(prep/'reference.npz')['expected']
rawtotal=hydro[:,18:36].copy();rawtotal[~frozen]+=grav[~frozen]
rawunmask=hydro[:,:18]+grav
expected=np.column_stack([rawtotal,B.reshape(-1,4),Kmasked.reshape(-1,4),incmasked.reshape(-1,4),grav,G.reshape(-1,4),Ag.reshape(-1,4),covAmasked.reshape(-1,4),np.zeros((N,3)),K.reshape(-1,4),incoming.reshape(-1,4),covA.reshape(-1,4),rawunmask,np.zeros(N)])
Qcombined=np.concatenate([Q.reshape(-1,8),Qmasked.reshape(-1,8)],axis=1)
np.savez(root/'reference.npz',expected=expected,Q=Qcombined,Kabsolute=Kabs,Qabsolute=Qabs)
activeg=grav[~frozen];ledger=np.r_[activeg[:,7:9].sum(0),activeg[:,17].sum(),(activeg[:,10]-activeg[:,11]).sum()]
np.save(root/'gravity_ledger.npy',ledger)
R=L@L.transpose(0,2,1)
print(json.dumps(dict(cells=N,incidences=len(inc),incoming_columns=4*len(inc),max_R_condition=float(np.max(np.linalg.cond(R))),max_Ag=float(abs(Ag).max()),max_G=float(abs(G).max()),gravity_thermal_max=float(abs(grav[:,17]-.5*(grav[:,13]+grav[:,16])).max()),gravity_origin_torque_max=float(abs(grav[:,10]-grav[:,11]).max()),scope=__doc__),indent=2))
