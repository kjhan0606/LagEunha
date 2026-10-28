"""Reuse completed node capture and raw state; no RHS evaluation or evolution."""
import hashlib,json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_hydro_cell_backends_20261018';root.mkdir(exist_ok=False)
capture=base/'Codex_kinetic_face_capture_20261018'
a=np.fromfile(capture/'inputs.bin',dtype='<f8').reshape(-1,33)
ff=np.fromfile(capture/'expected.bin',dtype='<f8').reshape(-1,18)
owners=np.load(capture/'owners.npy');N=1590
original=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
complete=np.load(base/'Codex_contact_prior_complete_rhs_20261018/rhs.npz')
assert np.array_equal(original['ids'],complete['ids'])
state=original['state'][:-4].reshape(-1,20)
raw=np.column_stack([original['mass'],state[:,3:]])
frozen=original['frozen'];c=raw[:,1:3]/raw[:,0,None]
gamma=float(np.float32(1.666667))
pv=(gamma-1)*(raw[:,17]-.5*(raw[:,13]+raw[:,16]))
rows=[[] for _ in range(N)]
dr=np.zeros((N,18));F=np.zeros((N,2));H=np.zeros((N,2,2));work=np.zeros(N);ext=np.zeros(4)
absface=np.zeros((N,18));absF=np.zeros((N,2));absH=np.zeros((N,2,2))
for k,((i,j),v,flux) in enumerate(zip(owners,a,ff)):
    traction=v[4]*v[31]*v[2:4]
    for receiver,sign in ((i,-1),(j,1)):
        if receiver<0:continue
        rows[receiver].append(sign*(k+1));dr[receiver]+=sign*flux
        F[receiver]+=sign*traction
        term=np.outer(v[:2]-c[receiver],traction)
        H[receiver]+=sign*term
        work[receiver]+=sign*v[4]*v[31]*v[32]
        absface[receiver]+=abs(flux);absF[receiver]+=abs(traction);absH[receiver]+=abs(term)
    if j<0:ext-=np.r_[flux[7:9],flux[17],flux[10]-flux[11]]
src=np.zeros((N,18));diag=np.zeros((N,8));pressure=np.zeros((N,18));S=np.zeros((N,2,2))
for i,r in enumerate(raw):
    p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2)
    S[i]=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/r[0]
    h=H[i]+pv[i]*np.eye(2)
    A=np.linalg.solve(S[i],h).T;a0=F[i]/r[0]-A@c[i]
    src[i,1:3]=p;src[i,3:7]=(D+D.T).ravel();src[i,9:13]=(Pi+pv[i]*np.eye(2)).ravel()
    src[i,13:17]=(np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T).ravel()
    pressure[i,7:9]=F[i];pressure[i,9:13]=(h+np.outer(c[i],F[i])).ravel()
    pressure[i,13:17]=src[i,13:17];pressure[i,17]=work[i]
    diag[i,:4]=A.ravel();diag[i,4:6]=a0;diag[i,7]=np.linalg.cond(S[i],np.inf)
unmasked=dr+src
frozenledger=-np.r_[unmasked[frozen,7:9].sum(axis=0),unmasked[frozen,17].sum(),(unmasked[frozen,10]-unmasked[frozen,11]).sum()]
hydro=unmasked.copy();hydro[frozen]=0.;pressure[frozen]=0.
legacy=np.column_stack([np.zeros(N),complete['rates'][:-4].reshape(-1,20)[:,3:]])
diff=float(abs(hydro[:,1:]-legacy[:,1:]).max())
pdiff=float(abs(pressure-complete['pressure_raw_rates']).max())
assert diff<1e-15 and pdiff<1e-15
expected=np.column_stack([unmasked,hydro,pressure,hydro-pressure,F,H.reshape(-1,4),src,diag,work,np.zeros(N)])
offsets=np.r_[0,np.cumsum([len(r) for r in rows])].astype('<i8')
arrays={'cells.bin':np.column_stack([raw,pv]),'owners.bin':owners.astype('<i8'),'offsets.bin':offsets,
        'incidence.bin':np.concatenate(rows).astype('<i8'),'frozen.bin':frozen.astype('<i8')}
for name,array in arrays.items():array.tofile(root/name)
np.savez(root/'reference.npz',expected=expected,ledger=np.r_[ext,frozenledger],raw=raw,pv=pv,
    spatial_covariance=S,absolute_face=absface,absolute_F=absF,absolute_H=absH,incidence_counts=np.diff(offsets))
report=dict(nodes=len(a),cells=N,gamma=gamma,PV_definition='(gamma-1)*(Ehyd-tr(Pi)/2), identical to P*area',
    full_hydro_replay_max=diff,pressure_replay_max=pdiff,maximum_S_condition=float(diag[:,7].max()),
    input_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in arrays},
    scope='Reuses captured inputs and completed raw state without new Python RHS, massroot or evolution.')
(root/'manifest.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
