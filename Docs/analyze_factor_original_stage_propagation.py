"""Propagate the ORIGINAL Python-stage discrepancy through the same factor map.

This is independent of CPU-vs-CUDA differences, and evaluates no new RHS.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
base=Path('/gpfs/kjhan/LagForce');stage=base/'Codex_moment_stage_capture_20261018'
root=base/'Codex_factor_update_backends_v2_20261018'
state=np.load(stage/'cell_inputs.npz');raw=state['raw'];L=state['factors'];frozen=state['frozen'];active=~frozen;mass=raw[:,0];N=len(raw)
original=np.load(stage/'reference.npz');st=original['expected'];Q=original['Q'][:,8:].reshape(-1,2)
offset=np.fromfile(root/'offsets.bin',dtype='<i8');cpu_ref=np.load(root/'cpu1_reference.npz')['expected']
cpu_stage=np.fromfile(stage/'cpu1.bin',dtype='<f8').reshape(N,94)
maps=[];rows=[]
for k,h in enumerate((5e-5,1e-4)):
    provisional=raw+h*st[:,:18];provisional[:,0]=raw[:,0];provisional[frozen]=raw[frozen]
    newL=L.copy();rebuilt=provisional.copy();R=L@L.transpose(0,2,1);thermal=np.zeros(N)
    for i in range(N):
        if frozen[i]:continue
        U,s,Vh=np.linalg.svd(L[i]);scaled=s/s.max()
        ratio=scaled[:,None]*scaled[None,:]/(scaled[:,None]**2+scaled[None,:]**2)
        K=st[i,22:26].reshape(2,2)
        loss=U@((Vh@K@Vh.T)*ratio)@U.T
        Ac=st[i,56:60].reshape(2,2)-loss
        columns=np.column_stack([expm(h*Ac)@L[i],np.sqrt(h)*(expm(.5*h*Ac)@Q[offset[i]:offset[i+1]].T)])
        _,r=np.linalg.qr(columns.T,mode='reduced');r=np.where(np.diag(r)<0,-1.,1.)[:,None]*r
        newL[i]=r.T;R[i]=newL[i]@newL[i].T
        rr=provisional[i];m=rr[0];q=rr[1:3];p=rr[7:9]
        S=rr[3:7].reshape(2,2)-np.outer(q,q)/m;C=rr[9:13].reshape(2,2)-np.outer(q,p)/m
        B=np.linalg.solve(S,C).T;Pi=np.outer(p,p)/m+B@S@B.T+R[i];Pi[1,0]=Pi[0,1]
        rebuilt[i,13:17]=Pi.ravel();thermal[i]=rr[17]-.5*np.trace(Pi)
        assert thermal[i]>0
    cpu=cpu_ref[k];dR=cpu[:,38:42].reshape(-1,2,2)-R;dPi=cpu[:,17:21]-rebuilt[:,13:17]
    dthermal=cpu[:,51]-thermal
    rows.append(dict(h=h,maximum_original_Python_stage_Ag_difference=float(abs(cpu_stage[:,52:56]-st[:,52:56]).max()),h_times_maximum_Ag_difference=float(h*abs(cpu_stage[:,52:56]-st[:,52:56]).max()),
        maximum_factor_difference=float(abs(cpu[active,:4]-newL[active].reshape(-1,4)).max()),
        maximum_R_difference=float(abs(dR[active]).max()),maximum_relative_R_Frobenius_difference=float(np.max(np.linalg.norm(dR[active],axis=(1,2))/np.linalg.norm(R[active],axis=(1,2)))),
        maximum_Pi_difference=float(abs(dPi[active]).max()),maximum_thermal_difference=float(abs(dthermal[active]).max()),
        maximum_specific_thermal_difference=float(np.max(abs(dthermal[active]/mass[active]))),
        mass_weighted_specific_thermal_RMS=float(np.sqrt(np.average((dthermal[active]/mass[active])**2,weights=mass[active]))),
        maximum_relative_thermal_difference=float(np.max(abs(dthermal[active]/thermal[active]))),
        unapplied_h_times_dm_maximum=float(h*abs(st[:,0]).max())))
    maps.append(dict(factors=newL,raw=rebuilt,R=R,thermal=thermal))
np.savez(root/'original_python_stage_factor_maps.npz',**{f'{key}_{k}':value for k,row in enumerate(maps) for key,value in row.items()})
report=dict(scope=__doc__,rows=rows,mass_contract='Fixed initial mass in both original-Python-stage and compiled-stage maps.',
    stage_reference_K='Original Python K retained, including normal last-bit asymmetry; its existing SVD loss formula is used without projection.',
    limits='Frozen-coefficient local discrepancy propagation, not actual midpoint evolution, cumulative-error bound, temporal order or long Kepler stability.',
    stage_reference_sha256=hashlib.sha256((stage/'reference.npz').read_bytes()).hexdigest(),cpu_map_reference_sha256=hashlib.sha256((root/'cpu1_reference.npz').read_bytes()).hexdigest())
(root/'original_stage_propagation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
