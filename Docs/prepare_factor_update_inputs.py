"""Frozen-coefficient factor maps, not new midpoint stages or evolution."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
base=Path('/gpfs/kjhan/LagForce');stage=base/'Codex_moment_stage_capture_20261018'
root=base/'Codex_factor_update_backends_v2_20261018';root.mkdir(exist_ok=False)
state=np.load(stage/'cell_inputs.npz');raw=state['raw'];L=state['factors'];frozen=state['frozen'];N=len(raw)
offsets=np.fromfile(base/'Codex_hydro_cell_backends_20261018/offsets.bin',dtype='<i8')*4
offsets.tofile(root/'offsets.bin')
references={};manifest=dict(cells=N,heights=[5e-5,1e-4],common_input='verified CPU1 stage raw/cov_A/K/Qfactor, captured Lbase=Llinear',scope=__doc__,inputs={})
for name in ('cpu1','omp4','cuda'):
    st=np.fromfile(stage/(name+'.bin'),dtype='<f8').reshape(N,94)
    Q=np.fromfile(stage/(name+'_Q.bin'),dtype='<f8').reshape(-1,16)[:,8:].reshape(-1,2).copy()
    Q.tofile(root/(name+'_Q.bin'))
    decks=[];expected=[];extra=[]
    for h in manifest['heights']:
        provisional=raw+h*st[:,:18]
        provisional[:,0]=raw[:,0] # Actual flat17 evolution holds mass separately.
        provisional[frozen]=raw[frozen]
        deck=np.column_stack([L.reshape(N,4),L.reshape(N,4),st[:,56:60],st[:,22:26],provisional,np.full(N,h),frozen.astype(float)])
        assert np.isfinite(deck).all() and np.array_equal(deck[:,13],deck[:,14])
        decks.append(deck);out=np.zeros((N,63));extras=np.zeros((N,3))
        for i in range(N):
            out[i,:4]=L[i].ravel();out[i,4:22]=provisional[i]
            if frozen[i]:continue
            baseL=L[i];linear=L[i];A=st[i,56:60].reshape(2,2);K=st[i,22:26].reshape(2,2)
            U,s,Vh=np.linalg.svd(linear);scaled=s/s.max()
            ratio=scaled[:,None]*scaled[None,:]/(scaled[:,None]**2+scaled[None,:]**2)
            loss=U@((Vh@K@Vh.T)*ratio)@U.T;Ac=A-loss
            eh=expm(.5*h*Ac);ef=expm(h*Ac)
            qi=Q[offsets[i]:offsets[i+1]].T
            columns=np.column_stack([ef@baseL,np.sqrt(h)*(eh@qi)])
            _,r=np.linalg.qr(columns.T,mode='reduced')
            signs=np.where(np.diag(r)<0,-1.,1.);newL=(signs[:,None]*r).T
            rn=newL@newL.T;pr=provisional[i];m=pr[0];q=pr[1:3];p=pr[7:9]
            S=pr[3:7].reshape(2,2)-np.outer(q,q)/m;C=pr[9:13].reshape(2,2)-np.outer(q,p)/m
            B=np.linalg.solve(S,C).T;pi=np.outer(p,p)/m+B@S@B.T+rn
            pi[1,0]=pi[0,1]
            rebuilt=pr.copy();rebuilt[13:17]=pi.ravel()
            before=pr[17]-.5*(pr[13]+pr[16]);after=pr[17]-.5*np.trace(pi)
            assert after>0
            out[i,:4]=newL.ravel();out[i,4:22]=rebuilt
            d=out[i,22:62];d[:4]=loss.ravel();d[4:8]=Ac.ravel();d[8:12]=eh.ravel();d[12:16]=ef.ravel();d[16:20]=rn.ravel();d[20:24]=B.ravel();d[24:28]=(pi-pr[13:17].reshape(2,2)).ravel();d[28:31]=before,after,after-before
            gram=columns@columns.T;d[39]=np.linalg.norm(gram,np.inf)
            extras[i]=np.linalg.cond(S,np.inf),np.linalg.cond(linear@linear.T,np.inf),len(columns.T)
        expected.append(out);extra.append(extras)
    arrays=np.stack(decks);arrays.tofile(root/(name+'_input.bin'))
    np.savez(root/(name+'_reference.npz'),expected=np.stack(expected),extra=np.stack(extra))
    manifest['inputs'][name]={key:hashlib.sha256((root/(name+key)).read_bytes()).hexdigest() for key in ('_input.bin','_Q.bin','_reference.npz')}
manifest['common_files']=['cpu1_input.bin','cpu1_Q.bin']
manifest['provisional_raw']='m exact raw_initial m; other raw moments initial+h*same_backend_stage_rate; frozen exact raw_initial. Fixed coefficients only, not an actual evaluated midpoint timestep.'
manifest['unapplied_h_times_dm_maximum']={name:{str(h):float(h*abs(np.fromfile(stage/(name+'.bin'),dtype='<f8').reshape(N,94)[:,0]).max()) for h in manifest['heights']} for name in ('cpu1','omp4','cuda')}
manifest['Q_layout']='column-major physical vectors stored as rows [column][2], original CSR node/sample order, frozen columns zero'
(root/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
