"""Persistent owner-vector residual: static actual-state and support-swap audit.

This tests a new operator, not the production HG operator. Support is copied
from saved HG neighborhoods; smooth changing support and evolution remain open.
"""
import json,struct
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive')
s=np.loadtxt(root/'raw_faces.dat.hgstate.rank0000');s=s[np.argsort(s[:,0])]
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
x=s[:,3:5]-3;m=s[:,5];U=s[:,8];volume=rhs[:,6];n=len(s)
assert np.array_equal(s[:,0],np.arange(n))
path=root/'raw_faces.dat.hgm.rank0000'
with path.open('rb') as f:magic,nm,ns,ncol=struct.unpack('<4q',f.read(32))
assert magic==20261017 and ncol==9
modes=np.memmap(path,mode='r',dtype='<f8',offset=32,shape=(nm,9))
support=np.memmap(path,mode='r',dtype='<f8',offset=32+nm*72,shape=(ns,3))
neighbors={}
for row in modes:
    owner=int(row[0])
    if np.linalg.norm(x[owner])>=.6 or owner in neighbors:continue
    start,length=map(int,row[2:4]);ids=support[start:start+length,0].astype(int)
    assert owner in ids and np.all(np.linalg.norm(x[ids],axis=1)<2.5)
    neighbors[owner]=ids
gamma=float(np.float32(1.666667));rng=np.random.default_rng(20261018)
z=rng.normal(size=(n,2));v=rng.normal(size=(n,2))
owners=np.array(sorted(neighbors));energy=float(.5*np.sum(z[owners]**2))
results=[]
for branch in ('original_support','drop_farthest'):
    force=np.zeros((n,2));zdot=np.zeros((n,2));maxmoment=0.;maxcond=0.;minimum_residual=1.
    for owner in owners:
        ids=neighbors[owner].copy();h=np.sqrt(volume[owner])
        if branch=='drop_farthest':
            farthest=np.argmax(np.linalg.norm(x[ids]-x[owner],axis=1))
            assert ids[farthest]!=owner
            ids=np.delete(ids,farthest)
        assert len(ids)>6
        dx=(x[ids]-x[owner])/h
        P=np.column_stack((np.ones(len(ids)),dx,dx[:,0]**2,dx[:,0]*dx[:,1],dx[:,1]**2))
        u,sv,_=np.linalg.svd(P,full_matrices=False)
        assert sv[-1]>1e-10*sv[0]
        maxcond=max(maxcond,float(sv[0]/sv[-1]))
        index=int(np.where(ids==owner)[0][0]);r=-u@u[index];r[index]+=1
        minimum_residual=min(minimum_residual,float(np.linalg.norm(r)))
        maxmoment=max(maxmoment,float(np.max(abs(P.T@r))/max(np.max(abs(P).T@abs(r)),1e-300)))
        # B has sqrt(mass)/time units; no division by small residual norms.
        sound=np.sqrt(gamma*(gamma-1)*U[owner]/m[owner])
        b=np.sqrt(m[owner])*sound/h*r
        zdot[owner]=b@v[ids]
        np.add.at(force,ids,-b[:,None]*z[owner])
    work=float(np.sum(force*v));aux=float(np.sum(z*zdot));den=np.sum(abs(force*v))+np.sum(abs(z*zdot))
    fnorm=np.linalg.norm(force,axis=1).sum()
    torque=float(np.sum(x[:,0]*force[:,1]-x[:,1]*force[:,0]))
    result={'branch':branch,'owners':len(owners),'auxiliary_energy':energy,
            'polynomial_moment_relative_max':maxmoment,'condition_number_max':maxcond,
            'residual_row_norm_min':minimum_residual,
            'work_closure_relative':float(abs(work+aux)/den),
            'net_force_relative':float(np.linalg.norm(force.sum(axis=0))/fnorm),
            'net_torque_relative':float(abs(torque)/(fnorm*.7))}
    assert maxmoment<1e-12 and result['work_closure_relative']<1e-12
    assert result['net_force_relative']<1e-12 and result['net_torque_relative']<1e-12
    results.append(result)
print(json.dumps({'scope':__doc__,'rows':results,
                  'state_identity':'persistent particle ID and Cartesian component; identical z in both branches',
                  'limitation':'Support swap preserves E_aux algebraically but force continuity, accumulated phase, shock thermodynamics and actual moving-state stability are unverified.'},indent=2))
