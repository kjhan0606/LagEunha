"""Static auxiliary-energy operator audit on saved C HG modes, inner region.

Mode coefficients are rescaled by a prescribed acoustic frequency, not the
HG dissipative coefficient. No actual trajectory or snapshot is modified.
"""
import json, struct, hashlib
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive')
state_path=root/'raw_faces.dat.hgstate.rank0000'
mode_path=root/'raw_faces.dat.hgm.rank0000'
state=np.loadtxt(state_path);state=state[np.argsort(state[:,0])]
n=len(state);assert np.array_equal(state[:,0],np.arange(n))
x=state[:,3:5]-3;mass=state[:,5];velocity=state[:,6:8];U=state[:,8]
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
volume=rhs[:,6];gamma=float(np.float32(1.666667))
with mode_path.open('rb') as f: magic,nm,ns,ncol=struct.unpack('<4q',f.read(32))
assert magic==20261017 and ncol==9
assert mode_path.stat().st_size==32+nm*72+ns*24
modes=np.memmap(mode_path,mode='r',dtype='<f8',offset=32,shape=(nm,9))
support=np.memmap(mode_path,mode='r',dtype='<f8',offset=32+nm*72,shape=(ns,3))
rng=np.random.default_rng(20261018);probe_velocity=rng.normal(size=(n,2))
force=np.zeros((n,2));aux_work=0.;actual_aux_work=0.;norm_force=0.;max_mom=0.
selected=[];owner_counts={};support_rmax=0.;freq2_sum=0.;probe_abs_work=0.
for row in modes:
    owner,key,start,length,kind,onesided=map(int,row[:6])
    if np.linalg.norm(x[owner])>=.6:continue
    assert kind==0 and onesided==0
    ss=support[start:start+length];ids=ss[:,0].astype(int);L=np.array(ss[:,1:3])
    assert np.all(np.linalg.norm(x[ids],axis=1)<2.5)
    support_rmax=max(support_rmax,float(np.linalg.norm(x[ids],axis=1).max()))
    h=np.sqrt(volume[owner]);alpha=np.sum(L**2/mass[ids,None])
    assert alpha>0 and U[owner]>0
    omega=np.sqrt(gamma*(gamma-1)*U[owner]/mass[owner])/h
    B=L*omega/np.sqrt(alpha)
    assert abs(np.sum(B**2/mass[ids,None])/omega**2-1)<1e-12
    # Relative local polynomial moment error of original vector row.
    dx=(x[ids]-x[owner])/h
    poly=np.column_stack((np.ones(len(ids)),dx,dx[:,0]**2,dx[:,0]*dx[:,1],dx[:,1]**2))
    moments=poly.T@B
    normalization=np.abs(poly).T@np.abs(B)
    max_mom=max(max_mom,float(np.max(np.abs(moments)/np.maximum(normalization,1e-300))))
    z=float(rng.normal());f=-z*B
    np.add.at(force,ids,f)
    power=z*np.sum(B*probe_velocity[ids]);aux_work+=power;probe_abs_work+=abs(power)
    actual_aux_work+=z*np.sum(B*velocity[ids]);norm_force+=np.linalg.norm(f,axis=1).sum()
    selected.append(key);owner_counts[owner]=owner_counts.get(owner,0)+1
    freq2_sum+=omega**2
kinetic_work=float(np.sum(force*probe_velocity));actual_kinetic=float(np.sum(force*velocity))
torque=float(np.sum(x[:,0]*force[:,1]-x[:,1]*force[:,0]))
force_error=float(np.linalg.norm(force.sum(axis=0))/norm_force)
torque_error=abs(torque)/(norm_force*support_rmax)
work_error=abs(kinetic_work+aux_work)/probe_abs_work
assert max_mom<1e-9 and force_error<1e-12 and torque_error<1e-12 and work_error<1e-12
def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()
out={'scope':__doc__,'mode_sha256':sha(mode_path),'state_sha256':sha(state_path),
     'source_time':6.000405917763673,'owner_radius_max':.6,'support_radius_max':support_rmax,
     'modes':len(selected),'owners':len(owner_counts),
     'modes_per_owner_range':[min(owner_counts.values()),max(owner_counts.values())],
     'polynomial_degree_2_relative_moment_max':max_mom,
     'net_force_relative':force_error,'net_torque_relative':torque_error,
     'random_velocity_work_relative_residual':work_error,
     'actual_velocity_work_absolute_residual':float(actual_kinetic+actual_aux_work),
     'frequency_squared_trace_bound':freq2_sum,
     'key_equals_owner_times_256_plus_local_face_slot':bool(all(k//256 in owner_counts for k in selected)),
     'decision':'Static operator contract passes only within audited interior support. Existing face-slot keys do not define persistent auxiliary states. No evolution adoption.'}
print(json.dumps(out,indent=2))
