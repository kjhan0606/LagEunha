from pathlib import Path
import numpy as np,json,hashlib
root=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017');n=36864
for case in ['control','operator']:assert (root/case/'exit.code').read_text().strip()=='0'
def read(case,filename):
 a=np.loadtxt(root/case/filename);return a[np.argsort(a[:,0])]
a=read('operator','rhs.dat');control=read('control','rhs.dat');assert np.array_equal(a,control),'export modified RHS'
cache=read('operator','raw_faces.dat.cache.rank0000');assert np.all(cache[:,5]==2),'incomplete QR coverage'
raw=np.fromfile(root/'operator/gradient_operator.bin',dtype='<f8');assert raw.size%6==0;rows=raw.reshape(-1,6);assert np.isfinite(rows).all()
i=rows[:,0].astype(np.int64);j=rows[:,1].astype(np.int64);assert np.all(rows[:,0]==i)&np.all(rows[:,1]==j);assert i.min()==j.min()==0 and i.max()==j.max()==n-1;assert np.array_equal(np.unique(i),np.arange(n))
wx,wy,dx,dy=rows[:,2:].T
sumrow=lambda q:np.bincount(i,weights=q,minlength=n)
trace=sumrow(wx*a[j,3]+wy*a[j,4]);stored=cache[:,7]+cache[:,10]
error=trace-stored;assert np.max(np.abs(error))<1e-10
moment={'constant_x':sumrow(wx),'constant_y':sumrow(wy),'d_x_of_x_minus_1':sumrow(wx*dx)-1,'d_x_of_y':sumrow(wx*dy),'d_y_of_x':sumrow(wy*dx),'d_y_of_y_minus_1':sumrow(wy*dy)-1,'rigid_rotation_divergence':sumrow(-wx*dy+wy*dx)}
for key,val in moment.items():assert np.max(np.abs(val))<1e-10,(key,np.max(np.abs(val)))
vol=a[:,6];bx=vol[i]*wx;by=vol[i]*wy
transpose=lambda p:np.column_stack((np.bincount(j,weights=bx*p[i],minlength=n),np.bincount(j,weights=by*p[i],minlength=n)))
p0=2.7777796321447566e-5;f0=transpose(np.full(n,p0));old0=-p0*a[:,15:17]
rng=np.random.default_rng(20261017);v=rng.normal(size=(n,2));p=rng.uniform(.5,1.5,n);f=transpose(p);vdot=sumrow(bx*v[j,0]+by*v[j,1]);forcework=float(np.sum(v*f));internalwork=-float(np.dot(p,vdot));scale=float(np.sum(np.abs(v*f))+np.sum(np.abs(p*vdot)));relative=abs(forcework+internalwork)/scale;assert relative<1e-12
r=np.hypot(a[:,1]-3,a[:,2]-3);regions={}
for name,sel in [('inner',r<.5),('disk_inner_edge',(r>=.5)&(r<.6)),('disk_core',(r>=.6)&(r<1.9))]:
 acc=np.linalg.norm(f0[sel],axis=1)/a[sel,7];regions[name]={'n':int(sel.sum()),'candidate_constant_pressure_force_L1':float(np.linalg.norm(f0[sel],axis=1).sum()),'original_constant_pressure_force_L1':float(np.linalg.norm(old0[sel],axis=1).sum()),'candidate_constant_pressure_acceleration_mass_weighted_RMS':float(np.sqrt(np.average(acc**2,weights=a[sel,7]))),'candidate_constant_pressure_acceleration_max':float(acc.max()),'divergence_reproduction_max_abs_error':float(np.max(np.abs(error[sel])))}
out={'scope':'Actual C QR coefficients exported at canonical t2 circular diagnostic state. B_ij=V_i times fitted divergence coefficients; test B transpose pressure without changing simulation. Linear-map checks do not prove a viable volume model.','records':len(rows),'all_C_RHS_unchanged_with_export':True,'gradient_trace_reproduction_max_abs_error':float(np.max(np.abs(error))),'local_moment_max_abs_error':{k:float(np.max(np.abs(v))) for k,v in moment.items()},'random_state_adjoint_energy_identity_relative_residual':relative,'constant_pressure_regions':regions,'global_B_transpose_one_L1':float(np.linalg.norm(f0/p0,axis=1).sum()),'binary_sha256':(root/'binary.sha256').read_text().strip(),'operator_sha256':hashlib.sha256(raw.tobytes()).hexdigest(),'decision':'Assess constant-pressure closure before any evolution. Adjoint energy matching alone does not ensure constant-pressure equilibrium.'}
(root/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
