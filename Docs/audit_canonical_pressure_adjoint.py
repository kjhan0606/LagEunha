from pathlib import Path
import numpy as np,json,hashlib
root=Path(__file__).parent
source=Path('/gpfs/kjhan/LagForce/Codex_t2_onset_rhs_20261017')
out={'scope':'Static original C RHS at t2; full precision exported real-face geometry. Split HG and finite-step RK errors excluded. Global virtual work includes frozen cells as diagnostic degrees of freedom, not physical boundary work.'}
for case in ['t2','circular']:
 c=np.loadtxt(source/case/'rhs.dat');c=c[np.argsort(c[:,0])];n=len(c)
 f=np.loadtxt(source/case/'raw_faces.dat.rank0000');i,j=f[:,:2].astype(np.int64).T
 assert np.array_equal(c[:,0],np.arange(n))
 dr=f[:,4:6];d=np.linalg.norm(dr,axis=1);e=dr/d[:,None];A=np.linalg.norm(f[:,8:10]-f[:,6:8],axis=1)
 t=.5*(f[:,6:8]+f[:,8:10])-dr*.5;t-=np.sum(t*e,axis=1)[:,None]*e
 def scatter(z,idx=i):
  a=np.zeros((n,)+z.shape[1:]);np.add.at(a,idx,z);return a
 def B(v):return scatter(A*(.5*np.sum((v[i]+v[j])*e,axis=1)+np.sum((v[i]-v[j])*t,axis=1)/d))
 def BT(p):return scatter(p[i,None]*A[:,None]*(.5*e+t/d[:,None]))+scatter(p[i,None]*A[:,None]*(.5*e-t/d[:,None]),j)
 P=c[:,5];v=c[:,3:5]
 normal=scatter(-A[:,None]*.5*(P[i]+P[j])[:,None]*e)
 tangent=scatter(-A[:,None]*(P[j]-P[i])[:,None]*t/d[:,None]);F=normal+tangent
 closure=scatter(A[:,None]*e);vd=B(v);bt=BT(P)
 q=np.loadtxt(source/case/'rhs.dat.faces');qi=q[:,0].astype(int)
 Fq=scatter(-q[:,4,None]*q[:,14,None]*q[:,6:8],qi)
 r=np.linalg.norm(c[:,1:3]-3,axis=1)
 res={ 'rhs_sha256':hashlib.sha256((source/case/'rhs.dat').read_bytes()).hexdigest(),'raw_geometry_sha256':hashlib.sha256((source/case/'raw_faces.dat.rank0000').read_bytes()).hexdigest(),'regions':{} }
 for name,sel in [('inner',r<.5),('disk_edge',(r>=.5)&(r<.6)),('disk_core',(r>=.6)&(r<1.9)),('active',r<2.5),('frozen',r>=2.5),('all',np.ones(n,dtype=bool))]:
  norm=lambda z:float(np.abs(z[sel]).sum())
  res['regions'][name]={'n':int(sel.sum()),'geometric_Vdot_relative_L1_error':norm(vd-c[:,12])/max(norm(c[:,12]),1e-300),'pressure_force_plus_Q_vs_C_relative_L1_error':norm(F+Fq-c[:,8:10])/max(norm(c[:,8:10]),1e-300),'adjoint_force_relative_L1_difference':norm(bt-F)/max(norm(F),1e-300),'constant_pressure_closure_L1':norm(closure),'internal_energy_identity_relative_L1_residual':norm(c[:,10]-np.sum(v*c[:,8:10],axis=1)+P*c[:,12]-c[:,11])/max(norm(P*c[:,12])+norm(c[:,11]),1e-300)}
 res['global_virtual_work']={'v_dot_pressure_force':float(np.sum(v*F)),'P_dot_geometric_Vdot':float(P@vd),'difference':float(np.sum(v*F)-P@vd)}
 if case=='t2':
  circ=np.loadtxt(source/'circular/rhs.dat');circ=circ[np.argsort(circ[:,0])];assert np.array_equal(c[:,[0,1,2,5,6,7]],circ[:,[0,1,2,5,6,7]])
  dv=v-circ[:,3:5]
  for name,sel in [('inner',r<.5),('disk_edge',(r>=.5)&(r<.6)),('disk_core',(r>=.6)&(r<1.9))]:
   res['regions'][name]['perturbation_pressure_work_normal']=float(np.sum(dv[sel]*normal[sel]))
   res['regions'][name]['perturbation_pressure_work_tangent']=float(np.sum(dv[sel]*tangent[sel]))
  res['perturbation_virtual_work']={'force_work':float(np.sum(dv*F)),'volume_work':float(P@B(dv)),'difference':float(np.sum(dv*F)-P@B(dv))}
 bad=np.linalg.norm(bt-F,axis=1)>1e-12
 res['adjoint_mismatch_locations']={'threshold_absolute_force':1e-12,'count':int(bad.sum()),'minimum_radius':float(r[bad].min()) if bad.any() else None,'maximum_radius':float(r[bad].max()) if bad.any() else None}
 out[case]=res
(root/'audit_v2.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
