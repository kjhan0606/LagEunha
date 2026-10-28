from pathlib import Path
import numpy as np,json,sys,re
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_hg_apply_20261017")
fields=['x','y','mass','vx','vy','ie','volume','density','pressure','csound']

def read(label,n):
 a=np.concatenate([np.loadtxt(root/label/f'raw_faces.dat.hgafter.rank{r:04d}') for r in range(n)])
 assert a.shape==(36864,13) and np.isfinite(a).all();a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(len(a))) and np.all(a[:,2]==1)
 return a
base=read('np1',1);before=np.loadtxt(root/'np1/raw_faces.dat.hgstate.rank0000');before=before[np.argsort(before[:,0])];assert np.array_equal(before[:,0],base[:,0]);frozen=np.hypot(before[:,3]-3,before[:,4]-3)>=2.5
out={'scope':'One HG operator tau=0.0005 at frozen t10 geometry; original hg_apply vs ordered distributed owner updates; not RK evolution','cases':{}}
for label,n in [('np1',1),('ordered1',1),('np2',2),('np4',4)]:
 assert (root/label/'exit.code').read_text().strip()=='0';a=read(label,n)
 assert np.array_equal(a[:,[3,4,5,9,10]],before[:,[3,4,5,9,10]])
 assert np.array_equal(a[frozen,6:9],before[frozen,6:9]),'frozen velocity or energy changed'
 assert np.all(a[:,8]>0) and np.all(a[:,11:]>0)
 diff=np.abs(a[:,3:]-base[:,3:]);bits=np.ascontiguousarray(a[:,3:]).view(np.uint64)!=np.ascontiguousarray(base[:,3:]).view(np.uint64)
 log=(root/label/'log').read_text();rows=re.findall(r'\[MPI_HG_APPLY\] rank=(\d+) tau=([^ ]+) global_modes=(\d+) blocks=(\d+) remote_state_returns=(\d+) bad=(\d+) seconds=([^\n]+)',log);assert len(rows)==n and all(int(r[5])==0 for r in rows) and '[ACCEPTED_STATE]' not in log
 mass=a[:,5];dk=.5*mass*(np.sum(a[:,6:8]**2,axis=1)-np.sum(before[:,6:8]**2,axis=1));du=a[:,8]-before[:,8];dp=mass[:,None]*(a[:,6:8]-before[:,6:8]);dl=(a[:,3]-3)*dp[:,1]-(a[:,4]-3)*dp[:,0]
 out['cases'][label]={'ranks':n,'different_particles':int(np.count_nonzero(np.any(diff!=0,axis=1))),'bitwise_different_particles':int(np.count_nonzero(np.any(bits,axis=1))),'max_abs_difference_by_field':dict(zip(fields,np.max(diff,axis=0).tolist())),'mass_position_volume_density_unchanged':True,'frozen_velocity_ie_unchanged':True,'delta_K':float(dk.sum()),'delta_U':float(du.sum()),'delta_E_raw':float((dk+du).sum()),'delta_Lz_raw':float(dl.sum()),'schedule':[{'rank':int(r[0]),'tau':float(r[1]),'global_modes':int(r[2]),'blocks':int(r[3]),'remote_returns':int(r[4]),'bad':int(r[5]),'seconds':float(r[6])} for r in rows]}
(root/'apply_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
