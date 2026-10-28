from pathlib import Path
import json,sys,re
import numpy as np
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_force_probe_20261017")

def read(n,face=False):
 stem='raw_faces.dat.rhs.faces.rank' if face else 'raw_faces.dat.rhs.rank'
 a=np.concatenate([np.loadtxt(root/f'np{n}'/f'{stem}{r:04d}') for r in range(n)])
 assert np.isfinite(a).all()
 if face:
  assert a.shape[1]==21
  return a[np.lexsort((a[:,7],a[:,6],a[:,1],a[:,0]))]
 assert a.shape==(36864,25) and len(np.unique(a[:,0]))==len(a)
 return a[np.argsort(a[:,0])]

def comparison(a,b,face=False):
 assert a.shape==b.shape and np.array_equal(a[:,:2] if face else a[:,0],b[:,:2] if face else b[:,0])
 diff=np.abs(a-b);scale=np.max(np.abs(a),axis=0)
 return {'rows':len(a),'different_rows':int(np.count_nonzero(np.any(diff!=0,axis=1))),'max_abs_difference_by_column':np.max(diff,axis=0).tolist(),'max_difference_divided_by_column_max':(np.max(diff,axis=0)/np.maximum(scale,1e-300)).tolist()}
base,basef=read(1),read(1,True);out={'scope':'Single pre-HG hydrodynamic RHS and face dissipation at t10; gravity, HG and RK evolution not applied','ranks':{}}
for n in [1,2,4]:
 run=root/f'np{n}';assert (run/'exit.code').read_text().strip()=='0';log=(run/'log').read_text();dt=re.findall(r'\[MPI_FORCE_ONLY\] rank=\d+ dt=([^ ]+) complete',log);assert len(dt)==n and '[ACCEPTED_STATE]' not in log
 out['ranks'][str(n)]={'cells':comparison(base,read(n)),'faces':comparison(basef,read(n,True),True),'globally_reduced_dt':list(map(float,dt))}
old=Path('/gpfs/kjhan/LagForce/Codex_TCF_rhs_t10_20261017/t10');a=np.loadtxt(old/'rhs.dat');a=a[np.argsort(a[:,0])];b=np.loadtxt(old/'rhs.dat.faces');b=b[np.lexsort((b[:,7],b[:,6],b[:,1],b[:,0]))]
out['original_single_rank_reference']={'cells':comparison(a,base),'faces':comparison(b,basef,True)}
(root/'rhs_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
