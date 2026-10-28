from pathlib import Path
import json,sys,collections
import numpy as np
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_cache_probe_20261017")
names=['cdx','cdy','qbp_ok','qbp_arc']+[f'g{k}' for k in range(4)]+[f'h{k}' for k in range(6)]+[f't{k}' for k in range(8)]

def read(n):
 rows=[];counts=[]
 for rank in range(n):
  a=np.loadtxt(root/f'np{n}/raw_faces.dat.cache.rank{rank:04d}')
  assert a.ndim==2 and a.shape[1]==25
  assert np.isfinite(a).all() and np.all(a[:,1]==rank) and np.all(a[:,2]==1)
  counts.append(len(a));rows.append(a)
 a=np.concatenate(rows);assert len(a)==36864 and len(np.unique(a[:,0]))==len(a)
 return a[np.argsort(a[:,0])],counts
base,_=read(1);out={'scope':'Owned centroid and QBP cache after actual MPI restart migration and geometry generation; no ghost cache exchange or evolution','ranks':{}}
for n in [1,2,4]:
 a,counts=read(n);assert np.array_equal(a[:,0],base[:,0])
 diff=np.abs(a[:,3:]-base[:,3:]);out['ranks'][str(n)]={'owned_per_rank':counts,'maximum_abs_difference':dict(zip(names,np.max(diff,axis=0).tolist())),'different_cells':int(np.count_nonzero(np.any(diff!=0,axis=1))),'qbp_status_counts':dict(collections.Counter(map(int,a[:,5]))),'arc_counts':dict(collections.Counter(map(int,a[:,6])))}
(root/'cache_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
