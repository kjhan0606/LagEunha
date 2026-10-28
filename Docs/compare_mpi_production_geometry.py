from pathlib import Path
import json,sys,numpy as np
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_geometry_probe_20261017")

def read(paths):
 cells={};faces={}
 for path in paths:
  with path.open() as f:
   for line in f:
    if line.startswith('# cell '):
     a=line.split()[1:];d=dict(zip(a[::2],a[1::2]));i=int(d.pop('cell'));rank=int(d.pop('rank'))
     assert i not in cells,('duplicate owner',i)
     cells[i]={k:float(v) for k,v in d.items()}
    elif not line.startswith('#') and line.strip():
     a=line.split();i,j=int(a[0]),int(a[1]);v=np.array([float(x) for x in a[2:10]])
     key=(i,j,round(v[2],8),round(v[3],8))
     assert key not in faces,('duplicate face/image',key)
     endpoints=sorted([tuple(v[4:6]),tuple(v[6:8])]);v[4:]=np.array(endpoints).ravel();faces[key]=v
 return cells,faces

basec,basef=read([root/'np1/raw_faces.dat.rank0000'])
assert len(basec)==36864
out={'scope':'Actual restart migration and per-rank cell construction; no force, HG or time evolution','ranks':{}}
for n in [1,2,4]:
 paths=[root/f'np{n}/raw_faces.dat.rank{r:04d}' for r in range(n)]
 assert (root/f'np{n}/exit.code').read_text().strip()=='0'
 log=(root/f'np{n}/log').read_text();assert log.count('complete; no evolution')==n
 assert '[ACCEPTED_STATE]' not in log
 c,f=read(paths);common=basef.keys()&f.keys()
 out['ranks'][str(n)]={'owned_cells':len(c),'missing_ids':len(basec.keys()-c.keys()),'extra_ids':len(c.keys()-basec.keys()),'raw_faces':len(f),'missing_face_keys':len(basef.keys()-f.keys()),'extra_face_keys':len(f.keys()-basef.keys()),'max_face_component_abs_difference':np.max(np.array([np.abs(f[k]-basef[k]) for k in common]),axis=0).tolist(),'cell_max_abs_difference':{field:max(abs(c[i][field]-basec[i][field]) for i in c.keys()&basec.keys()) for field in next(iter(basec.values()))}}
_,old=read([Path('/gpfs/kjhan/LagForce/Codex_shf_raw_geometry_20261017/audit/raw_faces.dat')]);common=old.keys()&basef.keys()
out['single_rank_vs_original']={'missing':len(old.keys()-basef.keys()),'extra':len(basef.keys()-old.keys()),'max_abs':float(max(np.max(np.abs(old[k]-basef[k])) for k in common))}
(root/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
