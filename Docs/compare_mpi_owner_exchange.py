from pathlib import Path
import json,sys,collections
import numpy as np
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_owner_exchange_20261017")
reference=Path('/gpfs/kjhan/LagForce/Codex_mpi_cache_probe_20261017/np1')
a=np.loadtxt(reference/'raw_faces.dat.cache.rank0000');cache={int(r[0]):r for r in a};faces=collections.defaultdict(list)
with (reference/'raw_faces.dat.rank0000').open() as f:
 for line in f:
  if line.startswith('#') or not line.strip():continue
  r=line.split();faces[int(r[0])].append([int(r[1])]+list(map(float,r[4:10])))
out={'scope':'Owner-directed runtime MPI exchange after distributed tessellation; imported packets checked against independent one-rank geometry/cache output; force reads not enabled','ranks':{}}
for n in [1,2,4]:
 run=root/f'np{n}';assert (run/'exit.code').read_text().strip()=='0'
 owner={};owned=[]
 for rank in range(n):
  rows=np.loadtxt(run/f'raw_faces.dat.cache.rank{rank:04d}');ids=set(map(int,rows[:,0]));assert not ids&owner.keys();owner.update({i:rank for i in ids});owned.append(ids)
 assert set(owner)==set(cache)
 reports=[]
 for rank in range(n):
  expected={int(f[0]) for i in owned[rank] for f in faces[i] if owner[int(f[0])]!=rank}
  received={};rf=collections.defaultdict(list)
  with (run/f'raw_faces.dat.exchange.rank{rank:04d}').open() as file:
   for line in file:
    if line.startswith('#'):continue
    s=line.split();i=int(s[1])
    if s[0]=='C':assert i not in received;received[i]=list(map(float,s[2:]))
    else:assert s[0]=='F';rf[i].append([int(s[2])]+list(map(float,s[3:])))
  assert set(received)==expected
  max_cache=0.;max_geom=0.;totalfaces=0
  for i,row in received.items():
   c=cache[i];assert row[:5]==[owner[i],1,len(faces[i]),int(c[5]),int(c[6])]
   target=np.r_[c[3:5],c[7:]];max_cache=max(max_cache,float(np.max(np.abs(np.array(row[5:])-target))))
   b=np.array(sorted(rf[i]));ref=np.array(sorted(faces[i]));assert b.shape==ref.shape
   max_geom=max(max_geom,float(np.max(np.abs(b-ref))));totalfaces+=len(b)
  assert max_cache==max_geom==0
  reports.append({'rank':rank,'requested_and_received_cells':len(received),'received_faces':totalfaces,'max_cache_abs_difference':max_cache,'max_geometry_abs_difference':max_geom})
 out['ranks'][str(n)]=reports
stale=root/'stale2';code=int((stale/'exit.code').read_text());log=(stale/'log').read_text();assert code!=0 and '[MPI_GEOMETRY_EXCHANGE] invalid_records=1' in log
out['stale_generation_control']={'ranks':2,'exit_code':code,'invalid_records':1,'rejected':True}
(root/'exchange_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
