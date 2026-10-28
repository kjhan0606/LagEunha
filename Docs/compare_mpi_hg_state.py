from pathlib import Path
import numpy as np,json,sys,collections
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_hg_state_20261017")
base=np.loadtxt(root/'np1/raw_faces.dat.hgstate.rank0000');assert base.shape==(36864,13) and np.isfinite(base).all();base=base[np.argsort(base[:,0])];assert np.array_equal(base[:,0],np.arange(36864))
adj=collections.defaultdict(set);original={}
with Path('/gpfs/kjhan/LagForce/Codex_mpi_geometry_probe_20261017/np1/raw_faces.dat.rank0000').open() as file:
 for line in file:
  if line.startswith('# cell '):
   a=line.split()[1:];d=dict(zip(a[::2],a[1::2]));original[int(d['cell'])]=d
  elif not line.startswith('#') and line.strip():
   a=line.split();adj[int(a[0])].add(int(a[1]))
fields=['x','y','mass','vx','vy','ie','volume','density','pressure','csound']
for i,row in enumerate(base):
 for k,key in enumerate(fields[:-1]):assert row[3+k]==float(original[i][key]),(i,key)
out={'scope':'Actual owner-request exchange and installed hg_ptr state for all raw two-hop neighbors; modes not constructed or applied','ranks':{}}
for n in [1,2,4]:
 assert (root/f'np{n}/exit.code').read_text().strip()=='0'
 owner={};owned=[]
 for rank in range(n):
  a=np.loadtxt(root/f'np{n}/raw_faces.dat.cache.rank{rank:04d}');ids=set(map(int,a[:,0]));assert not ids&owner.keys();owner.update({i:rank for i in ids});owned.append(ids)
 assert len(owner)==36864
 reports=[]
 for rank in range(n):
  direct=set().union(*(adj[i] for i in owned[rank]));second=set().union(*(adj[i] for i in direct));expected=owned[rank]|direct|second
  a=np.loadtxt(root/f'np{n}/raw_faces.dat.hgstate.rank{rank:04d}');assert a.shape[1]==13 and np.isfinite(a).all();ids=a[:,0].astype(int);assert len(set(ids))==len(ids) and set(ids)==expected
  assert np.all(a[:,2]==1) and all(int(row[1])==owner[int(row[0])] for row in a)
  diff=np.max(np.abs(a[:,3:]-base[ids,3:]),axis=0);assert np.max(diff)==0
  reports.append({'rank':rank,'owned':len(owned[rank]),'direct_remote':len(direct-owned[rank]),'all_remote_state':len(expected-owned[rank]),'additional_second_hop_remote':len(expected-owned[rank]-direct),'max_abs_difference_by_field':dict(zip(fields,diff.tolist()))})
 out['ranks'][str(n)]=reports
stale=root/'stale2';code=int((stale/'exit.code').read_text());log=(stale/'log').read_text();assert code!=0 and '[MPI_HG_STATE] invalid_records=1' in log
out['stale_state_negative_control']={'ranks':2,'exit_code':code,'invalid_records':1,'rejected':True}
(root/'state_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
