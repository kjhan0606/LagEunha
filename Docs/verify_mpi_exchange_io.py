from pathlib import Path
import json,hashlib,subprocess
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_mpi_exchange_io_20261017')
old=Path('/gpfs/kjhan/LagForce/Codex_mpi_rollback_20261017')
subprocess.run(['python3',str(root/'compare.py')],check=True)

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

out={'scope':'Exchange I/O separation; real migration and rollback diagnostic, not coupled MPI RK evolution.','cases':{}}
for name,n in [('control',1),('np1',1),('np2',2),('np4',4)]:
 d=root/name;files=list(d.iterdir())
 assert not any(any(part in p.name for part in ['.exchange.rank','.common.rank','.hgstate.rank']) for p in files)
 count=0
 for phase in [1,2]:
  for rank in range(n):
   for kind in ['cache','hgafter','hgm','rhs','rhs.faces']:
    fn=f'raw_faces.dat.repeat{phase:02d}.{kind}.rank{rank:04d}'
    if kind in ['rhs','rhs.faces']:
     a=np.loadtxt(d/fn);b=np.loadtxt(old/name/fn)
     assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
     keys=[0] if kind=='rhs' else [0,1]
     assert len(np.unique(a[:,keys],axis=0))==len(a)
     assert len(np.unique(b[:,keys],axis=0))==len(b)
     a=a[np.lexsort(tuple(a[:,k] for k in reversed(keys)))];b=b[np.lexsort(tuple(b[:,k] for k in reversed(keys)))]
     assert np.array_equal(a.view(np.uint64),b.view(np.uint64)),(name,fn)
    else:
     assert digest(d/fn)==digest(old/name/fn),(name,fn)
    count+=1
 out['cases'][name]={'exchange_dump_files':0,'reference_records_identical':count,'comparison':'RHS records keyed by cell or face IDs; cache, HG state and mode files byte-identical'}
for case,expected,marker in [('stale_geometry',3,'[MPI_GEOMETRY_EXCHANGE] invalid_records='),('stale_state',4,'[MPI_HG_STATE] invalid_records=')]:
 code=int((root/case/'exit.code').read_text());log=(root/case/'log').read_text()
 assert code==expected and marker in log,(case,code)
 out[case]={'exit_code':code,'expected_invalid_record_rejection':True}
(root/'io_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
