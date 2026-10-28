from pathlib import Path
import numpy as np,json,re,sys
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_hg_rebuild_20261017")
fields=['x','y','mass','vx','vy','ie','volume','density','pressure','csound']
out={'scope':'HG-only repeated operator with reconstruction on fixed t12 positions; no RK, gravity, particle migration or retry','cases':{}}

def load(label,n,stem):
 a=np.concatenate([np.loadtxt(root/label/f'{stem}.rank{r:04d}') for r in range(n)])
 assert np.isfinite(a).all();return a[np.argsort(a[:,0])]
for case,repeats in [('short',1),('repeat',2)]:
 for n in [1,4]:
  assert (root/f'{case}_np{n}/exit.code').read_text().strip()=='0'
  log=(root/f'{case}_np{n}/log').read_text();assert len(re.findall(r'\[MPI_PROBE_REBUILD\]',log))==n*repeats
  assert len(re.findall(r'\[MPI_HG_APPLY\]',log))==n*repeats and '[ACCEPTED_STATE]' not in log
 reports=[];previous=None;previous_cache=None
 for passno in range(1,repeats+1):
  prefix='raw_faces.dat'+(f'.repeat{passno:02d}' if repeats>1 else '')
  a=load(f'{case}_np1',1,prefix+'.hgafter');b=load(f'{case}_np4',4,prefix+'.hgafter');initial=load(f'{case}_np1',1,prefix+'.hgstate')
  assert a.shape==b.shape==initial.shape==(36864,13)
  assert np.array_equal(a[:,0],np.arange(36864)) and np.array_equal(b[:,0],a[:,0])
  assert np.all(a[:,2]==passno) and np.all(b[:,2]==passno) and np.all(initial[:,2]==passno)
  assert np.array_equal(a[:,[3,4,5,9,10]],initial[:,[3,4,5,9,10]])
  frozen=np.hypot(initial[:,3]-3,initial[:,4]-3)>=2.5
  assert np.array_equal(a[frozen,6:9],initial[frozen,6:9]) and np.array_equal(b[frozen,6:9],initial[frozen,6:9])
  diff=np.abs(a[:,3:]-b[:,3:]);bits=np.ascontiguousarray(a[:,3:]).view(np.uint64)!=np.ascontiguousarray(b[:,3:]).view(np.uint64)
  ca=load(f'{case}_np1',1,prefix+'.cache');cb=load(f'{case}_np4',4,prefix+'.cache');assert ca.shape==cb.shape==(36864,25) and np.array_equal(ca[:,0],cb[:,0])
  cd=np.max(np.abs(ca[:,3:]-cb[:,3:]));assert cd==0
  report={'pass':passno,'max_abs_state_difference':dict(zip(fields,np.max(diff,axis=0).tolist())),'bitwise_different_particles':int(np.count_nonzero(np.any(bits,axis=1))),'reconstruction_max_abs_difference':float(cd),'epoch':passno,'frozen_velocity_ie_unchanged':True}
  if previous is not None:
   assert np.array_equal(initial[:,6:9],previous[:,6:9]),'rebuild did not preserve prior applied velocity/internal energy'
   changed=int(np.count_nonzero(np.any(ca[:,7:]!=previous_cache[:,7:],axis=1)));assert changed>0
   report.update(previous_applied_state_preserved=True,changed_reconstruction_cells_from_previous_pass=changed)
  previous=a;previous_cache=ca;reports.append(report)
 out['cases'][case]=reports
(root/'rebuild_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
