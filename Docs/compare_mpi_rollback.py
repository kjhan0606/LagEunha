from pathlib import Path
import numpy as np,json,re
root=Path(__file__).parent

def read(case,n,stem,cols):
 a=np.concatenate([np.atleast_2d(np.loadtxt(root/case/f'raw_faces.dat.{stem}.rank{r:04d}')) for r in range(n)])
 assert a.shape==(36864,cols) and np.isfinite(a).all()
 a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(36864))
 return a

def same(a,b):
 return np.array_equal(np.ascontiguousarray(a).view(np.uint64),np.ascontiguousarray(b).view(np.uint64))

assert (root/'control/exit.code').read_text().strip()=='0'
ref=read('control',1,'repeat02.hgafter',13)
ref_cache=read('control',1,'repeat02.cache',25)
out={'scope':'Diagnostic position drift, actual migration, last-rank rejection, collective owned-state restore, cache rebuild and second HG application. Not production RK retry; no trial ledger or clock mutation is exercised.','cases':{}}
for n in [1,2,4]:
 case=f'np{n}';assert (root/case/'exit.code').read_text().strip()=='0'
 before=read(case,n,'repeat01.motion_before',37);after=read(case,n,'repeat01.motion_after',37);restored=read(case,n,'repeat01.motion_restored',37)
 assert same(before,restored),'state, ownership or RK checkpoint mismatch'
 changed=int(np.count_nonzero(before[:,1]!=after[:,1]));assert n==1 or changed>0
 assert same(before[:,4:],after[:,4:]),'migration modified non-position fields'
 rows=re.findall(r'\[MPI_PROBE_ROLLBACK\] rank=(\d+) requester=(\d+) local_reject=(\d+) global_reject=(\d+) restored_np=(\d+) bytes_ok=(\d+)',(root/case/'log').read_text())
 assert len(rows)==n and sorted(int(x[0]) for x in rows)==list(range(n))
 assert sum(int(x[2]) for x in rows)==1
 assert all(int(x[1])==n-1 and int(x[3])==1 and int(x[5])==1 for x in rows)
 final=read(case,n,'repeat02.hgafter',13);cache=read(case,n,'repeat02.cache',25)
 assert np.all(final[:,2]==2) and np.all(cache[:,2]==2)
 assert same(final[:,3:],ref[:,3:]),'post-restore HG differs from no-drift control'
 assert same(cache[:,3:],ref_cache[:,3:]),'rebuilt reconstruction differs from control'
 out['cases'][case]={'ownership_crossings_before_reject':changed,'single_requester':n-1,'all_ranks_rejected':True,'owned_state_and_25_RK_scalars_restored_bitwise':True,'original_ownership_restored':True,'rebuilt_cache_bitwise_matches_no_drift_control':True,'second_HG_state_bitwise_matches_no_drift_control':True}
(root/'comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
