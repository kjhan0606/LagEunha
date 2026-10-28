"""Compare completed TCF dumps after canonical ordering; output precision only."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser();ap.add_argument('root');ap.add_argument('--output',required=True)
args=ap.parse_args(); root=Path(args.root)
dirs=[root/'omp1',root/'omp8']
for d,n in zip(dirs,[1,8]):
    assert (d/'exit.code').read_text().strip()=='0'
    log=(d/'log').read_text()
    assert '[TCF] dump written' in log
    assert f"OMP_NUM_THREADS='{n}'" in log or f"OMP_NUM_THREADS='{n}'" in log.replace(' ','')
    assert "OMP_DYNAMIC='FALSE'" in log or "OMP_DYNAMIC='FALSE'" in log.replace(' ','')
assert (dirs[0]/'binary.sha256').read_text()==(dirs[1]/'binary.sha256').read_text()
result={'scope':'One TCF force evaluation at t6.002110493120864 on MPI1; CUDA nearest distance retained. Equality at printed precision, not binary memory equality or long evolution parity.','tables':{}}
for file,cols in [('rhs.dat',25),('rhs.dat.faces',21)]:
    arrays=[];hashes=[]
    for d in dirs:
        p=d/file; a=np.loadtxt(p)
        assert a.shape[1]==cols and np.isfinite(a).all()
        # Full lexicographic ordering handles repeated directed-face keys.
        a=a[np.lexsort(tuple(a[:,j] for j in reversed(range(cols))))]
        arrays.append(a);hashes.append(hashlib.sha256(p.read_bytes()).hexdigest())
    a,b=arrays; assert a.shape==b.shape
    different=a.view(np.uint64)!=b.view(np.uint64)
    result['tables'][file]={'rows':len(a),'columns':cols,'different_values':int(different.sum()),
                           'max_abs_difference':float(np.abs(a-b).max()),'raw_sha256':hashes}
assert all(t['different_values']==0 for t in result['tables'].values())
Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
