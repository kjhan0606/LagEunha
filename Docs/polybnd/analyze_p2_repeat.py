import json
from pathlib import Path
import numpy as np
ref=np.load('reference/reference.npz');records=[]
for mode in ('numpy','omp1','omp4'):
    for rep in (1,2):
        folder=Path(f'{mode}_{rep}');backend='numpy' if mode=='numpy' else 'openmp'
        z=np.load(folder/(backend+'.npz'));timing=json.loads((folder/(backend+'.timing.json')).read_text())
        exact=all(np.array_equal(ref[k],z[k]) for k in ref.files)
        close=all(ref[k].shape==z[k].shape and (np.array_equal(ref[k],z[k]) if k=='keys' else np.allclose(ref[k],z[k],rtol=1e-10,atol=1e-12)) for k in ref.files)
        thread_ok=mode=='numpy' or timing['omp_actual']==int(timing['omp_requested'])
        records.append(dict(mode=mode,rep=rep,seconds=timing['seconds'],exact_equal=exact,equivalence_pass=bool(close),thread_check=thread_ok))
medians={mode:float(np.median([r['seconds'] for r in records if r['mode']==mode])) for mode in ('numpy','omp1','omp4')}
base=json.loads(Path('reference/reference.timing.json').read_text())['seconds']
valid=all(r['equivalence_pass'] and r['thread_check'] for r in records)
fastest=min(medians,key=medians.get)
# Prefer the dependency-free single-core option when timings differ by <=5%.
choice='numpy' if medians['numpy']<=1.05*medians[fastest] else fastest
report=dict(records=records,median_seconds=medians,reference_seconds=base,speedup={k:base/v for k,v in medians.items()},all_pass=valid,recommended_backend=choice if valid else None,scope='Two repeats, n64 T=1 only; long-run equivalence untested.')
Path('repeat_summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if not valid:raise RuntimeError('repeat comparison failed')
