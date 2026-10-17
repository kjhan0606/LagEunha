import json
from pathlib import Path
import numpy as np
root=Path.cwd();ref=np.load(root/'reference/reference.npz');report={}
reference_time=json.loads((root/'reference/reference.timing.json').read_text())['seconds']
for folder,mode in [('reference','reference'),('numpy','numpy'),('omp1','openmp'),('omp4','openmp')]:
    z=np.load(root/folder/(mode+'.npz'));timing=json.loads((root/folder/(mode+'.timing.json')).read_text())
    errors={};passed=True
    for field in ref.files:
        errors[field]=float(np.max(abs(z[field]-ref[field]))) if z[field].shape==ref[field].shape else None
        passed &= z[field].shape==ref[field].shape and (np.array_equal(z[field],ref[field]) if field=='keys' else np.allclose(z[field],ref[field],rtol=1e-10,atol=1e-12))
    report[folder]=dict(timing=timing,reference_speedup=reference_time/timing['seconds'],equivalence_pass=bool(passed),max_abs_difference=errors)
report['recommended_measured_backend']=min((k for k in ('reference','numpy','omp1','omp4') if report[k]['equivalence_pass']),key=lambda k:report[k]['timing']['seconds'])
report['scope']='Single short run per backend; no long-term equivalence or scaling claim. Profile run excluded from timing ratio.'
Path('performance_summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if not all(report[k]['equivalence_pass'] for k in ('numpy','omp1','omp4')):raise RuntimeError('geometry optimization changed short-run results beyond tolerance')
