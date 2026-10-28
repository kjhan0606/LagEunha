"""One-shot comparison of complete snapshots shared by both running branches.

Never treats absent output or a partially written snapshot as run termination.
This file reports partial observations only; final gates remain in the runner.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path

import numpy as np

ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args()
runner=Path('/gpfs/kjhan/LagForce/run_canonical_av_branches_20261018.py')
spec=importlib.util.spec_from_file_location('branches',runner);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
manifest=json.loads((args.root/'manifest.json').read_text())
source=Path(manifest['source']);expected_size=source.stat().st_size
initial,t0=module.read(source)
ic,_=module.read(module.BASE/'Codex_canonical_initial_20261017/run/khout.000000.dat')
r0=np.hypot(ic['x']-3,ic['y']-3);inner=(r0<.5)&(ic['den']<.5)
names=[]
for case in manifest['cases']:
    d=args.root/case['name']
    names.append({p.name for p in d.glob('khout.*.dat') if p.stat().st_size==expected_size})
common=sorted(set.intersection(*names))
rows=[]
for name in common:
    row={'snapshot':name,'cases':{}};times=[]
    for case in manifest['cases']:
        path=args.root/case['name']/name
        a,t=module.read(path);times.append(t)
        assert np.array_equal(a['mass'],initial['mass']) and np.all(a['w2']==0)
        assert all(np.isfinite(a[k]).all() for k in module.FIELDS)
        assert np.all(a['pressure']>0)&np.all(a['ie']>0)&np.all(a['volume']>0)
        dx=np.array([a['x'][804]-a['x'][803],a['y'][804]-a['y'][803]])
        dv=np.array([a['vx'][804]-a['vx'][803],a['vy'][804]-a['vy'][803]])
        gap=float(np.linalg.norm(dx));h=float(np.sqrt(.5*(ic['volume'][803]+ic['volume'][804])))
        err=a['pressure'][inner]/ic['pressure'][inner]-1
        row['cases'][case['name']]={'sha256':module.sha(path),'gap_over_initial_h':gap/h,
                                   'separation_rate':float(dx@dv/gap),
                                   'pressure804_over_t0':float(a['pressure'][804]/ic['pressure'][804]),
                                   'inner_pressure_error_mass_RMS':float(np.sqrt(np.average(err*err,weights=a['mass'][inner]))),
                                   'inner_thermal_energy':float(a['ie'][inner].sum())}
    assert max(times)-min(times)<1e-12,(name,times)
    row['time']=times[0];rows.append(row)
out={'observed_utc':datetime.now(timezone.utc).isoformat(),'scope':'Partial common-snapshot observations; neither completion nor time-convergence claim.',
     'source_time':t0,'requested_horizon':manifest['horizon'],'source_sha256':manifest['source_sha256'],
     'pair':[803,804],'rows':rows}
with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps({'snapshots':len(rows),'latest':rows[-1] if rows else None},indent=2))
