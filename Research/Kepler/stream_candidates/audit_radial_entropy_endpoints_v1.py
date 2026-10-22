"""Signed endpoint identity from preserved states, with no inferred rates."""
from pathlib import Path
import csv, gzip, json
import numpy as np
from audit_snapshot import read_snapshot
root=Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
base=root/'radial_resolution_20261021_v2'
gamma=json.loads((base/'manifest.json').read_text())['gamma']
mask=np.uint64((1<<56)-1)
order=lambda p:np.argsort(p['u4if']&mask)
result={'scope':'saved positive endpoint identity; no interpolation of missing stage rates',
        'identity':'delta log A = delta log U + (gamma-1) delta log V at fixed mass','cases':{}}
for label in ('uniform','graded'):
    initial,_=read_snapshot(base/label/'khout.000001.dat');initial=initial[order(initial)]
    radius=np.hypot(initial['x']-12,initial['y']-12)
    manifest=json.loads((base/'native_v2'/label/'manifest.json').read_text())
    states=[]
    for stage in manifest['stages']:
        if 'physical_audit' in stage:states.append(Path(stage['checkpoint']))
        elif 'normal_step_start' in stage:states.append(Path(stage['normal_step_start']))
    rows=[];seen=set()
    for path in states:
        p,meta=read_snapshot(path);p=p[order(p)]
        if meta['t'] in seen:continue
        seen.add(meta['t'])
        if not np.array_equal(p['mass'],initial['mass']) or np.any(p['ie']<=0) or np.any(p['volume']<=0):
            raise ValueError('requires positive fixed-mass endpoints')
        du=np.log(p['ie'].astype(float)/initial['ie'])
        dv=(gamma-1)*np.log(p['volume'].astype(float)/initial['volume']);da=du+dv
        output=base/f'entropy_endpoints_{label}_t{meta["t"]:.12g}.csv.gz'
        with gzip.open(output,'wt',compresslevel=1,newline='') as f:
            writer=csv.writer(f);writer.writerow(['id','initial_radius','delta_log_U','(gamma-1)_delta_log_V','delta_log_A'])
            writer.writerows(zip((p['u4if']&mask).tolist(),radius.tolist(),du.tolist(),dv.tolist(),da.tolist()))
        row={'snapshot':str(path),'time':meta['t'],'ledger':str(output),'radial':[]}
        edges=[0,.047,.1,.2,.5,1,2,4,8,18]
        for lo,hi in zip(edges[:-1],edges[1:]):
            selected=(radius>=lo)&(radius<hi)
            if not selected.any():continue
            mass=p['mass'][selected];avg=lambda x:float(np.sum(mass*x[selected])/mass.sum())
            row['radial'].append(dict(Rlo=lo,Rhi=hi,n=int(selected.sum()),mean_delta_log_U=avg(du),
                mean_volume_log_term=avg(dv),mean_delta_log_A=avg(da),
                min_delta_log_A=float(da[selected].min()),max_delta_log_A=float(da[selected].max())))
        rows.append(row)
    result['cases'][label]=rows
(base/'entropy_endpoint_ledger.json').write_text(json.dumps(result,indent=2)+'\n')
print('Signed positive endpoint ledger complete; no new native RHS.')
