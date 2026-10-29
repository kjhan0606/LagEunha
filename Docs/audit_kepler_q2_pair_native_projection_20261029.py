"""Project actual native diagnostic pressure accelerations onto selected pairs."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK

root=Path('/gpfs/kjhan/LagForce/Kepler_Q2_PairRHS_Ada_20261029')
cfg=json.loads((root/'probe.json').read_text());audit=json.loads((root/'audit.json').read_text())
log=(root/'log').read_text(errors='replace')
if not audit['passed'] or '[WALL_RHS_PROBE] computed K1 for all 55351 owned sites; diagnostic only, no accepted evolution' not in log:raise ValueError('Selected RHS audit or no-evolution marker missing')
p,meta=ordered_global(Path(cfg['endpoint']));ids=p['u4if']&ID_MASK
if meta['t']!=cfg['time']:raise ValueError('Unexpected time')
with (root/'rk4_stage_r0.tsv').open() as f: rows={int(r['id']):r for r in csv.DictReader(f,delimiter='\t')}
pairs=[(31088,30249),(29637,23949),(25319,23161)];results=[]
for u,v in pairs:
    i,j=int(np.searchsorted(ids,u)),int(np.searchsorted(ids,v))
    if ids[i]!=u or ids[j]!=v or u not in rows or v not in rows:raise ValueError('Missing selected ID')
    delta=np.array([p['x'][j]-p['x'][i],p['y'][j]-p['y'][i]]);dist=float(np.linalg.norm(delta));n=delta/dist
    ai=np.array([float(rows[u]['ax_hydro']),float(rows[u]['ay_hydro'])])
    aj=np.array([float(rows[v]['ax_hydro']),float(rows[v]['ay_hydro'])])
    vel=np.array([p['vx'][j]-p['vx'][i],p['vy'][j]-p['vy'][i]])
    results.append({'id_i':u,'id_j':v,'separation':dist,'relative_speed_along_pair':float(np.dot(vel,n)),
                    'native_pressure_relative_acceleration_along_pair':float(np.dot(aj-ai,n))})
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
result={'passed':all(r['relative_speed_along_pair']<0 and r['native_pressure_relative_acceleration_along_pair']>0 for r in results),
        'time':meta['t'],'pairs':results,'audit_sha256':sha(root/'audit.json'),
        'native_tsv_sha256':sha(root/'rk4_stage_r0.tsv'),'input_sha256':sha(Path(cfg['endpoint'])),
        'script_sha256':sha(Path(__file__)),
        'scope':'Three deliberately selected close pairs, one frozen native K1 from wall-event v2 binary on Ada; no physical evolution. The acceleration is total pressure contribution, not a single face. No population-wide, causal or long-time conclusion.'}
with (root/'native_pair_projection.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
