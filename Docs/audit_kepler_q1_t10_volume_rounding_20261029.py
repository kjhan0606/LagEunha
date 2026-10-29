"""Independent polygon geometry bounds endpoint volume representation effects."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_checkpoint_geometry_20261029 import independent
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean
from audit_kepler_q1_thermal64_promotion_20261029 import payload

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/controls')
promotion=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/thermal64_init_t10_v1/khdiag.step020790.thermal64_initial.dat')
cfg=json.loads((root/'stage.json').read_text());source=Path(cfg['parent']['endpoint']);initial=Path(cfg['initial'])
out=root/'independent_volume_rounding.json'
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(source)!=cfg['parent']['sha256'] or sha(promotion)!='0456a8c6dafb3f7ee7e681089c8c475462a4fe7e4ee062a19cf19f4f2db52381':raise ValueError('Unexpected t10 state')
p0,m0=ordered_global(source);ids=p0['u4if']&ID_MASK
ic,_=ordered_global(initial,ids);prom,pm=ordered_global(promotion,ids);q0=payload(promotion,len(ids))
if m0['t']!=pm['t'] or m0['t']!=10. or not np.array_equal(p0['ie'],prom['ie']):raise ValueError('Mismatch promotion')
r0=np.hypot(ic['x']-3.,ic['y']-3.)
regions={'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.)}
cases={}
for name in ['coarse','half','quarter']:
    d=root/name;gate=json.loads((d/'checkpoint_gate.json').read_text());end=Path(gate['endpoint'])
    done=json.loads((d/'completed.json').read_text())
    if len(gate['checks'])!=28 or not all(gate['checks'].values()) or sha(end)!=done['sha256']:raise ValueError('Invalid endpoint '+name)
    p,m=ordered_global(end,ids)
    if m['t']!=10.1 or not np.array_equal(p0['mass'],p['mass']):raise ValueError('Invalid time or mass')
    vol,geometry=independent(np.column_stack((p['x'],p['y'])))
    volume_error=p['volume'].astype(float)/vol-1.
    if np.max(abs(volume_error))>=1e-6:raise ValueError('Endpoint independent geometry failed')
    native=np.log(p['ie'].astype(float)/p0['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(p['volume'].astype(float)/p0['volume'].astype(float))
    corrected=np.log(p['ie'].astype(float)/p0['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(vol/q0['v'])
    regions_out={}
    for key,mask in regions.items():
        mass=p0['mass'][mask].astype(float)
        regions_out[key]={'cells':int(mask.sum()),'native_abs_mass_mean':mean(abs(native[mask]),mass),
                          'independent_volume_abs_mass_mean':mean(abs(corrected[mask]),mass),
                          'volume_representation_delta_abs_mass_mean':mean(abs((native-corrected)[mask]),mass),
                          'volume_representation_delta_abs_max':float(np.max(abs((native-corrected)[mask])))}
    cases[name]={'endpoint_sha256':sha(end),'geometry':geometry,'native_volume_max_relative_error':float(np.max(abs(volume_error))),
                 'regions':regions_out}
result={'passed':True,'cases':cases,'source_sha256':sha(source),'promotion_sha256':sha(promotion),
        'scope':'Readout volume representation only: independent double geometry at saved endpoint and t10. It does not isolate accumulated per-step volume rounding, internal energy error, or physically incorrect trajectories.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
