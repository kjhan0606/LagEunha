"""Inspect one accepted bounded thermal64 step after forced wall rollback."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG
from audit_kepler_native_radial_20261026 import mean, potential

root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallSmoke_20261029')
start=root/'khout.000022.dat'; end=root/'khout.000023.dat'
p0,m0=ordered_global(start); p1,m1=ordered_global(end,p0['u4if']&ID_MASK)
q=payload(end,len(p1)); log=(root/'log').read_text()
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def energy(p,u):
    r=np.hypot(p['x']-3.,p['y']-3.)
    mass=p['mass'].astype(float); kinetic=.5*mass*(p['vx']**2+p['vy']**2)
    return float(np.sum(u+kinetic+mass*potential(r,'baseline'),dtype=np.longdouble))
e0=energy(p0,p0['ie'].astype(float)); e1=energy(p1,q['u'])
rho=p1['mass'].astype(float)/q['v']; gm1=float(np.float32(1.666667))-1.
pressure=gm1*q['u']/q['v']
checks={'expected_exact_times':m0['t']==0.0015284035789018555 and m1['t']==.0016,
        'same_fixed_masses_and_W0':bool(np.array_equal(p0['mass'],p1['mass']) and np.all(p1['w2']==0)),
        'all_tagged_positive_finite':bool(np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(np.isfinite(q['v'])) and np.all(q['u']>0) and np.all(q['v']>0)),
        'float_mirrors':bool(np.array_equal(p1['ie'],q['u'].astype(np.float32)) and np.array_equal(p1['volume'],q['v'].astype(np.float32))),
        'EOS_relative_1e-6':bool(np.max(abs(p1['den']/rho-1))<1e-6 and np.max(abs(p1['pressure']/pressure-1))<1e-6),
        'one_forced_exact_payload_and_cache_rollback':log.count('reason=verification_rejection payload_counts_dtold_ledgers_exact=1 device_cache_invalidated=1')==1,
        'one_accepted_two_attempts_no_impacts':len(re.findall(r'\[WALL_ACCEPT\].*attempts=2 impacts=0',log))==1,
        'no_floor_injection':bool('floor_cum=0.000e+00' in log and 'de_cum=0.000e+00' in log and 'sfl_cum=0.000e+00' in log),
        'total_energy_relative_1e-9':abs((e1-e0)/e0)<1e-9}
out={'passed':all(checks.values()),'checks':checks,'start_time':m0['t'],'end_time':m1['t'],
     'sites':len(p1),'energy_relative_change':(e1-e0)/abs(e0),
     'u64_minus_float_mirror_max_relative':float(np.max(abs(q['u']/p1['ie']-1))),
     'source_sha256':sha(start),'endpoint_sha256':sha(end),'native_log_sha256':sha(root/'log'),
     'scope':'One short no-impact step with a forced rejection; not a physical convergence, impact, restart-continuity, multi-rank or long-time validation.'}
path=root/'audit.json'
if path.exists():raise FileExistsError(path)
with path.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
if not out['passed']:raise SystemExit(1)
