"""Audit real native EOS overflow injection, rollback and successful recovery."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE

root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64EOSInjectSmoke_20261029')
start=root/'khout.000022.dat';end=root/'khout.000024.dat';out=root/'injection_audit.json'
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,m0=ordered_global(start);b,m1=ordered_global(end,a['u4if']&ID_MASK);q=payload(end,len(b))
log=(root/'log').read_text();rho=b['mass'].astype(float)/q['v'];p=GAMMA_MINUS_ONE*q['u']/q['v']
def energy(s,u):
    x=s['x']-3.;y=s['y']-3.;m=s['mass'].astype(float)
    return float(np.sum(u+.5*m*(s['vx']**2+s['vy']**2)-m/np.sqrt(x*x+y*y+.01175**2),dtype=np.longdouble))
e0=energy(a,a['ie'].astype(float));de=(energy(b,q['u'])-e0)/abs(e0)
checks={'native_exit_zero':(root/'exit_status.txt').read_text().strip()=='0',
        'one_test_injection':log.count('[EOS_INJECT_TEST]')==1,
        'one_exact_rollback_and_cache_invalidation':log.count('reason=nonpositive_or_nonfinite_trial payload_counts_dtold_ledgers_exact=1 device_cache_invalidated=1')==1,
        'two_accepted_steps':log.count('[WALL_ACCEPT]')==2,
        'exact_endpoint':m0['t']==.0015284035789018555 and m1['t']==.0016,
        'fixed_masses_and_W0':bool(np.array_equal(a['mass'],b['mass']) and np.all(b['w2']==0)),
        'valid_double_payload':bool(np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(q['u']>0) and np.all(np.isfinite(q['v'])) and np.all(q['v']>0)),
        'float_mirrors_and_eos':bool(np.array_equal(q['u'].astype(np.float32),b['ie']) and np.array_equal(q['v'].astype(np.float32),b['volume']) and np.max(abs(b['den']/rho-1))<1e-6 and np.max(abs(b['pressure']/p-1))<1e-6),
        'no_energy_injection':all(k+'=0.000e+00' in log for k in ('floor_cum','de_cum','sfl_cum')),
        'restart_energy_relative_1e-9':abs(de)<1e-9}
result={'passed':all(checks.values()),'checks':checks,'restart_energy_relative_change':de,
        'source_sha256':sha(start),'endpoint_sha256':sha(end),'native_log_sha256':sha(root/'log'),
        'scope':'Real CUDA bounded-Q1 smoke using a verification-only executable with a deliberately invalid trial once. Does not validate impact, all numerical failures, long-time physics, or production performance.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
