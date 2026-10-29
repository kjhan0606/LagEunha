"""Audit thermal64 native serial/OMP guards at the same physical wall impact."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from audit_kepler_wall_events_interleaved_20261029 import analyze
from audit_kepler_native_radial_20261026 import ordered, potential
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG
from audit_kepler_restart_bytes_20261029 import compare

root=Path(sys.argv[1]); cases={}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def energy(p,u):
    m=p['mass'].astype(float);x=p['x'].astype(float)-3;y=p['y'].astype(float)-3
    vx=p['vx'].astype(float);vy=p['vy'].astype(float)
    return float(np.sum(u+.5*m*(vx*vx+vy*vy)+m*potential(np.hypot(x,y),'baseline'),dtype=np.longdouble))
for label in ('thermal_impact_A','thermal_impact_B'):
    d=root/label;r=analyze(d);a,_=ordered(Path(r['initial']));b,_=ordered(Path(r['endpoint']))
    q=payload(Path(r['endpoint']),len(b));e0=energy(a,a['ie'].astype(float));e1=energy(b,q['u'])
    rho=b['mass'].astype(float)/q['v'];pressure=(float(np.float32(1.666667))-1)*q['u']/q['v']
    checks={'all_double_payload_tagged_finite_positive':bool(np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(np.isfinite(q['v'])) and np.all(q['u']>0) and np.all(q['v']>0)),
            'float_mirrors_exact':bool(np.array_equal(b['ie'],q['u'].astype(np.float32)) and np.array_equal(b['volume'],q['v'].astype(np.float32))),
            'eos_relative_1e-6':bool(np.max(abs(b['den']/rho-1))<1e-6 and np.max(abs(b['pressure']/pressure-1))<1e-6),
            'double_energy_relative_1e-7':abs((e1-e0)/abs(e0))<1e-7,
            'thermal_initialization_log':'[THERMAL64_INIT] accepted U promoted exactly' in (d/'log').read_text(errors='replace')}
    r['checks'].update(checks);r['passed']=all(r['checks'].values());r['energy64_relative_change']=(e1-e0)/abs(e0)
    r['log_sha256']=sha(d/'log');cases[label]=r
parity=compare(*[Path(cases[label]['endpoint']) for label in ('thermal_impact_A','thermal_impact_B')])
checks={'both_physical_and_thermal_gates':all(v['passed'] for v in cases.values()),
        'both_real_impacts':all(len(v['impacts'])>0 for v in cases.values()),
        'both_real_rollbacks':all(v['rollbacks']>0 for v in cases.values()),
        'full_record_footer_parity':parity['passed']}
result={'passed':all(checks.values()),'checks':checks,'cases':cases,'parity':parity,
        'audit_sha256':sha(Path(__file__)),
        'scope':'Same thermal64 executable and legacy pre-impact input, initialized natively at t48.3753, serial vs OpenMP wall guard; one physical collision. Shared Ada GPU: no timing qualification, independent conversion snapshot, MPI, temporal convergence or long-time stability claim.'}
with (root/'thermal_guard_impact_audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
