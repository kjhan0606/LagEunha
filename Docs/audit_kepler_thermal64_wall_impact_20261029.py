"""Audit native Q2 collision replay and authoritative double thermal energy."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_kepler_wall_events_20261029 import analyze
from audit_kepler_native_radial_20261026 import ordered, potential
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Thermal64WallImpact_20261029')
r=analyze(root)
a,_=ordered(Path(r['initial'])); b,_=ordered(Path(r['endpoint']))
q=payload(Path(r['endpoint']),len(b))
def energy(p,u):
    mass=p['mass'].astype(float);x=p['x'].astype(float)-3.;y=p['y'].astype(float)-3.
    vx=p['vx'].astype(float);vy=p['vy'].astype(float)
    return float(np.sum(u+.5*mass*(vx*vx+vy*vy)+mass*potential(np.hypot(x,y),'baseline'),dtype=np.longdouble))
e0=energy(a,a['ie'].astype(float)); e1=energy(b,q['u'])
rho=b['mass'].astype(float)/q['v']; pressure=(float(np.float32(1.666667))-1)*q['u']/q['v']
checks={'all_double_payload_tagged_finite_positive':bool(np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(np.isfinite(q['v'])) and np.all(q['u']>0) and np.all(q['v']>0)),
        'float_mirrors_exact':bool(np.array_equal(b['ie'],q['u'].astype(np.float32)) and np.array_equal(b['volume'],q['v'].astype(np.float32))),
        'eos_relative_1e-6':bool(np.max(abs(b['den']/rho-1))<1e-6 and np.max(abs(b['pressure']/pressure-1))<1e-6),
        'double_energy_relative_1e-7':abs((e1-e0)/abs(e0))<1e-7}
r['checks'].update(checks);r['passed']=all(r['checks'].values())
r['restart_double_energy_relative_change']=(e1-e0)/abs(e0)
r['native_log_sha256']=hashlib.sha256((root/'log').read_bytes()).hexdigest()
r['scope']='Actual localized Q2 wall collision with thermal64, not Q1 or long-time stability; momentum and energy audits do not establish temporal convergence.'
out=root/'thermal_impact_audit.json'
if out.exists():raise FileExistsError(out)
with out.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
print(json.dumps(r,indent=2))
if not r['passed']:raise SystemExit(1)
