"""Compare restarted one-step evolution with the retained continuous trajectory."""
import json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_quadratic_restart_evolution_20261018'
summary=root/'summary.json'
if not summary.exists():
    print(json.dumps(dict(status='pending',scope=__doc__)));sys.exit(2)
s=json.loads(summary.read_text())
if s['status']!='completed':
    print(json.dumps(dict(status=s['status'],events=s['events'],scope=__doc__)));sys.exit(1)
continuous=base/'Codex_continuous_quadratic_evolution_dt5e5_20261018'
original=np.load(continuous/'step0001.npz');loaded=np.load(root/'initial.npz')
for key in ('state','factors','mass','ids','time'):
    assert original[key].tobytes()==loaded[key].tobytes(),key
a=np.load(continuous/'step0002.npz');b=np.load(root/'step0001.npz')
assert float(a['time'])==float(b['time'])==s['completed_time']
assert np.array_equal(a['mass'],b['mass']) and np.array_equal(a['ids'],b['ids'])
rows={}
for key in ('state','factors'):
    assert np.isfinite(b[key]).all()
    error=abs(b[key]-a[key]);tolerance=1e-14+1e-10*abs(a[key])
    rows[key]=dict(bitwise_equal=a[key].tobytes()==b[key].tobytes(),
        maximum_absolute_difference=float(error.max()),maximum_tolerance_fraction=float((error/tolerance).max()))
mass=a['mass'];aa=a['state'][:-4].reshape(-1,20);bb=b['state'][:-4].reshape(-1,20)
dv=(bb[:,9:11]-aa[:,9:11])/mass[:,None]
report=dict(scope=__doc__,status='completed',time=float(a['time']),initial_bitwise_equal=True,comparisons=rows,
    velocity_mass_RMS=float(np.sqrt(np.average(np.sum(dv*dv,axis=1),weights=mass))),
    limits='Legacy checkpoint lacked nonlinear seed. Numerical equivalence tolerance 1e-14+1e-10*abs(reference); not bitwise equivalence unless explicitly reported. No production restart or long-term claim.')
print(json.dumps(report,indent=2))
assert all(row['maximum_tolerance_fraction']<1 for row in rows.values())
