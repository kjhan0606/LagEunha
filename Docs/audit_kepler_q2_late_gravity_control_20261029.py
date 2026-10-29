"""Gravity-only RK4 attribution on saved Q2 t30 to preimpact state."""
import hashlib
import json
from pathlib import Path
import re
import socket
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_q1_core_rk4_control_20261029 import advance
from audit_kepler_native_radial_20261026 import mean

initial=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/t001_v2/khout.000000.dat')
start=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to30_a10_v1/khout.005247.dat')
end=Path('/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/control_capture/khout.011580.dat')
log=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/log')
out=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs/KEPLER_Q2_LATE_GRAVITY_CONTROL_2026-10-29.json')
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run ODE control on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
clock=time.monotonic();a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
b,mb=ordered_global(start,ids);c,mc=ordered_global(end,ids)
if ma['t']!=0 or mb['t']!=30. or abs(mc['t']-48.37534965270793)>1e-12 or not (np.array_equal(a['mass'],b['mass']) and np.array_equal(a['mass'],c['mass'])) or np.any(c['w2']!=0):raise ValueError('Unexpected native states')
rows=[]
for line in log.open(errors='replace'):
    q=re.search(r'Time=\s*\S+\s+step=\s*(\d+)\s+dt=\s*(\S+)',line)
    if q:rows.append((int(q[1]),float(q[2])))
if len(rows)!=11579 or [q[0] for q in rows]!=list(range(2,11581)):raise ValueError('Incomplete timestep log')
hs=np.array([q[1] for q in rows]);interval=mc['t']-mb['t'];printed_sum=float(hs.sum())
if np.any(hs<=0) or abs(printed_sum-interval)>1e-5:raise ValueError('Printed timestep sequence inconsistent with interval')
hs*=interval/printed_sum
r0=np.hypot(a['x']-3.,a['y']-3.);mask=(r0>=.5)&(r0<.75);w=a['mass'][mask].astype(float)
q0=np.column_stack((b['x'][mask]-3.,b['y'][mask]-3.,b['vx'][mask],b['vy'][mask]))
native=np.column_stack((c['x'][mask]-3.,c['y'][mask]-3.))
results=[];states=[]
for n in (1,2):
    q=q0.copy()
    for h in hs:
        for _ in range(n):q=advance(q,h/n)
    if not np.isfinite(q).all():raise ValueError('Invalid gravity solution')
    states.append(q)
    diff=np.linalg.norm(q[:,:2]-native,axis=1)
    results.append({'substeps':n,'mass_mean_position_difference_from_native':mean(diff,w),
                    'cell_p99_position_difference_from_native':float(np.percentile(diff,99))})
diff=np.linalg.norm(states[0][:,:2]-states[1][:,:2],axis=1)
result={'time_start':mb['t'],'time_end':mc['t'],'cells':int(mask.sum()),'logged_steps':len(hs),
        'printed_timestep_sum_minus_exact_interval':printed_sum-interval,'sum_rescale_factor':interval/printed_sum,
        'cases':results,'one_vs_two_substeps_mass_mean_position_difference':mean(diff,w),
        'one_vs_two_substeps_p99_position_difference':float(np.percentile(diff,99)),
        'source_sha256':{str(p):sha(p) for p in (initial,start,end,log,Path(__file__))},
        'elapsed_seconds':time.monotonic()-clock,
        'scope':'Gravity-only continuations from the actual t30 native state. Printed h values are rounded and timing alignment is approximate; sequence is globally rescaled to the exact saved interval. Difference from native includes pressure and all other dynamics plus approximation, not an isolated pressure force. One-vs-two substeps assesses gravity RK4 sensitivity, not a rigorous error bound or a full coupled gas calculation.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
