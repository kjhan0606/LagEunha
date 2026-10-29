"""Rebase H200 thermal64 entropy on an H200 zero-time promoted state."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Coarse_H200_20261029/init')
half=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Half_H200_20261029')
ada=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/thermal64_init_t10_v1/khdiag.step020790.thermal64_initial.dat')
initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
parent=root/'khout.020790.dat';h200=root/'khdiag.step020790.thermal64_initial.dat'
comp=json.loads((half/'arithmetic_comparison.json').read_text())
output=root/'initial_rebase.json'
if output.exists():raise FileExistsError(output)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if (root/'exit_status.txt').read_text().strip()!='0' or sha(parent)!=comp['source_sha256'] or sha(ada)!=comp['promotion_sha256']:raise ValueError('Invalid source')
log=(root/'log').read_text()
if '[THERMAL64_INIT_ONLY] zero physical time advanced; step=20790 time=10' not in log:raise ValueError('Missing zero-time native marker')
p0,m0=ordered_global(parent);ids=p0['u4if']&ID_MASK;ic,_=ordered_global(initial,ids)
pa,ma=ordered_global(ada,ids);ph,mh=ordered_global(h200,ids);qa=payload(ada,len(ids));qh=payload(h200,len(ids))
if not (m0['t']==ma['t']==mh['t']==10.):raise ValueError('Mismatched times')
for p,q in [(pa,qa),(ph,qh)]:
    if not all(np.array_equal(p0[k],p[k]) for k in ('x','y','vx','vy','mass','w2','ie')):raise ValueError('Physical state changed during initialization')
    if not (np.all(q['tag']==TAG) and np.array_equal(q['u'],p0['ie'].astype(float)) and np.all(np.isfinite(q['v'])) and np.all(q['v']>0)):raise ValueError('Invalid thermal payload')
files=[p for p in half.glob('khout.*.dat') if p.name!='khout.020790.dat']
if len(files)!=1 or sha(files[0])!=comp['thermal64_endpoint_sha256']:raise ValueError('Invalid half endpoint')
end,me=ordered_global(files[0],ids);qe=payload(files[0],len(ids))
if me['t']!=10.1 or not np.all(qe['tag']==TAG):raise ValueError('Invalid half payload')
r0=np.hypot(ic['x']-3.,ic['y']-3.);vol_delta=np.log(qh['v']/qa['v']);weights=p0['mass'].astype(float)
base_a=np.log(qe['u']/qa['u'])+GAMMA_MINUS_ONE*np.log(qe['v']/qa['v'])
base_h=np.log(qe['u']/qh['u'])+GAMMA_MINUS_ONE*np.log(qe['v']/qh['v'])
regions={}
for name,mask in {'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.)}.items():
    w=weights[mask]
    regions[name]={'cells':int(mask.sum()),'ada_h200_initial_log_volume_abs_mass_mean':mean(abs(vol_delta[mask]),w),
      'ada_h200_initial_log_volume_abs_max':float(np.max(abs(vol_delta[mask]))),
      'half_entropy_abs_mass_mean_ada_reference':mean(abs(base_a[mask]),w),
      'half_entropy_abs_mass_mean_h200_reference':mean(abs(base_h[mask]),w),
      'half_entropy_signed_mass_mean_h200_reference':mean(base_h[mask],w)}
result={'passed':True,'regions':regions,'initial_h200_sha256':sha(h200),'initial_ada_sha256':sha(ada),
        'half_endpoint_sha256':sha(files[0]),'scope':'Zero-time initial comparison across Ada/H200 GPUs. H200 reference uses same candidate source and hardware type, not necessarily identical node/GPU UUID; no physical evolution occurs in the initialization run.'}
with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
