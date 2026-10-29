"""Radial localization of t10-to-10.1 entropy increments in existing controls."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean

root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/controls')
cfg=json.loads((root/'stage.json').read_text());source=Path(cfg['parent']['endpoint']);initial=Path(cfg['initial'])
output=root/'entropy_locality.json'
if output.exists():raise FileExistsError(output)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(source)!=cfg['parent']['sha256']:raise ValueError('Changed parent')
p0,m0=ordered_global(source);ids=p0['u4if']&ID_MASK;ic,_=ordered_global(initial,ids)
r0=np.hypot(ic['x']-3.,ic['y']-3.);m=p0['mass'].astype(float);states={};hashes={}
for name in ('coarse','half','quarter'):
    d=root/name;gate=json.loads((d/'checkpoint_gate.json').read_text());done=json.loads((d/'completed.json').read_text());end=Path(gate['endpoint'])
    if len(gate['checks'])!=28 or not all(gate['checks'].values()) or sha(end)!=done['sha256']:raise ValueError('Invalid '+name)
    p,meta=ordered_global(end,ids)
    if meta['t']!=10.1 or not np.array_equal(p0['mass'],p['mass']):raise ValueError('Wrong time or mass')
    states[name]=np.log(p['ie'].astype(float)/p0['ie'].astype(float))+GAMMA_MINUS_ONE*np.log(p['volume'].astype(float)/p0['volume'].astype(float))
    hashes[name]=sha(end)
def weighted_corr(x,y,w):
    x=x-mean(x,w);y=y-mean(y,w)
    den=np.sqrt(mean(x*x,w)*mean(y*y,w))
    return mean(x*y,w)/den if den>0 else None
out=[]
for lo,hi in [(0,.1),(.1,.5),(.5,.75),(.75,1.),(1.,1.25),(1.25,1.5),(1.5,1.75),(1.75,2.),(2.,2.65)]:
    mask=(r0>=lo)&(r0<hi);w=m[mask]
    if not len(w):raise ValueError('Empty radial band')
    vals={}
    for name in states:
        a=states[name][mask]
        vals[name]={'absolute_mass_mean':mean(abs(a),w),'signed_mass_mean':mean(a,w),
                    'positive_mass_fraction':mean((a>0).astype(float),w),'cell_abs_p99':float(np.percentile(abs(a),99))}
    a,b=states['half'][mask],states['quarter'][mask]
    out.append({'r0_lower':lo,'r0_upper':hi,'cells':int(mask.sum()),'cases':vals,
                'half_quarter_weighted_correlation':weighted_corr(a,b,w),
                'half_quarter_same_sign_mass_fraction':mean((np.sign(a)==np.sign(b)).astype(float),w)})
result={'source_sha256':sha(source),'endpoint_sha256':hashes,'regions':out,
        'scope':'Initial-ID radial bands and observed signed log-entropy increments. Correlations are descriptive, not proof of rounding, physical heating, causation or temporal convergence.'}
with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
