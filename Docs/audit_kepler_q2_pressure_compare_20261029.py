"""Compare independently completed Q2 pressure-sensitivity array elements."""
import hashlib
import json
import re
import sys
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global,ID_MASK
from audit_kepler_native_radial_20261026 import energy

root=Path(sys.argv[1]);out=Path(sys.argv[2]);contract=json.loads((root/'contract.json').read_text())
original=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/t001_v2/khout.000000.dat')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
if sha(original)!='a4d557ae5837e933c8c2c8c335f91e5c97d36f9a3bdb8acd5568b0499140b665':raise ValueError('original IC changed')
orig,meta=ordered_global(original);ids=orig['u4if']&ID_MASK
xy0=np.column_stack((orig['x']-3,orig['y']-3));r0=np.linalg.norm(xy0,axis=1)
phi=np.arctan2(xy0[:,1],xy0[:,0]);rad=np.round(r0,10);sel=(r0>=.5)&(r0<.75)
rings=[]
for r in np.unique(rad[sel]):
    ix=np.flatnonzero(sel&(rad==r));ix=ix[np.argsort(phi[ix])]
    if len(ix)<24:raise ValueError('invalid original ring')
    rings.append(ix)
L0=xy0[:,0]*orig['vy']-xy0[:,1]*orig['vx']
def band_power(p):
    x=p['x']-3;y=p['y']-3;L=x*p['vy']-y*p['vx']-L0;power=0.
    for ix in rings:
        z=L[ix]-np.mean(L[ix]);f=np.fft.rfft(z)/len(ix);q=abs(f)**2
        q[1:]*=2
        if len(ix)%2==0:q[-1]*=.5
        freq=np.fft.rfftfreq(len(ix));mass=float(np.sum(orig['mass'][ix]))
        power+=mass*np.sum(q[(freq>=.2)&(freq<.3)])
    return float(power)
def find_endpoint(d):
    if (d/'exit_status.txt').read_text().strip()!='0':raise ValueError('solver failed '+str(d))
    found=[]
    for path in d.glob('khout.*.dat'):
        if path.name=='khout.000001.dat':continue
        if not re.fullmatch(r'khout\.\d{6}\.dat',path.name):continue
        p,m=ordered_global(path,ids)
        if abs(m['t']-contract['target_time'])<1e-9:found.append((path,p,m))
    if len(found)!=1:raise ValueError('expected one complete target checkpoint in '+str(d))
    return found[0]
cases=[]
for factor in (1,4):
    d=root/f'x{factor}';start,sm=ordered_global(d/'khout.000001.dat',ids);path,p,mt=find_endpoint(d)
    if sm['t']!=30 or not np.array_equal(start['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('mass/W0 invariant')
    for field in ('x','y','vx','vy','ie','pressure','den','volume'):
        if not np.all(np.isfinite(p[field])):raise ValueError('nonfinite '+field)
    if np.any(p['ie']<=0) or np.any(p['pressure']<=0) or np.any(p['den']<=0):raise ValueError('nonpositive state')
    a=band_power(start);b=band_power(p)
    if not a>0 or not b>0:raise ValueError('no mode power')
    gamma=float(np.float32(1.666667));K=p['pressure'].astype(float)/p['den'].astype(float)**gamma
    K0=start['pressure'].astype(float)/start['den'].astype(float)**gamma
    logK=np.abs(np.log(K[sel]/K0[sel]));mass=p['mass'][sel].astype(float)
    e0=energy(start,'baseline');e1=energy(p,'baseline')
    cases.append({'pressure_scale':factor,'endpoint':str(path),'endpoint_sha256':sha(path),
                  'time':mt['t'],'band_power_initial':a,'band_power_final':b,
                  'band_log_amplitude_growth_per_time':float(.5*np.log(b/a)/(mt['t']-sm['t'])),
                  'inner_abs_log_entropy_change_mass_mean':float(np.average(logK,weights=mass)),
                  'relative_total_energy_change':float((e1-e0)/abs(e0))})
result={'cases':cases,'source_contract_sha256':sha(root/'contract.json'),'audit_script_sha256':sha(Path(__file__)),
        'scope':'Two independent fixed-mass W0 native runs from identical t30 geometry/velocities, with pressure/internal energy scaled at restart. x4 is a deliberate counterfactual, not the original Kepler physics. Short-mode response is not a long-time stability or convergence proof. Compare own energy references; threshold acceptance requires additional wall/geometry/log audits.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
