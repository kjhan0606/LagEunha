"""Validate early thermal64 Q1 against the preserved eta=.06 legacy control."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import GAMMA_MINUS_ONE, mean
from audit_kepler_q1_checkpoint_geometry_20261029 import audit as geometry_audit
from audit_kepler_q1_thermal64_promotion_20261029 import payload, TAG

sys.path.insert(0,'/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
from run_kepler_q1_orbital_cap_chain_20261029 import core_screen
root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64Early_Ada_20261029')
legacy_root=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
ic=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
cfg=json.loads((root/'contract.json').read_text());policy=json.loads((legacy_root/'contract.json').read_text())
out=root/'comparison.json'
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
seed=Path(cfg['source']);legacy=Path(cfg['legacy_reference'])
if sha(seed)!=cfg['source_sha256'] or sha(legacy)!=cfg['legacy_reference_sha256'] or sha(root/'eunha2.exe')!=cfg['binary_sha256'] or (root/'exit_status.txt').read_text().strip()!='0':raise ValueError('Frozen inputs or native exit invalid')
ends=[p for p in root.glob('khout.*.dat') if p.name!=seed.name]
if len(ends)!=1:raise ValueError('Expected one final checkpoint')
end=ends[0];a,_=ordered_global(ic);ids=a['u4if']&ID_MASK;l,lm=ordered_global(legacy,ids);p,pm=ordered_global(end,ids);q=payload(end,len(ids))
if lm['t']!=.4125 or pm['t']!=.4125 or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Changed IDs, time, mass or weights')
if not (np.all(q['tag']==TAG) and np.all(np.isfinite(q['u'])) and np.all(q['u']>0) and np.all(np.isfinite(q['v'])) and np.all(q['v']>0)):raise ValueError('Invalid double payload')
gate_path=root/'checkpoint_gate.json'
if gate_path.exists():raise FileExistsError(gate_path)
proc=subprocess.run([sys.executable,str(Path(__file__).with_name('audit_kepler_q1_bounded_gate_20261029.py')),str(ic),str(root),'--time','0.4125','--output',str(gate_path)],capture_output=True,text=True)
(root/'checkpoint_gate.log').write_text(proc.stdout+proc.stderr)
gate=json.loads(gate_path.read_text())
geometry=geometry_audit(ic,end)
with (root/'independent_geometry.json').open('x') as f:json.dump(geometry,f,indent=2);f.write('\n')
core=core_screen(ic,end,geometry,policy)
with (root/'core_screen.json').open('x') as f:json.dump(core,f,indent=2);f.write('\n')
r0=np.hypot(a['x']-3.,a['y']-3.);mass=a['mass'].astype(float)
regions={'core':r0<.1,'adjacent':(r0>=.1)&(r0<.5),'disk':(r0>=.5)&(r0<2.)}
eps=.01175;phase=.4125*(r0*r0+eps*eps)**(-.75);xx=a['x']-3.;yy=a['y']-3.
xref=3.+xx*np.cos(phase)-yy*np.sin(phase);yref=3.+xx*np.sin(phase)+yy*np.cos(phase)
def metrics(state,u,v,mask):
    w=mass[mask];dx=state['x'][mask]-3.;dy=state['y'][mask]-3.;l0=np.sum(w*(xx[mask]*a['vy'][mask]-yy[mask]*a['vx'][mask]),dtype=np.longdouble)
    lt=np.sum(w*(dx*state['vy'][mask]-dy*state['vx'][mask]),dtype=np.longdouble)
    ent=np.log(u[mask]/a['ie'][mask].astype(float))+GAMMA_MINUS_ONE*np.log(v[mask]/a['volume'][mask].astype(float))
    return {'normalized_orbit_error_mass_mean':mean(np.hypot(state['x'][mask]-xref[mask],state['y'][mask]-yref[mask])/np.maximum(r0[mask],eps),w),
            'density_change_abs_mass_mean':mean(abs(state['den'][mask].astype(float)/a['den'][mask].astype(float)-1),w),
            'entropy_change_abs_mass_mean':mean(abs(ent),w),'entropy_change_signed_mass_mean':mean(ent,w),
            'angular_momentum_relative_change':float(lt/l0-1) if abs(l0)>0 else None}
data={}
for name,mask in regions.items():
    data[name]={'cells':int(mask.sum()),'legacy':metrics(l,l['ie'].astype(float),l['volume'].astype(float),mask),
                'thermal64':metrics(p,q['u'],q['v'],mask),
                'position_difference_mass_mean':mean(np.hypot(p['x'][mask]-l['x'][mask],p['y'][mask]-l['y'][mask]),mass[mask])}
def energy(state,u):
    x=state['x']-3.;y=state['y']-3.;return float(np.sum(u+.5*mass*(state['vx']**2+state['vy']**2)-mass/np.sqrt(x*x+y*y+eps*eps),dtype=np.longdouble))
e0=energy(a,a['ie'].astype(float));el=energy(l,l['ie'].astype(float));et=energy(p,q['u'])
result={'passed':bool(proc.returncode==0 and len(gate['checks'])==28 and all(gate['checks'].values()) and geometry['passed'] and core['passed']),
        'checks':{'legacy_style_28_gates':bool(proc.returncode==0 and len(gate['checks'])==28 and all(gate['checks'].values())),
                  'independent_geometry':bool(geometry['passed']),'core_accuracy_screen':bool(core['passed'])},
        'regions':data,'energy_relative_to_initial':{'legacy':(el-e0)/abs(e0),'thermal64':(et-e0)/abs(e0)},
        'candidate_endpoint_sha256':sha(end),'legacy_endpoint_sha256':sha(legacy),
        'scope':'Two native eta=.06 runs from the same accepted early restart to t=.4125, both on Ada; candidate and legacy binaries differ. Not a t6.5 or long-time test. Thermal entropy reference uses legacy t0 float volume, so initial representation is retained.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
