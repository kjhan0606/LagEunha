"""Compare original and eta=.06 native Q1 states at the same physical time."""
import hashlib
import json
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_native_radial_20261026 import mean, energy
from audit_kepler_q1_restart_dt_bundle_20261029 import region_metrics
from audit_kepler_q1_checkpoint_geometry_20261029 import audit as geometry_audit

bounded=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029')
chain=bounded/'chain_h200_v2';capped=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029')
original=bounded/'ic_v1/khout.000001.dat';ref=chain/'to0p825/khout.001218.dat'
out=capped/'to0p825/eta006_vs_uncapped.json'
if out.exists():raise FileExistsError(out)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
state=json.loads((chain/'chain_state.json').read_text());exp=json.loads((capped/'chain_state.json').read_text())
if state['completed'][0]['time']!=.825 or sha(ref)!=state['completed'][0]['sha256']:raise ValueError('Changed uncapped reference')
stage=capped/'to0p825';receipt=json.loads((stage/'restart_receipt.json').read_text())
gate=json.loads((stage/'checkpoint_gate.json').read_text());core=json.loads((stage/'core_screen.json').read_text())
if (stage/'exit_status.txt').read_text().strip()!='0' or len(gate['checks'])!=28 or not all(gate['checks'].values()):raise ValueError('Native run or original checkpoint gate failed')
if sha(Path(receipt['source']))!=receipt['sha256'] or receipt['target_time']!=.825:raise ValueError('Invalid capped restart provenance')
ends=[p for p in stage.glob('khout.*.dat') if p.name!=Path(receipt['source']).name]
if len(ends)!=1:raise ValueError('Ambiguous capped endpoint')
end=ends[0]
entries=[v for v in exp['completed'] if v['time']==.825]
if entries and (len(entries)!=1 or sha(end)!=entries[0]['sha256']):raise ValueError('Changed certified endpoint')
a,m0=ordered_global(original);ids=a['u4if']&ID_MASK;b,mb=ordered_global(ref,ids);c,mc=ordered_global(end,ids)
if m0['t']!=0 or mb['t']!=mc['t'] or mc['t']!=.825 or not np.array_equal(a['mass'],b['mass']) or not np.array_equal(a['mass'],c['mass']):raise ValueError('Different physical time, IDs or masses')
if np.any(b['w2']!=0) or np.any(c['w2']!=0):raise ValueError('Expected W0')
r0=np.hypot(a['x']-3.,a['y']-3.);m=a['mass'].astype(float)
geo_reference=geometry_audit(original,ref)
geo_capped=json.loads((capped/'to0p825/independent_geometry.json').read_text())
if not geo_reference['passed'] or not geo_capped['passed']:raise ValueError('Independent geometry failed')
regions={}
for name,lo,hi in [('core',0.,.1),('adjacent',.1,.5),('disk',.5,2.)]:
    mask=(r0>=lo)&(r0<hi);w=m[mask]
    def geometry(g):return next(row for row in g['regions'] if row['rlo']==lo and row['rhi']==hi)['measures']['native_density_relative_to_circular_geometry']
    regions[name]={'cells':int(mask.sum()),'uncapped':region_metrics(a,a,b,mask,.825),'eta006':region_metrics(a,a,c,mask,.825),
                   'uncapped_density_vs_circular_geometry':geometry(geo_reference),'eta006_density_vs_circular_geometry':geometry(geo_capped),
                   'endpoint_position_difference_mass_mean':mean(np.hypot(b['x'][mask]-c['x'][mask],b['y'][mask]-c['y'][mask]),w)}
e0=energy(a,'baseline')
result={'passed':True,'time':.825,'eta006_core_screen_passed':bool(core['passed']),'regions':regions,
        'total_energy_relative_to_IC':{'uncapped':(energy(b,'baseline')-e0)/abs(e0),'eta006':(energy(c,'baseline')-e0)/abs(e0)},
        'input_sha256':sha(original),'uncapped_endpoint_sha256':sha(ref),'eta006_endpoint_sha256':sha(end),
        'scope':'Same initial Q1 IDs and .825 physical time, comparing native global steps; eta006 resumed from an independently validated early checkpoint. Circular geometry is kinematic reference, not an exact finite-cell gas solution. No t6.5/t600 stability claim.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
