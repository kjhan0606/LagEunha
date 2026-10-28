"""Compare two dt paths at t=.001 against their own actual-cell references.

Two dt values provide an endpoint comparison, not a convergence-order proof.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
for k in ('coarse','fine','output'):parser.add_argument('--'+k,type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False);reports=[json.loads((p/'report.json').read_text()) for p in (args.coarse,args.fine)];assert all(r['status']=='completed_physical_state_comparison' and r['time']==.001 for r in reports)
a,b=[np.load(p/'fields.npz') for p in (args.coarse,args.fine)];assert np.array_equal(a['ids'],b['ids'])
w=a['actual_mass'];assert np.array_equal(w,b['actual_mass']);n=len(w);ld=np.longdouble
regions={k:a[k]|b[k] for k in ('all','inner','cavity','disk_core','inner_contact','outer_contact','dynamic_buffer','wall')}
def stats(delta,mask):
 count=int(mask.sum())
 if not count:return None
 d=np.asarray(delta[mask],dtype=ld).reshape(count,-1);ww=w[mask].astype(ld)
 return dict(maximum_absolute=float(abs(d).max()),mass_weighted_RMS=float(np.sqrt(np.sum(ww*np.sum(d*d,axis=1))/ww.sum())))
result={}
for name,mask in regions.items():
 row=dict(cells=int(mask.sum()),coarse_cells=int(a[name].sum()),fine_cells=int(b[name].sum()),fields={})
 for key in ('mass','density','centroid','velocity','radial_velocity','azimuthal_velocity','intrinsic_ell','origin_L','R','thermal','specific_thermal','pressure','Psi'):
  av=a['actual_'+key];bv=b['actual_'+key];ae=av-a['reference_'+key];be=bv-b['reference_'+key]
  row['fields'][key]=dict(actual_state_difference=stats(bv-av,mask),physical_error_difference=stats(be-ae,mask))
 c=stats(a['relative_pressure_error'],mask);f=stats(b['relative_pressure_error'],mask);row['relative_pressure_error']=dict(coarse=c,fine=f,fine_to_coarse_RMS_ratio=f['mass_weighted_RMS']/c['mass_weighted_RMS'] if c and c['mass_weighted_RMS'] else None)
 result[name]=row
report=dict(status='completed_full_disk_equal_time_comparison',scope=__doc__,time=.001,coarse_step=.001,fine_step=.0005,fine_steps=2,regions=result,global_values=[r['global_values'] for r in reports],source_sha256={str(p/'report.json'):hashlib.sha256((p/'report.json').read_bytes()).hexdigest() for p in (args.coarse,args.fine)},provider_calls=0,evolution_steps=0,limits='Both paths use evolved states, no reset, on their own moving meshes. Stable-ID actual differences and actual-cell physical-error differences are reported separately. Two dt values do not establish time order or long-term Kepler stability.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],pressure={k:v['relative_pressure_error'] for k,v in result.items()}),indent=2))
