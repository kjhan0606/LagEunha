"""Physical errors at actual saved times; safe while new runs are writing."""
import argparse,hashlib,json,re,struct,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--original-baseline-only',action='store_true')
parser.add_argument('--case-only',help='Analyze only this manifest case, without historical folders')
parser.add_argument('--reference-snapshot',type=Path,help='Analyze exactly one existing external snapshot with the same physical definitions')
parser.add_argument('--root',type=Path,default=b/'Codex_smooth_qr_physical_to_t10_20261018')
args=parser.parse_args();root=args.root
if args.reference_snapshot and (args.case_only or args.original_baseline_only):
 parser.error('--reference-snapshot cannot be combined with a case selector')
manifest=json.loads((root/'manifest.json').read_text())
ic,tic=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
assert abs(tic)<1e-12
x0=np.column_stack((ic['x']-3,ic['y']-3));r0=np.linalg.norm(x0,axis=1)
theta0=np.arctan2(x0[:,1],x0[:,0]);om0=np.sqrt(gravity_coefficient(r0));mass=ic['mass'];gamma=float(np.float32(1.666667))
# Same masks/radius calculation as analyze_canonical_from_t0_t12.py.
masks={'inner_background':(r0<.5)&(ic['den']<.5),
 'inner_disk_edge':(r0<.6)&(ic['den']>.5),
 'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5),
 'inner_shear_band':(r0>=.3)&(r0<.5)}
folders={'legacy_history_to_t4':b/'Codex_coupled_onset_to_t4_20261018/reference',
         'legacy_continue':b/'Codex_smooth_qr_physical_to_t10_20261018/legacy_continue',
         **{name:root/name for name in manifest['cases']}}
if args.original_baseline_only:folders={'original_baseline':b/'Codex_canonical192_accepted_20261017/omp8_t20'}
if args.case_only:
 assert args.case_only in manifest['cases']
 folders={args.case_only:root/args.case_only}
if args.reference_snapshot:
 args.reference_snapshot=args.reference_snapshot.resolve(strict=True)
 key='reference_'+hashlib.sha256(str(args.reference_snapshot.parent).encode()).hexdigest()[:16]
 folders={key:args.reference_snapshot.parent}
summary={}
for name,folder in folders.items():
 dest=root/'physical_frames'/name;dest.mkdir(parents=True,exist_ok=True)
 rows=[];pending=[]
 for p in ([args.reference_snapshot] if args.reference_snapshot else sorted(folder.glob('khout.*.dat'))):
  with p.open('rb') as f:
   head=f.read(16)
   if len(head)!=16:pending.append(p.name);continue
   magic,n,size,stress=struct.unpack('<4i',head)
   assert magic==-1 and n==len(ic) and size==624
   if p.stat().st_size!=16+n*size+52:pending.append(p.name);continue
  target=dest/(p.name+'.json')
  if target.exists():rows.append(json.loads(target.read_text()));continue
  a,t=read(p);assert np.array_equal(a['id'] & np.uint64(72057594037927935),ic['id'] & np.uint64(72057594037927935))
  assert all(np.isfinite(a[k]).all() for k in a.dtype.names)
  assert np.array_equal(a['mass'],mass) and np.all(a['w2']==0)
  assert np.all(a['pressure']>0) and np.all(a['volume']>0)
  x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));r=np.linalg.norm(x,axis=1)
  vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
  phase=np.angle(np.exp(1j*(np.arctan2(x[:,1],x[:,0])-theta0-om0*t)))
  co=np.cos(om0*t);si=np.sin(om0*t);xe=np.column_stack((co*x0[:,0]-si*x0[:,1],si*x0[:,0]+co*x0[:,1]))
  logK=np.log(a['pressure']/ic['pressure'])+gamma*np.log(a['volume']/ic['volume'])
  row={'snapshot':str(p),'sha256':sha(p),'time':t,'mass_unchanged':True,'weights_zero':True,'regions':{}}
  for label,mask in masks.items():
   def rms(z):return float(np.sqrt(np.average(z[mask]**2,weights=mass[mask])))
   row['regions'][label]={'cells':int(mask.sum()),'pressure_error_mass_RMS':rms(a['pressure']/ic['pressure']-1),
    'phase_error_mass_RMS_deg':rms(phase)*180/np.pi,'radial_displacement_mass_RMS':rms(r-r0),
    'circular_velocity_error_mass_RMS':rms(np.linalg.norm(v-vc,axis=1)),
    'exact_orbit_position_error_mass_RMS':rms(np.linalg.norm(x-xe,axis=1)),
    'log_entropy_ratio_mass_RMS':rms(logK)}
  row['scope']='Unperturbed physical trajectory relative to original IC at actual snapshot time. Phase is wrapped modulo2pi; position/velocity errors provide additional checks. Not a stability pass or inferred exact-time cross-run comparison.'
  with target.open('x') as f:json.dump(row,f,indent=2)
  rows.append(row)
 summary[name]={'complete_frames':len(rows),'incomplete_skipped':pending,'latest':rows[-1] if rows else None}
print(json.dumps(summary,indent=2))
