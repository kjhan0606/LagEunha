"""Compare equal-time closures against each endpoint's own current polygons."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--old-reference',type=Path,required=True);parser.add_argument('--new-reference',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
roots=[args.old_reference,args.new_reference];reports=[json.loads((r/'report.json').read_text()) for r in roots];data=[np.load(r/'reference.npz') for r in roots]
assert reports[0]['time']==reports[1]['time']==1e-4 and reports[0]['physical']==reports[1]['physical']
checkpoints=[np.load(report['checkpoint']) for report in reports];initial=[np.load(report['initial']) for report in reports]
for key in ('state','factors','mass','ids','time','weight_rate_seed'):assert np.array_equal(initial[0][key],initial[1][key]),key
for k in ('ids','mass','frozen','initial_radius'):assert np.array_equal(data[0][k],data[1][k]),k
metrics=['pressure_relative_P0','velocity','mass_relative_reference','ell','factor_R','angular_momentum_origin']
regions=['active','inner','inner_transition','contact_band','outer_interior','active_edge','contact_cells_active','active_frozen_boundary_cells_active']
comparison={}
for metric in metrics[:-1]:
 comparison[metric]={}
 for region in regions:
  rows=[report['numerical_minus_q32_selected_else_q16'][metric][region] for report in reports]
  quadrature=[report['q32_minus_q16_selected'][metric][region] for report in reports]
  old,new=[row['mass_weighted_RMS'] for row in rows]
  comparison[metric][region]=dict(old_RMS=old,new_RMS=new,new_over_old=new/old if old else None,old_maximum=rows[0]['maximum'],new_maximum=rows[1]['maximum'],old_maximum_id=rows[0]['maximum_id'],new_maximum_id=rows[1]['maximum_id'],quadrature_selected_RMS=[q.get('mass_weighted_RMS') for q in quadrature])
# Origin angular momentum is independently derived without reusing ell.
mask=~data[0]['frozen'];mass=data[0]['mass'];angular=[]
for d in data:
 ref=d['reference_q32_selected'];raw=d['raw'];delta=(raw[:,10]-raw[:,11])-(ref[:,10]-ref[:,11]);angular.append(dict(active_RMS=float(np.sqrt(np.average(delta[mask]**2,weights=mass[mask]))),active_maximum=float(abs(delta[mask]).max())))
paths=[]
for r,report in zip(roots,reports):paths.extend([r/'report.json',r/'reference.npz',Path(report['checkpoint']),Path(report['initial'])])
report=dict(status='completed_comparison',scope=__doc__,time=1e-4,comparison=comparison,angular_momentum_origin=angular,raw_endpoint_path_difference_max=float(abs(data[1]['raw']-data[0]['raw']).max()),generator_path_difference_max=float(abs(checkpoints[1]['state'][:-4].reshape(-1,20)[:,:2]-checkpoints[0]['state'][:-4].reshape(-1,20)[:,:2]).max()),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),*paths]},limits='Old/new closure on identical initial raw state and horizon; each reference integrates steady Kepler field on its own moving cell. Not time convergence, material orbit error, long-time stability, or a closure adoption verdict. q32 is independent only on saved selected masks.')
args.output.write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],active={m:comparison[m]['active'] for m in metrics[:-1]}),indent=2))
