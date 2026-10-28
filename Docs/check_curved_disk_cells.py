"""True curved-cell boundary/positive-volume/clip fixtures, no fluid RHS."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_circular_sparse_modules_20261019';sys.path[:0]=[str(modules),str(base/'Codex_psi_coupled_modules_20261019')]
from curved_disk_cell import DiskCell
from circular_power_reference import cells as reference_cells
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter();rows={};R=2.5
fixture_root=base/'Codex_circular_sparse_support_fixtures_20261019'
for name,xy,W in [('single',np.array([[.125,-.25]]),np.zeros(1)),('pair',np.array([[-.5,.125],[.5,-.125]]),np.array([1.,-1.])/1024),('weighted_sparse_seed',np.array([(x,y) for y in range(-3,4) for x in range(-3,4)],float)/4,((np.arange(49)*7)%13).astype(float)/65536)]:
 native=json.loads((fixture_root/(name+'.json')).read_text());reference=reference_cells(xy,W,R);curved=[DiskCell.from_native(c,R) for c in native['cells']];boundary=[];quad=[];orders=(12,24)
 for i,(cell,old) in enumerate(zip(curved,reference)):
  V,q,Q=cell.moments();vo,qo,Qo=old.moments();err=max(abs(float(V)-vo),float(np.max(abs(q-qo))),float(np.max(abs(Q-Qo))))
  assert err<2e-11
  nodes=[];perorder=[]
  for order in orders:
   p,w=cell.positive_quadrature(order,(.25,.5,2.));pl=p.astype(np.longdouble);wl=w.astype(np.longdouble)
   value=np.r_[wl.sum(),wl@pl,np.einsum('n,ni,nj->ij',wl,pl,pl).ravel()]
   target=np.r_[V,q,Q.ravel()];difference=float(np.max(abs(value-target)));assert difference<2e-11
   perorder.append(difference);nodes.append(len(p))
  boundary.append(err);quad.append(dict(id=native['cells'][i]['owner'],errors=perorder,nodes=nodes,lines=len(cell.lines),arcs=len(cell.arcs)))
 rows[name]=dict(status='passed',maximum_reference_moment_error=max(boundary),maximum_positive_moment_error=max(max(v['errors']) for v in quad),cells=quad)
 print(json.dumps(dict(fixture=name,status='passed',maximum_reference_moment_error=max(boundary),maximum_positive_moment_error=rows[name]['maximum_positive_moment_error'])),flush=True)
 if name=='single':circle=curved[0]
half=circle.clip_halfplane([1.,0.],0.);quarter=half.clip_halfplane([0.,1.],0.)
for name,cell,fraction in [('half',half,.5),('quarter',quarter,.25)]:
 V,q,Q=cell.moments();p,w=cell.positive_quadrature(24,(.25,.5,2.))
 assert abs(float(V)-fraction*np.pi*R**2)<2e-11 and abs(w.sum()-float(V))<2e-11
 rows[name]=dict(status='passed',area=float(V),target_area=fraction*np.pi*R**2,positive_nodes=len(p),arcs=len(cell.arcs),lines=len(cell.lines))
report=dict(status='passed',scope=__doc__,results=rows,seconds=time.perf_counter()-started,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'curved_disk_cell.py']},limits='Small fixtures only; full 19945-cell disk mesh/IC and closed nonlinear root/dynamics are not covered.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status='passed',seconds=report['seconds'])))
