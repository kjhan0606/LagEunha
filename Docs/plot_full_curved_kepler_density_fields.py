"""Render cell mean density and actual-cell errors; no particle scatter.

Arc sampling here is solely for image rendering, never solver geometry.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,Normalize
parser=argparse.ArgumentParser(description=__doc__)
for k in ('comparison','geometry','output'):parser.add_argument('--'+k,type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False);base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base/'Codex_full_disk_arc_roots_modules_v7_20261019'))
from curved_disk_cell import DiskCell
fields=np.load(args.comparison/'fields.npz');report=json.loads((args.comparison/'report.json').read_text());assert report['status']=='completed_physical_state_comparison';mesh=np.load(args.geometry/'mesh.npz');native=json.loads((args.geometry/'full_exact_supports.json').read_text());assert np.array_equal(fields['ids'],mesh['ids']);polygons=[]
centers=mesh['moments'][:,1:3]/mesh['moments'][:,0,None]
for i,record in enumerate(native['cells']):
 assert record['owner']==int(fields['ids'][i]);cell=DiskCell.from_native(record,2.5);points=[]
 for line in cell.lines:points.extend([list(map(float,line.A)),list(map(float,line.B))])
 for arc in cell.arcs:
  theta=np.linspace(arc.lo,arc.hi,max(4,int(np.ceil((arc.hi-arc.lo)*400)))+1);points.extend((arc.radius*np.column_stack([np.cos(theta),np.sin(theta)])).tolist())
 xy=np.unique(np.asarray(points),axis=0);delta=xy-centers[i];polygons.append(xy[np.argsort(np.arctan2(delta[:,1],delta[:,0]))])
rho=fields['actual_density'];drho=fields['relative_density_reference_error'];dp=fields['relative_pressure_error'];assert np.all(rho>0) and np.isfinite(drho).all() and np.isfinite(dp).all()
fig,axes=plt.subplots(1,3,figsize=(14,4.5),constrained_layout=True)
for ax,data,title,cmap,norm in [(axes[0],rho,'Cell mean density (m / V)','viridis',LogNorm(vmin=min(.01,float(rho.min())),vmax=max(1.,float(rho.max())))),(axes[1],drho,'Density relative error on actual cells','RdBu_r',Normalize(vmin=-float(abs(drho).max() or 1),vmax=float(abs(drho).max() or 1))),(axes[2],dp,'Pressure relative error','RdBu_r',Normalize(vmin=-float(abs(dp).max() or 1),vmax=float(abs(dp).max() or 1)))]:
 collection=PolyCollection(polygons,array=data,cmap=cmap,norm=norm,edgecolors='none',rasterized=True);ax.add_collection(collection);ax.set(xlim=(-2.55,2.55),ylim=(-2.55,2.55),xlabel='x',ylabel='y',title=title);ax.set_aspect('equal');fig.colorbar(collection,ax=ax,shrink=.8)
fig.suptitle(f"Kepler full disk, t = {report['time']:g}; actual moving-cell reference")
path=args.output/'density_fields.png';fig.savefig(path,dpi=180);plt.close(fig)
manifest=dict(scope=__doc__,time=report['time'],image=str(path),cells=len(polygons),reference='Steady field integrated on each actual moving cell; no reference injection into state',density='Cell mean m/V, not a scatter plot and not the full intra-cell fitted profile',rendering='True arcs sampled for display only; exact/curved solver geometry and data are unchanged',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),args.comparison/'report.json',args.geometry/'report.json']})
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
