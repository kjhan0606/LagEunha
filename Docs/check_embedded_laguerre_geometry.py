"""Check moving patch geometry against retained cells and face swept areas."""
import argparse,importlib.util,json,sys
from pathlib import Path
import numpy as np
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--topology',action='store_true');args=parser.parse_args()
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
geom=load('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3]);ids=h['ids']
y=h['state'][:-4].reshape(-1,20);d=h['rates'][:-4].reshape(-1,20)
builder=EmbeddedLaguerreGeometry(allx,ids,clip.clip)
polys,edges=builder(y[:,:2],y[:,2],0.)
reference=geom.polygons(allx,np.zeros(len(allx)),ids,clip.clip)
assert all(np.array_equal(a,b) for a,b in zip(polys,reference))
rate=np.zeros(len(ids))
for i,j,a,b,n,length in edges:
    if j<0:
        assert h['frozen'][i] and np.all(d[i,:3]==0)
        continue
    distance=np.linalg.norm(y[j,:2]-y[i,:2]);point=(a+b)/2
    wn=((y[j,:2]-y[i,:2])@d[j,:2]-(d[j,:2]-d[i,:2])@(point-y[i,:2]))/distance+(d[i,2]-d[j,2])/(2*distance)
    rate[i]+=length*wn;rate[j]-=length*wn
rows=[]
def pairs(edges):return {(i,j) for i,j,a,b,n,length in edges if j>=0 and length>1e-12}
initial_pairs=pairs(edges)
for eps in (1e-4,5e-5,2.5e-5,1.25e-5):
    values=[];changed=set();events=[]
    for sign in (-1,1):
        trial=y+sign*eps*d
        p,e=builder(trial[:,:2],trial[:,2],sign*eps)
        values.append(np.array([geom.area(v) for v in p]))
        if args.topology:
            current=pairs(e);added=current-initial_pairs;removed=initial_pairs-current
            changed.update(k for pair in added|removed for k in pair)
            events.append(dict(sign=sign,added_pairs=[[int(ids[i]),int(ids[j])] for i,j in sorted(added)],
                               removed_pairs=[[int(ids[i]),int(ids[j])] for i,j in sorted(removed)]))
    derivative=(values[1]-values[0])/(2*eps)
    error=float(np.max(abs(derivative-rate)))
    rows.append(dict(epsilon=eps,maximum_area_rate_error=error,relative_L2_error=float(np.linalg.norm(derivative-rate)/np.linalg.norm(rate))))
    if args.topology:
        mask=np.zeros(len(ids),dtype=bool);mask[list(changed)]=True;difference=abs(derivative-rate)
        rows[-1].update(changed_cells=int(mask.sum()),events=events,
            maximum_error_changed=float(np.max(difference[mask],initial=0)),
            maximum_error_unchanged=float(np.max(difference[~mask],initial=0)),
            largest_error_cell=int(ids[np.argmax(difference)]),
            largest_error_cell_changed=bool(mask[np.argmax(difference)]))
        np.savez(f'topology_eps{eps:g}.npz',ids=ids,changed=mask,area_rate=rate,finite_difference=derivative)
    print(json.dumps(rows[-1]),flush=True)
report=dict(scope=__doc__,cells=len(ids),faces=len(edges),initial_polygons_exactly_equal=True,results=rows,
    limits='Directional geometric conservation check only. Fixed external generators; reservoir cells can change shape. No evolved density, hydro, GPU or long-term stability claim.')
Path('report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
assert rows[-1]['maximum_area_rate_error']<1e-8
