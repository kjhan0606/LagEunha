"""Associate equal-time step-size differences with endpoint neighbor changes."""
import argparse,importlib.util,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--active-radius',type=float,choices=(.55,.60),default=.55)
parser.add_argument('--extended-time',type=float)
args=parser.parse_args()
sys.path[:0]=[str(base/'Codex_continuous_quadratic_evolution_modules_20261018'),str(base)]
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
paths=[base/'Codex_continuous_quadratic_evolution_dt1e4_20261018/step0001.npz',
       base/'Codex_continuous_quadratic_evolution_dt5e5_20261018/step0002.npz']
if args.active_radius==.60:
    paths=[base/'Codex_quadratic_r060_dt1e4_20261018/step0001.npz',
           base/'Codex_quadratic_r060_dt5e5_20261018/step0002.npz']
if args.extended_time is not None:
    if args.active_radius != .60 or not np.isfinite(args.extended_time):
        parser.error('extended time requires active radius .60 and a finite time')
    paths=[]
    for tag,dt in [('dt1e4',1e-4),('dt5e5',5e-5)]:
        step=round((args.extended_time-1e-4)/dt)
        if step<1 or abs(1e-4+step*dt-args.extended_time)>1e-12:
            parser.error('requested time must be a shared accepted-step boundary')
        paths.append(base/f'Codex_quadratic_extended_{tag}_20261018'/f'step{step:04d}.npz')
saved=[np.load(p) for p in paths];ids=saved[0]['ids'];mass=saved[0]['mass']
assert np.array_equal(ids,saved[1]['ids']) and np.array_equal(mass,saved[1]['mass'])
assert float(saved[0]['time'])==float(saved[1]['time'])
if args.extended_time is not None:
    assert abs(float(saved[0]['time'])-args.extended_time)<1e-12
initial=np.load(paths[0].parent/'initial.npz')['state'][:-4].reshape(-1,20)
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz');active=~h['frozen']
# Keep the original cells in both boundary experiments.
assert np.all(np.linalg.norm(initial[active,:2],axis=1)<.55)
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
def geometry(a):
    polys,edges=builder(a[:,:2],a[:,2],0.)
    areas=np.array([.5*np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)) for p in polys])
    assert np.all(areas>0)
    return {(i,j) for i,j,aa,bb,n,length in edges if j>=0 and length>1e-12},areas
old,_=geometry(initial);states=[s['state'][:-4].reshape(-1,20) for s in saved]
geometries=[geometry(a) for a in states]
graphs=[g[0] for g in geometries];changed=set();events=[]
for graph in graphs:
    added=graph-old;removed=old-graph
    changed.update(k for pair in added|removed for k in pair)
    events.append(dict(added=[[int(ids[i]),int(ids[j])] for i,j in sorted(added)],
                       removed=[[int(ids[i]),int(ids[j])] for i,j in sorted(removed)]))
mask=np.zeros(len(ids),dtype=bool);mask[list(changed)]=True
dv=(states[1][:,9:11]-states[0][:,9:11])/mass[:,None]
square=np.sum(dv*dv,axis=1);total=float(np.sum(mass[active]*square[active]));rows=[]
for name,selected in [('changed',active&mask),('unchanged',active&~mask)]:
    contribution=float(np.sum(mass[selected]*square[selected]))
    rows.append(dict(region=name,cells=int(selected.sum()),velocity_mass_RMS=float(np.sqrt(np.average(square[selected],weights=mass[selected]))) if np.any(selected) else None,
        fraction_of_active_squared_error=contribution/total if total else None))
fields=[]
for a,(_,area) in zip(states,geometries):
    C=a[:,11:15].reshape(-1,2,2)-a[:,3:5,None]*a[:,None,9:11]/mass[:,None,None]
    U=a[:,19]-.5*(a[:,15]+a[:,18])
    assert np.all(U>0)
    fields.append(dict(internal_angular_momentum=C[:,0,1]-C[:,1,0],
                       cell_average_pressure=(float(np.float32(1.666667))-1)*U/area))
localization={}
radius=np.linalg.norm(h['state'][:-4].reshape(-1,20)[:,:2],axis=1)
regions={'changed':active&mask,'unchanged':active&~mask,
         'inner_r_lt_040':active&(radius<.4),
         'r040_048':active&(radius>=.4)&(radius<.48),
         'r048_052':active&(radius>=.48)&(radius<.52),
         'r052_055':active&(radius>=.52)}
for key in fields[0]:
    diff=fields[1][key]-fields[0][key];sq=mass*diff**2;total=sq[active].sum()
    entries=[]
    for name,selected in regions.items():
        entries.append(dict(region=name,cells=int(selected.sum()),
            mass_weighted_RMS=float(np.sqrt(sq[selected].sum()/mass[selected].sum())) if np.any(selected) else None,
            fraction_of_squared_difference=float(sq[selected].sum()/total) if total else None))
    top=np.where(active)[0][np.argsort(sq[active])[-10:][::-1]]
    localization[key]=dict(regions=entries,largest_contributors=[dict(id=int(ids[k]),radius=float(radius[k]),
        difference=float(diff[k]),changed=bool(mask[k]),fraction=float(sq[k]/total) if total else None) for k in top])
print(json.dumps(dict(scope=__doc__,active_radius=args.active_radius,comparison_region='original active cells r<.55',time=float(saved[0]['time']),endpoint_graphs_equal=graphs[0]==graphs[1],
    events=events,results=rows,field_localization=localization,limits='Pressure is cell-average (gamma-1)U/area, not face pressure. Differences are between time steps, not errors against an exact solution. Endpoint neighbor sets only: does not detect changes that reverse within a step, nor establish causality or temporal convergence. Face-length threshold 1e-12.'),indent=2))
