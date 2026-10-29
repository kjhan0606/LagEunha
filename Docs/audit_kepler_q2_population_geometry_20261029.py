"""Test small ring-population perturbations on prescribed Kepler orbits."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
    raise RuntimeError('run on a compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
source=Path(prior['initial'])
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(source)!=prior['initial_sha256']:raise ValueError('unexpected source hash')
p,meta=ordered_global(source)
if meta['t']!=0 or np.any(p['w2']!=0):raise ValueError('invalid initial state')
xy=np.column_stack((p['x'].astype(float)-3,p['y'].astype(float)-3))
radius=np.linalg.norm(xy,axis=1);rr=np.round(radius,10)
disk=(radius>=.5)&(radius<2);rings=np.unique(rr[disk])
if len(rings)!=68:raise ValueError('unexpected ring layout')
new=[xy[~disk]]; changes=[]
for k,r in enumerate(rings):
    ix=np.flatnonzero(disk & (rr==r));angles=np.sort(np.arctan2(xy[ix,1],xy[ix,0]));n=len(ix)
    gaps=np.diff(np.r_[angles,angles[0]+2*np.pi])
    if np.max(abs(gaps/(2*np.pi/n)-1))>1e-4:raise ValueError('nonuniform input ring')
    nn=n+(1 if k%2==0 else -1)
    angles=angles[0]+np.arange(nn)*2*np.pi/nn
    # Keep the actual ring radius and one generator's original angular anchor.
    ractual=float(np.mean(radius[ix]));new.append(ractual*np.column_stack((np.cos(angles),np.sin(angles))))
    changes.append({'radius':float(r),'old_count':n,'new_count':nn})
altered=np.vstack(new)
if len(altered)!=len(xy):raise ValueError('net population changed')

def measure(points,guard,time):
    r=np.linalg.norm(points,axis=1)
    near=points[r<guard];r=r[r<guard]
    selected=np.flatnonzero((r>=.5)&(r<.75));lookup=np.full(len(near),-1,int);lookup[selected]=np.arange(len(selected))
    theta=time*(r*r+.01175**2)**(-.75)
    ct=np.cos(theta);st=np.sin(theta)
    pos=np.column_stack((near[:,0]*ct-near[:,1]*st,near[:,0]*st+near[:,1]*ct))
    area,info,tri,centers=power_areas(pos,np.zeros(len(pos)),selected,lookup,return_geometry=True)
    if info['hidden_active'] or np.any(area<=0):raise ValueError('geometry failed')
    # All polygon vertices must stay well inside the supplied radial support.
    touching=tri[np.any(lookup[tri]>=0,axis=1)]
    if len(touching)==0 or np.max(np.linalg.norm(pos[touching],axis=2))>=guard-.05:
        raise ValueError('active polygons too close to support cutoff')
    return selected,near,r,area

cases={}
for name,points in [('original',xy),('alternating_population',altered)]:
    selected,near,r,a0=measure(points,1.2,0)
    _,_,_,a1=measure(points,1.2,6.5)
    _,_,_,a0check=measure(points,1.4,0)
    _,_,_,a1check=measure(points,1.4,6.5)
    guarderr=max(float(np.max(abs(a0/a0check-1))),float(np.max(abs(a1/a1check-1))))
    if guarderr>1e-8:raise ValueError('guard support disagreement')
    if name=='original':
        oldix=np.flatnonzero((radius<1.2)&(radius>=.5)&(radius<.75))
        err=float(np.max(abs(a0/p['volume'][oldix].astype(float)-1)))
        if err>1e-6:raise ValueError('source area mismatch')
    else:err=None
    change=a0/a1-1 # fixed masses and initial uniform density
    ring_rows=[];rsel=np.round(r[selected],10)
    for rr0 in np.unique(rsel):
        loc=np.flatnonzero(rsel==rr0);ang=np.arctan2(near[selected[loc],1],near[selected[loc],0]);loc=loc[np.argsort(ang)]
        values=a0[loc]/np.mean(a0[loc])-1
        ps=abs(np.fft.rfft(values))**2
        ring_rows.append({'n':len(loc),'initial_area_peak_m':int(np.argmax(ps[1:])+1),
                          'initial_area_rms':float(np.sqrt(np.mean(values**2)))})
    # Uniform initial density implies mass weighting proportional to a0.
    cases[name]={'sites':len(selected),'max_guard_relative_error':guarderr,'source_area_max_relative_error':err,
                 'ballistic_abs_density_change_mass_mean':float(np.sum(a0*abs(change))/np.sum(a0)),
                 'ballistic_abs_density_change_p99':float(np.quantile(abs(change),.99)),
                 'rings':ring_rows}
result={'initial_sha256':sha(source),'ring_counts':changes,'time':6.5,'epsilon':.01175,'cases':cases,
        'script_sha256':sha(Path(__file__)),
        'scope':'Exploratory fixed-position IC modification changing each disk ring count by alternating +1/-1. Inner cells are assessed using independently rebuilt bounded-support Voronoi geometry and prescribed softened circular trajectories, not evolved hydro, force, entropy, conservation or stability. Actual full-domain wall effects and native execution are untested; no generator steering or spatial periodic boundary.'}
with (docs/'KEPLER_Q2_POPULATION_GEOMETRY_2026-10-29.json').open('x') as f:
    json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
