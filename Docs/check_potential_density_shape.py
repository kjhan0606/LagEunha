"""Sample physical density distortion of the added-potential fit."""
import importlib.util,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
from potential_moment_density import PotentialDensity
import moment_ale_evolve as model
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
data=np.load(base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz');ids=data['ids'];q=data['moments'];cross=data['crosses_inner_edge']
co=np.load(base/'Codex_potential_moment_density_20261018/coefficients.npz')['coefficients'];assert np.array_equal(co[:,0],ids)
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');x=np.column_stack([ic['x']-3,ic['y']-3])
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']='12';results=[]
for local,(i,a,c) in enumerate(zip(ids,q,co)):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]];vv=ff[:,6:8]
    poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    old=model.Density(poly,a[0],a[1:3],np.array([[a[3],a[4]],[a[4],a[5]]]))
    new=PotentialDensity.__new__(PotentialDensity);new.origin=c[1:3];new.scale=c[3];new.beta=c[4:10];new.potential_scale=c[10];new.theta=c[11:17];new.logZ=c[17];new.mass=a[0];new.potential=geom.potential
    pts=np.concatenate([poly,.5*(poly+np.roll(poly,-1,axis=0))]);exact=np.where(np.linalg.norm(pts,axis=1)<.5,.01,1.)
    ro=old.value(pts);rn=new.value(pts)
    results.append(dict(id=int(i),crossing=bool(cross[local]),old_min=float(ro.min()),new_min=float(rn.min()),old_max=float(ro.max()),new_max=float(rn.max()),old_relative_error_rms=float(np.sqrt(np.mean((ro/exact-1)**2))),new_relative_error_rms=float(np.sqrt(np.mean((rn/exact-1)**2)))))
summary={}
for key in (False,True):
    rr=[r for r in results if r['crossing']==key]
    summary[str(key)]=dict(cells=len(rr),old_min=min(r['old_min'] for r in rr),new_min=min(r['new_min'] for r in rr),old_max=max(r['old_max'] for r in rr),new_max=max(r['new_max'] for r in rr),old_cells_below_half_floor=sum(r['old_min']<.005 for r in rr),new_cells_below_half_floor=sum(r['new_min']<.005 for r in rr))
print(json.dumps(dict(scope=__doc__,summary=summary,worst_new_min=sorted(results,key=lambda r:r['new_min'])[:5],limits='Vertices and face midpoints only, not extrema proof. Sharp contact cannot be inferred exactly from a finite list of moments. No state repair applied.'),indent=2))
