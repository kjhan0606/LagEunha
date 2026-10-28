"""Equal-time endpoints of two independently driven one-step paths.

No RHS or new evolution; factor gauges are compared through L L.T.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_compiled_midpoint_once_20261018'
status=json.loads((root/'summary.json').read_text());assert status['status']=='completed'
sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from contact_prior_density import polygon_area
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py');clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
paths=[base/'Codex_contact_prior_evolution_dt1e4_20261018/step0001.npz',root/'accepted.npz']
checkpoints=[np.load(path) for path in paths];initial=np.load(root/'initial.npz');mass=initial['mass'];ids=initial['ids'];frozen=initial['frozen']
assert all(float(a['time'])==1e-4 and np.array_equal(a['ids'],ids) and np.array_equal(a['mass'],mass) for a in checkpoints)
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3]);builder=EmbeddedLaguerreGeometry(allx,ids,clip.clip)
fields=[];diagnostics=[]
for saved in checkpoints:
    a=saved['state'][:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]]);q=raw[:,1:3];p=raw[:,7:9]
    assert np.array_equal(a[frozen],initial['state'][:-4].reshape(-1,20)[frozen])
    S=raw[:,3:7].reshape(-1,2,2)-q[:,:,None]*q[:,None,:]/mass[:,None,None]
    C=raw[:,9:13].reshape(-1,2,2)-q[:,:,None]*p[:,None,:]/mass[:,None,None]
    B=np.linalg.solve(S,C).transpose(0,2,1);L=saved['factors'];R=L@L.transpose(0,2,1)
    Pi=raw[:,13:17].reshape(-1,2,2);thermal=raw[:,17]-.5*np.trace(Pi,axis1=1,axis2=2)
    polygons,_=builder(a[:,:2],a[:,2],float(saved['time']));area=np.array([polygon_area(poly) for poly in polygons])
    pressure=(float(np.float32(1.666667))-1)*thermal/area
    assert all(np.isfinite(v).all() for v in (raw,L,area,pressure)) and np.all(area>0)
    fields.append(dict(raw=raw,generator=a[:,:2],weight=a[:,2],centroid=q/mass[:,None],velocity=p/mass[:,None],ell=C[:,0,1]-C[:,1,0],R=R,thermal=thermal,specific_thermal=thermal/mass,pressure=pressure))
    rebuilt=p[:,:,None]*p[:,None,:]/mass[:,None,None]+B@S@B.transpose(0,2,1)+R
    diagnostics.append(dict(minimum_thermal=float(thermal.min()),minimum_R_eigenvalue=float(np.linalg.eigvalsh(R).min()),maximum_Pi_factor_mismatch=float(abs(Pi-rebuilt).max()),minimum_pressure=float(pressure.min())))
regions={'active':~frozen,'inner':(~frozen)&(np.linalg.norm(initial['state'][:-4].reshape(-1,20)[:,:2],axis=1)<.4),'frozen':frozen}
comparison={}
for name,mask in regions.items():
    comparison[name]=dict(cells=int(mask.sum()),fields={})
    for key in fields[0]:
        delta=(fields[1][key]-fields[0][key])[mask].reshape(int(mask.sum()),-1)
        comparison[name]['fields'][key]=dict(maximum_absolute_difference=float(abs(delta).max()),mass_weighted_RMS=float(np.sqrt(np.average((delta*delta).sum(1),weights=mass[mask]))))
report=dict(scope=__doc__,time=1e-4,paths=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],diagnostics=diagnostics,comparison=comparison,
    pressure='(gamma-1)*(Ehyd-trPi/2)/actual moving polygon area; no density fit or RHS rerun.',
    factor_gauge='Compared L L.T, not factor entries; independent paths may use opposite QR column signs.',
    frozen_pressure_note='Frozen raw and factors are held exactly, but neighboring motion may change frozen polygon area and pressure. A pressure difference alone is not a frozen-state violation.',
    limits='Endpoint difference between independent Python and CUDA-driven midpoint paths; not a common-stage shadow, time-order result, full GPU geometry/MPI or long Kepler validation.')
(root/'endpoint_comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
