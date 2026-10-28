"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
from steady_kepler_ale_reference import boundary_derivative,field_moments,exact_fixed_mass_weights
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
hydro_root=base/'Codex_continuous_initial_rhs_quadratic_split_20261018'
path=hydro_root/'rhs_r0.55_f8.npz';h=np.load(path)
gravity_root=base/'Codex_continuous_quadratic_gravity_20261018'
gr=json.loads((gravity_root/'report.json').read_text())
assert gr['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
g=np.load(gravity_root/'gravity_q24.npz')
assert np.array_equal(g['ids'],h['ids'])
ids=h['ids'];mass=h['mass'];a=h['state'][:-4].reshape(-1,20)
raw=np.column_stack([mass,a[:,3:]])
projection=np.load(base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz')
assert np.array_equal(projection['ids'],ids)
q=projection['moments']
projected_raw=np.column_stack([q[:,:3],q[:,3],q[:,4],q[:,4],q[:,5],q[:,6:12],q[:,12],q[:,13],q[:,13],q[:,14],q[:,15]])
assert np.array_equal(raw,projected_raw)
hr=h['rates'][:-4].reshape(-1,20)
combined=g['combined_rates'][:-4].reshape(-1,20)
active=~h['frozen']
assert np.allclose(combined[active,3:],hr[active,3:]+g['source'][active,1:],rtol=0,atol=1e-20)
actual=np.column_stack([np.zeros(len(ids)),combined[:,3:]])
hydro=np.column_stack([np.zeros(len(ids)),hr[:,3:]])
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3]);assert np.array_equal(a[:,:2],allx[ids])
env=json.loads((base/'Codex_canonical_initial_20261017/run/environment.json').read_text())
assert env['EUNHA_KEPLER_BG_ROTATE']=='1' and env['EUNHA_KEPLER_PROFILE']=='pakmor'
rin=float(env['EUNHA_KEPLER_RIN']);rout=float(env['EUNHA_KEPLER_ROUT']);floor=float(env['EUNHA_KEPLER_FLOOR'])
def field(points):
    return field_moments(points,gravity_coefficient,geom.potential,float(ic['pressure'][0]),float(np.float32(1.666667)),rin,rout,floor)
rc=float(env['EUNHA_KEPLER_RC']);assert rc==.5
field.radii=tuple(sorted(set((rin,rout,rc/2,rc))))
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0];lookup={int(k):i for i,k in enumerate(ids)};edges=[]
for local,i in enumerate(ids):
    for row in faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]:
        other=int(row[1]);j=lookup.get(other,-1)
        if j>=0 and local>=j:continue
        aa=allx[i]+row[6:8];bb=allx[i]+row[8:10];delta=allx[other]-allx[i]
        edges.append((local,j,aa,bb,delta/np.linalg.norm(delta),float(np.linalg.norm(bb-aa))))
geometry_errors=[]
for i,j,aa,bb,n,length in edges:
    geometry_errors.extend([abs(float(n@n)-1),abs(float((bb-aa)@n)),abs(float(np.linalg.norm(bb-aa))-length)])
    if j>=0:
        delta=a[j,:2]-a[i,:2]
        geometry_errors.extend(abs(2*((np.array([aa,bb])-a[i,:2])@delta)-delta@delta-a[i,2]+a[j,2]).tolist())
assert max(geometry_errors)<1e-10

root=base/'Codex_steady_ale_exact_mass_20261018';root.mkdir(exist_ok=False)
rows=[];previous=None
for nq in (16,32):
    z,diagnostic=exact_fixed_mass_weights(a[:,:2],hr[:,:2],edges,field,h['frozen'],nq)
    ref=boundary_derivative(a[:,:2],hr[:,:2],z,edges,field,nq)
    # Reassemble with higher-order quadrature, independent of the sparse solve.
    check=boundary_derivative(a[:,:2],hr[:,:2],z,edges,field,48)
    np.savez(root/f'exact_mass_q{nq}.npz',ids=ids,weight_rates=z,reference=ref,check_reference=check,active=active)
    row=dict(order=nq,**diagnostic,exact_boundary_mass_rate_max=float(abs(check[active,0]).max()),
             actual_minus_exact_weight_rate_max=float(abs(hr[active,2]-z[active]).max()),
             actual_minus_exact_weight_rate_RMS=float(np.sqrt(np.average((hr[active,2]-z[active])**2,weights=mass[active]))),
             frozen_weights_zero=bool(np.all(z[~active]==0)))
    if previous is not None:row['successive_weight_rate_max']=float(abs(z-previous).max())
    previous=z;rows.append(row)
    assert row['linear_residual_max']<1e-12 and row['exact_boundary_mass_rate_max']<1e-12
report=dict(scope='Exact continuous Kepler mass-flux linear solve at unchanged numerical generator velocities; frozen Dirichlet z=0.',results=rows,
 limits='Reference-only weight rates; numerical solver was not changed or reevaluated at these rates. No evolution or production backend validation. Actual-w and exact-z boundary derivatives are both preserved.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
