"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
from steady_kepler_ale_reference import boundary_derivative,field_moments
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
root=base/'Codex_steady_ale_initial_balance_20261018';root.mkdir(exist_ok=False)
radius=np.linalg.norm(a[:,:2],axis=1)
masks={'active':active,'inner_r_lt_040':active&(radius<.4),'contact_r048_055':active&(radius>=.48)}
def summarize(values,mask):
    b=values[mask];return dict(max_abs=float(abs(b).max()),mass_weighted_RMS=float(np.sqrt(np.average(np.sum(b.reshape(len(b),-1)**2,axis=1),weights=mass[mask]))))
def ell_rate(rates):
    q=raw[:,1:3];p=raw[:,7:9];dq=rates[:,1:3];dp=rates[:,7:9]
    cross=lambda x,y:x[:,0]*y[:,1]-x[:,1]*y[:,0]
    return rates[:,10]-rates[:,11]-(cross(dq,p)+cross(q,dp))/mass+cross(q,p)*rates[:,0]/mass**2
previous=None;rows=[]
for nq in (8,16,32):
    exact=boundary_derivative(a[:,:2],hr[:,:2],hr[:,2],edges,field,nq)
    assert np.isfinite(exact).all()
    np.savez(root/f'reference_q{nq}.npz',ids=ids,reference=exact,numerical=actual,
             hydro=hydro,gravity=g['source'],mass=mass,active=active)
    regions={}
    for name,mask in masks.items():
        regions[name]=dict(cells=int(mask.sum()),
            momentum_per_mass=summarize((actual[:,7:9]-exact[:,7:9])/mass[:,None],mask),
            D_per_mass=summarize((actual[:,9:13]-exact[:,9:13])/mass[:,None],mask),
            Pi_per_mass=summarize((actual[:,13:17]-exact[:,13:17])/mass[:,None],mask),
            energy_per_mass=summarize((actual[:,17]-exact[:,17])/mass,mask),
            exact_mass_rate=summarize(exact[:,0],mask),
            internal_angular_rate_defect=summarize(ell_rate(actual)-ell_rate(exact),mask),
            hydro_internal_angular_rate=summarize(ell_rate(hydro),mask),
            gravity_internal_torque=summarize(ell_rate(g['source']),mask),
            exact_internal_angular_rate=summarize(ell_rate(exact),mask))
    row=dict(order=nq,regions=regions)
    if previous is not None:row['reference_successive_max_abs']=float(abs(exact-previous).max())
    previous=exact;rows.append(row)
report=dict(scope=__doc__,hydro_source=str(path),hydro_sha256=gr['sha256'],results=rows,
    initial_raw_projection_max_difference=float(abs(raw-projected_raw).max()),
    circular_split_radii=field.radii,maximum_geometry_defect=max(geometry_errors),
    potential_transport_status='unavailable',
    limits='Initial continuous projection only. Frozen-exterior j<0 normal speed is 0 by this experiment contract. Fixed actual computed face speed; exact-field mass rate need not vanish. Numerical dm is represented as 0 under the solved net-mass constraint; true solver residual was reported separately. Gravity uses fitted density and reconstructed velocity; independent exact-discontinuous-volume gravity and exact-flux closure check are pending. Potential boundary derivative is saved but no numerical potential transport exists to compare. This is not an evolved-state exact solution or a long-time/backend pass.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
