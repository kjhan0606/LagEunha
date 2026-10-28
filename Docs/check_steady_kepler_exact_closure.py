"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
from steady_kepler_ale_reference import boundary_derivative,field_moments,polygon_integral,exact_face_rates
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

root=base/'Codex_steady_ale_exact_closure_20261018';root.mkdir(exist_ok=False)
pressure=float(ic['pressure'][0]);polys=[]
for i in ids:
    rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vertices=rows[:,6:8];polys.append(vertices[np.argsort(np.arctan2(vertices[:,1],vertices[:,0]))]+allx[i])
def integrand(points):
    density,velocity=field(points);rho=density[:,0]
    acceleration=-gravity_coefficient(np.linalg.norm(points,axis=1))[:,None]*points
    force=rho[:,None]*acceleration;gravity=np.zeros_like(density)
    gravity[:,7:9]=force
    gravity[:,9:13]=(points[:,:,None]*force[:,None,:]).reshape(-1,4)
    VA=velocity[:,:,None]*force[:,None,:]
    gravity[:,13:17]=(VA+VA.transpose(0,2,1)).reshape(-1,4)
    gravity[:,17]=np.sum(velocity*force,axis=1)
    gravity[:,18]=-gravity[:,17]
    return np.column_stack([density,gravity,np.ones(len(points))])
volume_results=[];report_rows=[]
for nq in (8,16):
    values=[];errors=[]
    for index,poly in enumerate(polys):
        value,error=polygon_integral(poly,integrand,field.radii,nq)
        values.append(value);errors.append(float(error))
        if (index+1)%100==0:
            with (root/'events.jsonl').open('a') as f:f.write(json.dumps(dict(order=nq,cells=index+1))+'\n')
    values=np.array(values);moments=values[:,:19];gravity=values[:,19:38];areas=values[:,38]
    volume=np.zeros_like(moments)
    volume[:,1:3]=moments[:,7:9]
    D=moments[:,9:13].reshape(-1,2,2)
    volume[:,3:7]=(D+D.transpose(0,2,1)).reshape(-1,4)
    volume[:,9:13]=moments[:,13:17]
    volume[:,9]+=pressure*areas;volume[:,12]+=pressure*areas
    np.savez(root/f'volume_q{nq}.npz',ids=ids,moments=moments,gravity=gravity,volume=volume,areas=areas)
    checks=[]
    for label,z in [('actual',hr[:,2]),('exact_mass',np.load(base/'Codex_steady_ale_exact_mass_20261018/exact_mass_q32.npz')['weight_rates'])]:
        adv,pface=exact_face_rates(a[:,:2],hr[:,:2],z,edges,field,pressure,32)
        reference=boundary_derivative(a[:,:2],hr[:,:2],z,edges,field,32)
        total=adv+pface+volume+gravity;defect=total-reference
        np.savez(root/f'closure_{label}_q{nq}.npz',ids=ids,advective=adv,pressure=pface,volume=volume,gravity=gravity,reference=reference,total=total,defect=defect)
        checks.append(dict(face_motion=label,active_max_abs=float(abs(defect[active]).max()),all_cells_max_abs=float(abs(defect).max()),
            active_column_max_abs=np.max(abs(defect[active]),axis=0).tolist(),
            shared_boundary_global_sum_max=float(abs(reference.sum(axis=0)).max())))
    c=moments[:,1:3]/moments[:,0,None]
    torque=gravity[:,10]-gravity[:,11]
    center_force=c[:,0]*gravity[:,8]-c[:,1]*gravity[:,7]
    internal=torque-center_force
    row=dict(volume_order=nq,maximum_quad_vec_error_estimate=max(errors),checks=checks,
        initial_raw_projection_max_difference=float(abs(moments[:,:18]-raw).max()),
        central_gravity_origin_torque_max=float(abs(torque).max()),
        gravity_internal_identity_max=float(abs(internal+center_force).max()),
        circular_gravity_work_max=float(abs(gravity[:,17]).max()),
        model_minus_exact_gravity_max=float(abs(g['source']-gravity[:,:18]).max()))
    if volume_results:row['successive_volume_max_abs']=float(abs(values-volume_results[-1]).max())
    volume_results.append(values);report_rows.append(row)
report=dict(scope='Independent exact-field ALE face flux plus discontinuity-split volume/gravity closure; constant-pressure Pi force source=0 explicitly.',results=report_rows,
 limits='Reference identity, not model improvement or production/backend/long-time validation. Frozen cells are not forced to zero in reference. Numerical potential evolution remains unavailable. Exact-mass and actual face motions share volume integrals.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
