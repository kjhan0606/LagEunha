"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_density_trace_audit_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from steady_kepler_ale_reference import boundary_derivative,field_moments,interface_cuts
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

import os,time
import moment_ale_evolve as model
import fixedmass_hllc
from kinetic_factor_rhs import rhs
from quadratic_upwind_face_nodes import QuadraticUpwindNodes
os.environ.update(LAGFORCE_DENSITY_QUADRATURE_ORDER='12',LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_SHARED_FACE_VELOCITY='1')
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
polys=[]
for i in ids:
    rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vertices=rows[:,6:8];polys.append(vertices[np.argsort(np.arctan2(vertices[:,1],vertices[:,0]))]+allx[i])
def geometry(xx,ww,tt):
    assert np.array_equal(xx,a[:,:2]) and np.array_equal(ww,a[:,2]) and tt==0
    return polys,edges
model.geometry=geometry;model.boundary_velocity=lambda point:np.zeros(2)
def face_density(point,i,j):
    radius=np.linalg.norm(point)
    value=1. if rin<=radius<=rout else floor
    return value,value
root=base/'Codex_density_trace_audit_20261018';root.mkdir(exist_ok=False)
L=h['covariance_factors'];flat=h['state'];radius=np.linalg.norm(a[:,:2],axis=1)
initial_copy=flat.copy();factor_copy=L.copy()
results=[];baseline_z=None
cross=lambda x,y:x[:,0]*y[:,1]-x[:,1]*y[:,0]
for case in ('baseline','trace_fixed_z','trace_solved_z'):
    seed=hr[:,2].copy() if baseline_z is None else baseline_z.copy()
    history=[];start=time.monotonic()
    for iteration in range(10):
        provider=QuadraticUpwindNodes(a[:,:2],mass,L,seed,8,model.boundary_velocity,extra_cuts=lambda aa,bb:interface_cuts(aa,bb,field.radii))
        d,diag=rhs(0,flat,mass,covariance_factors=L,frozen_mask=~active,face_order=8,face_node_provider=provider,
            audit_face_density=None if case=='baseline' else face_density,
            audit_weight_rates=baseline_z if case=='trace_fixed_z' else None)
        numeric=d[:-4].reshape(-1,20);newz=numeric[:,2].copy();change=float(abs(newz-seed).max());seed=newz
        event=dict(case=case,iteration=iteration,weight_rate_change=change,mass_flux=float(diag['mass_flux']))
        history.append(event)
        with (root/'events.jsonl').open('a') as stream:stream.write(json.dumps(event)+'\n')
        if case=='trace_fixed_z' or change<1e-11:break
    else:raise RuntimeError('density trace mass/split outer iteration failed')
    assert np.array_equal(numeric[:,:2],hr[:,:2])
    assert np.array_equal(flat,initial_copy) and np.array_equal(L,factor_copy)
    if case!='trace_fixed_z':assert diag['mass_flux']<1e-11
    if case=='baseline':baseline_z=seed.copy()
    actual_raw=np.column_stack([diag['mass_rates'],numeric[:,3:]])
    actual_raw[active]+=g['source'][active]
    ref=boundary_derivative(a[:,:2],numeric[:,:2],seed,edges,field,32)[:,:18]
    delta=actual_raw-ref
    v=raw[:,7:9]/mass[:,None];c=raw[:,1:3]/mass[:,None]
    U=raw[:,17]-.5*(raw[:,13]+raw[:,16]);u=U/mass
    values=dict(bulk_mean_velocity_rate=(delta[:,7:9]-v*delta[:,0,None])/mass[:,None],
        centroid_velocity=(delta[:,1:3]-c*delta[:,0,None])/mass[:,None],
        specific_thermal_rate=(delta[:,17]-.5*(delta[:,13]+delta[:,16])-u*delta[:,0])/mass,
        internal_angular_rate=delta[:,10]-delta[:,11]-(cross(delta[:,1:3],raw[:,7:9])+cross(raw[:,1:3],delta[:,7:9]))/mass+cross(raw[:,1:3],raw[:,7:9])*delta[:,0]/mass**2,
        pressure_rate=diag['pressure_rates'])
    metrics={};largest={}
    for name,mask in [('active',active),('inner_r_lt_040',active&(radius<.4)),('r040_048',active&(radius>=.4)&(radius<.48)),('contact_r048_055',active&(radius>=.48))]:
        metrics[name]={}
        for key,val in values.items():
            q=val[mask].reshape(int(mask.sum()),-1)
            metrics[name][key]=float(np.sqrt(np.average(np.sum(q*q,axis=1),weights=mass[mask])))
    for key,val in values.items():
        magnitude=np.linalg.norm(val.reshape(len(ids),-1),axis=1);top=np.flatnonzero(active)[np.argsort(magnitude[active])[-5:][::-1]]
        largest[key]=[dict(id=int(ids[i]),radius=float(radius[i]),norm=float(magnitude[i])) for i in top]
    np.savez(root/f'{case}.npz',ids=ids,rates=d,mass_rates=diag['mass_rates'],reference=ref,raw_derivative=actual_raw,**values)
    results.append(dict(case=case,history=history,metrics=metrics,largest=largest,seconds=time.monotonic()-start,
        baseline_rate_max_difference=float(abs(d-h['rates']).max()) if case=='baseline' else None))
report=dict(results=results,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
 scope='Counterfactual face-density-only audit. Decoded volume density, raw state, quadratic velocity and pressure reconstruction, generators/factors and gravity unchanged. Pressure F/H/Pi source recalculated. Each case uses its own face motion exact reference.',
 limits='Analytic density trace is audit-only, not a moment-consistent closure or evolution method. C retains kinetic residual samples, unlike monokinetic exact-z. No time integration, production/backend or long-term validation.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
