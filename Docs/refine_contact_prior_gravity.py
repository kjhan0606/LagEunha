"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_contact_prior_rhs_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
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
import quadratic_velocity_reconstruction as quad
from kinetic_factor_rhs import rhs,ALPHA
from quadratic_upwind_face_nodes import QuadraticUpwindNodes
from contact_prior_density import ContactPriorDensity,local_contact_priors
from moment_gravity_source import source as gravity_source
os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0',LAGFORCE_SHARED_FACE_VELOCITY='1',LAGFORCE_CAPTURE_PRESSURE_RATES='1')
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
root=base/'Codex_contact_prior_gravity_refinement_20261018';root.mkdir(exist_ok=False)
polys=[]
for i in ids:
    rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vertices=rows[:,6:8];polys.append(vertices[np.argsort(np.arctan2(vertices[:,1],vertices[:,0]))]+allx[i])
priors=local_contact_priors(polys,raw,edges)
densities=[ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),meta,12) for p,r,meta in zip(polys,raw,priors)]
(root/'priors.json').write_text(json.dumps(priors,indent=2))
def geometry(xx,ww,tt):
    assert tt==0 and np.array_equal(xx,a[:,:2]) and np.array_equal(ww,a[:,2])
    return polys,edges
model.geometry=geometry;model.boundary_velocity=lambda point:np.zeros(2)
L=h['covariance_factors'];flat=h['state'];initial_copy=flat.copy();factor_copy=L.copy()
states=[model.decode(p,r,density=d) for p,r,d in zip(polys,raw,densities)]
for state,factor in zip(states,L):state['R']=factor@factor.T
reconstruction=quad.prepare(states,polys,edges,mass,L)
gen=np.array([model.velocity(state,x) for state,x in zip(states,a[:,:2])]);gen[~active]=0.
gravity=[];moment_defects=[]
for i,(state,density) in enumerate(zip(states,densities)):
    points=density.points;dm=density.weighted_mass
    samples,_=quad.velocity_samples(state,L[i],mass[i],points,ALPHA)
    velocity=samples.mean(axis=1)
    momentum=dm@velocity;D=points.T@(dm[:,None]*velocity)
    Pi=np.einsum('n,nki,nkj->ij',dm/len(ALPHA),samples,samples)
    moment_defects.append(np.r_[momentum-raw[i,7:9],D.ravel()-raw[i,9:13],Pi.ravel()-raw[i,13:17]])
    acc=-gravity_coefficient(np.linalg.norm(points,axis=1))[:,None]*points
    gravity.append(gravity_source(points,dm,velocity,acc))
gravity=np.array(gravity);moment_defects=np.array(moment_defects)
assert abs(moment_defects).max()<1e-11

from steady_kepler_ale_reference import polygon_integral
saved=np.load(base/'Codex_contact_prior_complete_rhs_20261018/rhs.npz')
assert np.array_equal(saved['ids'],ids) and np.array_equal(saved['gravity'],gravity)
rows=[];previous=gravity.copy();np.savez(root/'source_q12.npz',ids=ids,source=gravity,active=active)
for nq in (24,48):
    refined=np.zeros_like(gravity);errors=[]
    for i in np.flatnonzero(active):
        state=states[i];density=densities[i]
        def integrand(points):
            rho=density.value(points)
            velocity=state['v']+(points-state['c'])@state['B'].T+quad.correction(state,points)
            acceleration=-gravity_coefficient(np.linalg.norm(points,axis=1))[:,None]*points
            force=rho[:,None]*acceleration;result=np.zeros((len(points),18))
            result[:,7:9]=force;result[:,9:13]=(points[:,:,None]*force[:,None,:]).reshape(-1,4)
            va=velocity[:,:,None]*force[:,None,:]
            result[:,13:17]=(va+va.transpose(0,2,1)).reshape(-1,4)
            result[:,17]=np.sum(velocity*force,axis=1)
            return result
        for part in density.parts:
            value,error=polygon_integral(part,integrand,(rc/2,rc),nq)
            refined[i]+=value;errors.append(float(error))
        if (i+1)%100==0:
            with (root/'events.jsonl').open('a') as f:f.write(json.dumps(dict(order=nq,cell_index=int(i)))+'\n')
    np.savez(root/f'source_q{nq}.npz',ids=ids,source=refined,active=active)
    delta=refined-previous;bulk=delta[:,7:9]/mass[:,None]
    torque=delta[:,10]-delta[:,11]-(raw[:,1]*delta[:,8]-raw[:,2]*delta[:,7])/mass
    thermal=(delta[:,17]-.5*(delta[:,13]+delta[:,16]))/mass
    rms=lambda a:float(np.sqrt(np.average(np.sum(a[active].reshape(int(active.sum()),-1)**2,axis=1),weights=mass[active])))
    rows.append(dict(order=nq,previous_order=12 if nq==24 else 24,
        active_source_max_difference=float(abs(delta[active]).max()),bulk_rate_difference_RMS=rms(bulk),
        internal_torque_difference_RMS=rms(torque),specific_thermal_rate_difference_RMS=rms(thermal),
        maximum_thermal_source=float(abs(refined[active,17]-.5*(refined[active,13]+refined[active,16])).max()),
        maximum_origin_gravity_torque=float(abs(refined[active,10]-refined[active,11]).max()),
        maximum_adaptive_error_estimate=max(errors)))
    previous=refined
report=dict(results=rows,active_cells=int(active.sum()),
 scope='Fixed fitted density and fixed quadratic mean velocity, independent split-volume gravity integration. Active cells only; frozen source arrays are zero placeholders, not evaluated values.',
 limits='No density/basis refit across gravity quadrature orders; no RHS rerun or time evolution. Adaptive outer integration and circle cuts resolve the spline branches. Potential closure remains unresolved.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
