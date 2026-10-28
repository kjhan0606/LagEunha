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
root=base/'Codex_contact_prior_complete_rhs_20261018';root.mkdir(exist_ok=False)
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
class ContactNodes(QuadraticUpwindNodes):
    def __call__(self,edge_id,edge,states,gen):
        i,j,aa,bb,*_=edge
        def cuts(a,b):
            values=[0.,1.]
            for k in (i,j):
                if k<0 or priors[k]['flat']:continue
                meta=priors[k];n=np.asarray(meta['normal']);den=float((b-a)@n)
                if den!=0:
                    t=(meta['offset']-a@n)/den
                    if 0<t<1:values.append(float(t))
            return np.unique(values)
        self.extra_cuts=cuts
        return super().__call__(edge_id,edge,states,gen)
z=hr[:,2].copy();history=[];start=time.monotonic()
for iteration in range(10):
    nodes=ContactNodes(a[:,:2],mass,L,z,8,model.boundary_velocity)
    d,diag=rhs(0,flat,mass,covariance_factors=L,frozen_mask=~active,face_order=8,face_node_provider=nodes,density_provider=lambda i,p,r:densities[i])
    rates=d[:-4].reshape(-1,20);newz=rates[:,2].copy();change=float(abs(newz-z).max());z=newz
    event=dict(iteration=iteration,weight_rate_change=change,mass_flux=float(diag['mass_flux']))
    history.append(event)
    with (root/'events.jsonl').open('a') as stream:stream.write(json.dumps(event)+'\n')
    if change<1e-11:break
else:raise RuntimeError('contact prior RHS outer solve did not converge')
assert diag['mass_flux']<1e-11 and np.array_equal(flat,initial_copy) and np.array_equal(L,factor_copy)
assert np.array_equal(rates[:,:2],gen)
actual=np.column_stack([diag['mass_rates'],rates[:,3:]]);actual[active]+=gravity[active]
reference=boundary_derivative(a[:,:2],rates[:,:2],z,edges,field,32)[:,:18]
delta=actual-reference;cross=lambda x,y:x[:,0]*y[:,1]-x[:,1]*y[:,0]
v=raw[:,7:9]/mass[:,None];c=raw[:,1:3]/mass[:,None];u=(raw[:,17]-.5*(raw[:,13]+raw[:,16]))/mass
areas=np.array([density.area for density in densities])
values=dict(bulk_mean_velocity_rate=(delta[:,7:9]-v*delta[:,0,None])/mass[:,None],centroid_velocity=(delta[:,1:3]-c*delta[:,0,None])/mass[:,None],specific_thermal_rate=(delta[:,17]-.5*(delta[:,13]+delta[:,16])-u*delta[:,0])/mass,internal_angular_rate=delta[:,10]-delta[:,11]-(cross(delta[:,1:3],raw[:,7:9])+cross(raw[:,1:3],delta[:,7:9]))/mass+cross(raw[:,1:3],raw[:,7:9])*delta[:,0]/mass**2,pressure_rate=diag['pressure_rates']+(model.GAMMA-1)*(gravity[:,17]-.5*(gravity[:,13]+gravity[:,16]))/areas)
radius=np.linalg.norm(a[:,:2],axis=1);metrics={};largest={}
for name,mask in [('active',active),('inner_r_lt_040',active&(radius<.4)),('r040_048',active&(radius>=.4)&(radius<.48)),('contact_r048_055',active&(radius>=.48))]:
    metrics[name]={}
    for key,val in values.items():
        q=val[mask].reshape(int(mask.sum()),-1)
        metrics[name][key]=float(np.sqrt(np.average(np.sum(q*q,axis=1),weights=mass[mask])))
for key,val in values.items():
    magnitude=np.linalg.norm(val.reshape(len(ids),-1),axis=1);top=np.flatnonzero(active)[np.argsort(magnitude[active])[-5:][::-1]]
    largest[key]=[dict(id=int(ids[i]),radius=float(radius[i]),norm=float(magnitude[i])) for i in top]
np.savez(root/'rhs.npz',ids=ids,mass=mass,rates=d,gravity=gravity,raw_derivative=actual,reference=reference,pressure_raw_rates=diag['pressure_raw_rates'],transport_kinematic_raw_rates=diag['transport_kinematic_raw_rates'],moment_defects=moment_defects,**values)
report=dict(history=history,metrics=metrics,largest=largest,reconstruction=reconstruction,
 maximum_reconstructed_momentum_D_Pi_defect=float(abs(moment_defects).max()),seconds=time.monotonic()-start,
 nonflat_priors=sum(not p['flat'] for p in priors),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
 scope='Complete initial contact-prior candidate: new density/basis/samples/generator/gravity, unchanged raw/factors, converged kinetic mass root; own moving-boundary exact reference.',
 limits='Research initial RHS only. No evolution/backend/long-time pass. New-density gravity uses cut quadrature order12; independent gravity refinement and potential closure remain separate gates.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
