"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_contact_prior_covariance_modules_20261018'),str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
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

from contact_prior_density import ContactPriorDensity,local_contact_priors,features
from density_quadrature import rule
polys=[]
for i in ids:
    rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vertices=rows[:,6:8];polys.append(vertices[np.argsort(np.arctan2(vertices[:,1],vertices[:,0]))]+allx[i])

priors=local_contact_priors(polys,raw,edges)
angle=.73;R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]);shift=np.array([2.1,-1.7])
transformed=[];rr=raw.copy()
for i,(poly,r) in enumerate(zip(polys,raw)):
    transformed.append(poly@R.T+shift);rq=R@r[1:3]
    rr[i,1:3]=rq+r[0]*shift
    rr[i,3:7]=(R@r[3:7].reshape(2,2)@R.T+np.outer(rq,shift)+np.outer(shift,rq)+r[0]*np.outer(shift,shift)).ravel()
newpriors=local_contact_priors(transformed,rr,edges)
classification_mismatch=[int(ids[i]) for i in range(len(ids)) if priors[i]['flat']!=newpriors[i]['flat']]
assert not classification_mismatch
normal_error=0.;offset_error=0.
for old,new in zip(priors,newpriors):
    if old['flat']:continue
    expected=R@np.asarray(old['normal'])
    normal_error=max(normal_error,float(np.linalg.norm(expected-np.asarray(new['normal']))))
    offset_error=max(offset_error,abs(new['offset']-old['offset']-expected@shift))
chosen=np.flatnonzero([not p['flat'] for p in priors]);chosen=chosen[[0,len(chosen)//2,-1]]
rows=[]
for i in chosen:
    d=ContactPriorDensity(polys[i],raw[i,0],raw[i,1:3],raw[i,3:7].reshape(2,2),priors[i],12)
    e=ContactPriorDensity(transformed[i],rr[i,0],rr[i,1:3],rr[i,3:7].reshape(2,2),newpriors[i],12)
    actual=d.value(d.points);changed=e.value(d.points@R.T+shift)
    error=float(np.max(abs(actual-changed)/actual))
    assert error<5e-7
    rows.append(dict(id=int(ids[i]),relative_density_transform_error=error,theta_norm=float(np.linalg.norm(d.theta))))
assert normal_error<1e-9 and offset_error<1e-9
# An exactly aligned cut uses the cell-interior trace rather than an epsilon shift.
p=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
cut=.5;low=1.;high=4.;mass=2.5;q=np.array([1.625,1.25]);Q=np.array([[1.2083333333333333,.8125],[.8125,.8333333333333334]])
d=ContactPriorDensity(p,mass,q,Q,dict(flat=False,low=low,high=high,normal=[1.,0.],offset=cut),12)
point=np.array([cut,.4]);left=float(d.value_trace(point,[-1.,0.])[0]);right=float(d.value_trace(point,[1.,0.])[0])
assert abs(left-low)<1e-10 and abs(right-high)<1e-10
print(json.dumps(dict(passed=True,cells=len(ids),classification_mismatch=classification_mismatch,maximum_normal_error=normal_error,maximum_offset_error=offset_error,selected_density_checks=rows,aligned_trace=[left,right],scope='Inferred priors under rigid coordinate transforms; selected nonzero-theta density fits and one-sided trace.'),indent=2))
