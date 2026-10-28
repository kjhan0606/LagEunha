"""Full initial moment-rate comparison against moving-cell steady Kepler.

Reuses completed converged hydro and gravity evaluations; no evolution.
Exact boundary quadrature splits both circular density interfaces.
"""
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
sys.path[:0]=[str(base/'Codex_contact_prior_density_modules_20261018'),str(base/'Codex_steady_ale_reference_modules_20261018'),str(base)]
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
root=base/'Codex_contact_prior_actual_fit_20261018';root.mkdir(exist_ok=False)
(root/'priors.json').write_text(json.dumps(priors,indent=2))
results=[]
for nq,ncheck in ((12,24),(24,32)):
    records=[];failures=[]
    for i,(poly,r,meta) in enumerate(zip(polys,raw,priors)):
        try:
            d=ContactPriorDensity(poly,r[0],r[1:3],r[3:7].reshape(2,2),meta,nq)
            ps=[];ws=[]
            for part in d.parts:
                p,w=rule(part,ncheck);ps.append(p);ws.append(w)
            p=np.concatenate(ps);w=np.concatenate(ws);rho=d.value(p);dm=w*rho
            normalized_error=float(abs((dm@features(p,d.origin,d.scale)-d.weighted_mass@features(d.points,d.origin,d.scale))/r[0]).max())
            recovered=np.r_[dm.sum(),dm@p,(p.T@(dm[:,None]*p)).ravel()]
            moment_error=float(abs(recovered-r[:7]).max())
            sample=np.vstack([p,poly,(poly+np.roll(poly,-1,axis=0))/2]);values=d.value(sample)
            low=meta.get('value',meta['low']);high=meta.get('value',meta['high'])
            records.append(dict(id=int(ids[i]),fit_residual=d.residual,iterations=d.iterations,condition=d.condition,
                independent_normalized_moment_error=normalized_error,independent_raw_moment_error=moment_error,
                min_density=float(values.min()),max_density=float(values.max()),
                ratio_to_prior_low=float(values.min()/low),ratio_to_prior_high=float(values.max()/high),
                theta=d.theta.tolist(),logZ=d.logZ,origin=d.origin.tolist(),scale=d.scale))
        except Exception as exc:failures.append(dict(id=int(ids[i]),error=repr(exc)))
        if (i+1)%100==0:
            with (root/'events.jsonl').open('a') as f:f.write(json.dumps(dict(order=nq,cells=i+1,failures=len(failures)))+'\n')
    (root/f'fits_q{nq}.json').write_text(json.dumps(dict(records=records,failures=failures),indent=2))
    results.append(dict(order=nq,independent_order=ncheck,successful=len(records),failures=failures,
        max_fit_residual=max((r['fit_residual'] for r in records),default=None),
        maximum_independent_normalized_moment_error=max((r['independent_normalized_moment_error'] for r in records),default=None),
        maximum_independent_raw_moment_error=max((r['independent_raw_moment_error'] for r in records),default=None),
        minimum_density=min((r['min_density'] for r in records),default=None),
        maximum_density=max((r['max_density'] for r in records),default=None),
        maximum_prior_high_overshoot=max((r['ratio_to_prior_high'] for r in records),default=None)))
report=dict(cells=len(ids),nonflat_priors=sum(not p['flat'] for p in priors),results=results,
 scope='Current raw moments and exactly two-hop local means only; no analytic interface/phase information used in inference or fit.',
 limits='Density candidate fit, not full velocity/gravity reconstruction, RHS, evolution or backend validation. Prior is not a hard density bound. Failures are explicit, no fallback/clipping.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
