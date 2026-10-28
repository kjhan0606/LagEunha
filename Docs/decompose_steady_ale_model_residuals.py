"""Attribute initial model-minus-exact ALE rate defects to hydro/gravity.

Uses completed exact closure quadrature; no solver change or time evolution.
"""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
g=np.load(base/'Codex_continuous_quadratic_gravity_20261018/gravity_q24.npz')
c=np.load(base/'Codex_steady_ale_exact_closure_20261018/closure_actual_q16.npz')
assert np.array_equal(h['ids'],c['ids']) and np.array_equal(g['ids'],h['ids'])
ids=h['ids'];m=h['mass'];a=h['state'][:-4].reshape(-1,20);r=np.column_stack([m,a[:,3:]])
hydro_model=np.column_stack([np.zeros(len(m)),h['rates'][:-4].reshape(-1,20)[:,3:]])
hydro_exact=(c['advective']+c['pressure']+c['volume'])[:,:18]
dh=hydro_model-hydro_exact;dg=g['source']-c['gravity'][:,:18]
active=~h['frozen'];radius=np.linalg.norm(a[:,:2],axis=1)
v=r[:,7:9]/m[:,None];centroid=r[:,1:3]/m[:,None];U=r[:,17]-.5*(r[:,13]+r[:,16]);u=U/m
cross=lambda x,y:x[:,0]*y[:,1]-x[:,1]*y[:,0]
def derived(d):
    return dict(bulk_mean_velocity_rate=(d[:,7:9]-v*d[:,0,None])/m[:,None],
        centroid_velocity=(d[:,1:3]-centroid*d[:,0,None])/m[:,None],
        specific_thermal_rate=(d[:,17]-.5*(d[:,13]+d[:,16])-u*d[:,0])/m,
        internal_angular_rate=d[:,10]-d[:,11]-(cross(d[:,1:3],r[:,7:9])+cross(r[:,1:3],d[:,7:9]))/m+cross(r[:,1:3],r[:,7:9])*d[:,0]/m**2)
fields={name:derived(d) for name,d in [('hydro',dh),('gravity',dg),('total',dh+dg)]}
regions={}
for name,mask in [('active',active),('inner_r_lt_040',active&(radius<.4)),
                 ('r040_048',active&(radius>=.4)&(radius<.48)),('contact_r048_055',active&(radius>=.48))]:
    rows={}
    for key in fields['total']:
        rows[key]={}
        for component,values in fields.items():
            q=values[key][mask].reshape(int(mask.sum()),-1)
            rows[key][component]=dict(max_abs=float(abs(q).max()),mass_weighted_RMS=float(np.sqrt(np.average(np.sum(q*q,axis=1),weights=m[mask]))))
    regions[name]=dict(cells=int(mask.sum()),rates=rows)
largest={}
for key in fields['total']:
    magnitude=np.linalg.norm(fields['total'][key].reshape(len(ids),-1),axis=1)
    indices=np.flatnonzero(active)[np.argsort(magnitude[active])[-10:][::-1]]
    largest[key]=[dict(id=int(ids[i]),radius=float(radius[i]),
        components={k:np.asarray(f[key][i]).tolist() for k,f in fields.items()}) for i in indices]
root=base/'Codex_steady_ale_residual_decomposition_20261018';root.mkdir(exist_ok=False)
np.savez(root/'defects.npz',ids=ids,mass=m,hydro=dh,gravity=dg,active=active,raw=r)
report=dict(regions=regions,largest=largest,
 scope='Initial model-minus-independent-exact rates at the same actual face motion. Hydro includes mass-constraint and reconstruction effects. Components are signed rate defects, not causal attribution of long-time error. Frozen excluded from metrics; no numerical Psi evolution.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
