"""Read immutable completed volume quadratures, separating frozen cells."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_steady_ale_exact_closure_20261018'
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
g=np.load(base/'Codex_continuous_quadratic_gravity_20261018/gravity_q24.npz')
rows=[]
for nq in (8,16):
    p=root/f'volume_q{nq}.npz'
    if not p.exists():continue
    v=np.load(p);assert np.array_equal(v['ids'],h['ids'])
    difference=g['source']-v['gravity'][:,:18]
    regions={}
    for name,mask in [('active',~h['frozen']),('frozen',h['frozen'])]:
        regions[name]=dict(cells=int(mask.sum()),evaluated_model_minus_exact_gravity_max=float(abs(difference[mask]).max()),
                           source_applied_in_evolution=name=='active')
    cases=[]
    for motion in ('actual','exact_mass'):
        p=root/f'closure_{motion}_q{nq}.npz'
        if not p.exists():continue
        c=np.load(p);defect=c['defect'];active=~h['frozen']
        cases.append(dict(motion=motion,active_max=float(abs(defect[active]).max()),
            frozen_max=float(abs(defect[~active]).max()),column_active_max=np.max(abs(defect[active]),axis=0).tolist()))
    rows.append(dict(order=nq,regions=regions,closure=cases))
print(json.dumps(dict(results=rows,ready_for_two_order_evaluation=len(rows)==2 and all(len(r['closure'])==2 for r in rows),
 note='Saved gravity source was evaluated for all cells; evolution applies it only to active cells. Frozen evaluated-source differences are diagnostic only. Internal torque algebra is not an independent source validation.'),indent=2))
