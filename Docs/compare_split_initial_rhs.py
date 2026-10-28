"""Compare converged split-face orders and the saved unsplit order32 RHS."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');split=base/'Codex_continuous_initial_rhs_split_20261018'
report=json.loads((split/'report.json').read_text())
assert all(r['status']=='completed' and r['split_history'][-1]['weight_rate_change']<1e-11 for r in report['results'])
cases={'split4':np.load(split/'rhs_r0.55_f4.npz'),'split8':np.load(split/'rhs_r0.55_f8.npz'),
       'unsplit32':np.load(base/'Codex_continuous_initial_rhs_faces_20261018/rhs_r0.55_f32.npz')}
rows=[]
for a,b in [('split4','split8'),('unsplit32','split8')]:
    left,right=cases[a],cases[b];assert np.array_equal(left['state'],right['state'])
    active=~left['frozen'];ar=left['rates'][:-4].reshape(-1,20);br=right['rates'][:-4].reshape(-1,20)
    row=dict(cases=[a,b])
    for name,u,v in [('pressure_rate',left['pressure_rates'],right['pressure_rates']),('momentum_rate',ar[:,9:11],br[:,9:11]),('energy_rate',ar[:,19],br[:,19]),('weight_rate',ar[:,2],br[:,2])]:
        diff=(v-u)[active];norm=np.linalg.norm(v[active])
        row[name]=dict(relative_l2=float(np.linalg.norm(diff)/norm),max_abs=float(np.max(abs(diff))))
    rows.append(row)
print(json.dumps(dict(scope=__doc__,report=report,comparisons=rows,
    limits='Instantaneous affine hydro at t0 only. No boundary independence, gravity, or time evolution established.'),indent=2))
