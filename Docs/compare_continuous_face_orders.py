"""Saved RHS quadrature convergence for identical initial state and boundary."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
old=base/'Codex_continuous_initial_rhs_20261018'
new=base/'Codex_continuous_initial_rhs_faces_20261018'
data={k:np.load((old if k<16 else new)/f'rhs_r0.55_f{k}.npz') for k in (4,8,16,32)}
rows=[]
for a,b in ((4,8),(8,16),(16,32)):
    left,right=data[a],data[b]
    assert np.array_equal(left['state'],right['state']) and np.array_equal(left['mass'],right['mass'])
    active=~left['frozen'];ar=left['rates'][:-4].reshape(-1,20);br=right['rates'][:-4].reshape(-1,20)
    row=dict(orders=[a,b],active_cells=int(active.sum()))
    for name,u,v in [('pressure_rate',left['pressure_rates'],right['pressure_rates']),('momentum_rate',ar[:,9:11],br[:,9:11]),('energy_rate',ar[:,19],br[:,19]),('weight_rate',ar[:,2],br[:,2])]:
        diff=(v-u)[active];norm=np.linalg.norm(v[active])
        row[name]=dict(difference_l2=float(np.linalg.norm(diff)),relative_l2=float(np.linalg.norm(diff)/norm),max_abs=float(np.max(abs(diff))))
    rows.append(row)
print(json.dumps(dict(scope=__doc__,comparisons=rows,limits='Successive quadrature differences, no continuum or time convergence. Fixed artificial boundary at radius .55, no gravity.'),indent=2))
