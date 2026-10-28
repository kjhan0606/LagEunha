"""Compare quadratic split quadrature and affine hydro on identical moments."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_continuous_initial_rhs_quadratic_split_20261018'
report=json.loads((root/'report.json').read_text())
cases={'quadratic4':np.load(root/'rhs_r0.55_f4.npz'),'quadratic8':np.load(root/'rhs_r0.55_f8.npz'),
       'affine8':np.load(base/'Codex_continuous_initial_rhs_split_20261018/rhs_r0.55_f8.npz')}
rows=[]
for a,b in [('quadratic4','quadratic8'),('affine8','quadratic8')]:
    left,right=cases[a],cases[b];assert np.array_equal(left['state'],right['state'])
    ar=left['rates'][:-4].reshape(-1,20);br=right['rates'][:-4].reshape(-1,20)
    x=left['state'][:-4].reshape(-1,20)[:,:2];active=~left['frozen']
    row=dict(cases=[a,b],regions={})
    for region,mask in [('active',active),('inner',np.linalg.norm(x,axis=1)<.4)]:
        out={}
        for name,u,v in [('pressure_rate',left['pressure_rates'],right['pressure_rates']),('momentum_rate',ar[:,9:11],br[:,9:11]),('energy_rate',ar[:,19],br[:,19]),('weight_rate',ar[:,2],br[:,2])]:
            un=np.linalg.norm(u[mask]);vn=np.linalg.norm(v[mask]);dn=np.linalg.norm((v-u)[mask])
            out[name]=dict(first_l2=float(un),second_l2=float(vn),difference_l2=float(dn),relative_difference=float(dn/vn) if vn else None)
        row['regions'][region]=out
    rows.append(row)
print(json.dumps(dict(scope=__doc__,report=report,comparisons=rows,
    limits='Instantaneous hydro only. Smaller pressure rate is not a long-time accuracy proof. Generator paths differ because reconstruction changes their material velocity; no steering imposed.'),indent=2))
