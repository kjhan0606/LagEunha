"""Compare completed instantaneous RHS cases on common cells, no evolution."""
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_continuous_initial_rhs_20261018')
report=json.loads((root/'report.json').read_text())
cases={}
for row in report['results']:
    if row['status']=='completed':
        key=(row['cutoff'],row['face_order'])
        cases[key]=np.load(root/f'rhs_r{key[0]:g}_f{key[1]}.npz')
comparison=[]
def compare(a,b,mask,label):
    aa,bb=cases[a],cases[b];assert np.array_equal(aa['ids'],bb['ids'])
    ar=aa['rates'][:-4].reshape(-1,20);br=bb['rates'][:-4].reshape(-1,20)
    out=dict(label=label,cases=[a,b],cells=int(mask.sum()))
    for name,va,vb in [('momentum_rate',ar[:,9:11],br[:,9:11]),('hydro_energy_rate',ar[:,19],br[:,19]),('pressure_rate',aa['pressure_rates'],bb['pressure_rates']),('weight_rate',ar[:,2],br[:,2])]:
        diff=np.linalg.norm((vb-va)[mask]);norm=np.linalg.norm(vb[mask])
        out[name]=dict(difference_l2=float(diff),relative_to_second_l2=float(diff/norm) if norm else None,maximum_absolute_difference=float(np.max(abs((vb-va)[mask]))))
    comparison.append(out)
for cutoff in (.45,.5,.55):
    a,b=(cutoff,4),(cutoff,8)
    if a in cases and b in cases:compare(a,b,~cases[a]['frozen'],'face_quadrature')
for cutoff in (.45,.5):
    a,b=(cutoff,8),(.55,8)
    if a in cases and b in cases:
        x=cases[a]['state'][:-4].reshape(-1,20)[:,:2]
        compare(a,b,np.linalg.norm(x,axis=1)<.4,'frozen_boundary_common_inner')
print(json.dumps(dict(scope=__doc__,source_report=report,comparisons=comparison,
    limits='Face order4/8 difference is not an established asymptotic convergence order. Common inner r<.4 comparison does not assess physical outer-boundary evolution. Weight gauge anchored to frozen cells.'),indent=2))
