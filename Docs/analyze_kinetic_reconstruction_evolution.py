"""Compare actual kinetic moment evolution with/without projected curvature."""
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
results=[]
initial_hashes=set()
for mode in ('affine','projected'):
    folder=root/('results_kinetic_'+mode);cases=[];states=[]
    for dt in (.002,.001,.0005):
        tag=f'alpha0_dt{dt:g}'
        meta=json.loads((folder/(tag+'.json')).read_text());cases.append(meta)
        initial_hashes.add(meta['initial_hash'])
        if meta['status']=='completed':
            saved=np.load(folder/(tag+'.npz'));assert float(saved['time'])==.05
            states.append(saved['state'].copy())
    row=dict(mode=mode,cases=cases)
    if len(states)==3:
        errors=[float(np.max(abs(states[i]-states[i+1]))) for i in (0,1)]
        row.update(max_state_differences=errors,ratio=errors[0]/errors[1],order=float(np.log2(errors[0]/errors[1])))
    results.append(row)
assert len(initial_hashes)==1
print(json.dumps(dict(scope=__doc__,results=results,
    limitations='Existing exponential factor/midpoint integrator; order measured, not assumed. Low-pressure perturbed shear, no gravity, no analytic perturbed reference. No long Kepler pass.'),indent=2))
