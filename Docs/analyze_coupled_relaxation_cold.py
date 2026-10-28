"""Low-pressure moving-cell experiment; does not establish Kepler stability."""
import json
from pathlib import Path
import numpy as np

root=Path('/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018/results_cold')
out=[]
for alpha in (0,1):
    states=[];cases=[]
    for dt in (.002,.001,.0005):
        tag=f'alpha{alpha}_dt{dt:g}'
        result=json.loads((root/(tag+'.json')).read_text())
        history=json.loads((root/(tag+'.energy.json')).read_text())
        data=np.load(root/(tag+'.npz'))
        assert result['status']=='completed' and float(data['time'])==.05
        if states:assert np.array_equal(initial,data['initial'])
        initial=data['initial'].copy();states.append(data['state'].copy())
        errors=[abs(r['energy']['decomposition_residual']) for r in history]
        assert max(errors)<1e-12
        delta={k:history[-1]['energy'][k]-history[0]['energy'][k] for k in history[0]['energy']}
        cases.append(dict(dt=dt,final=result,energy_change=delta,
            energy_initial=history[0]['energy'],energy_final=history[-1]['energy'],
            max_energy_decomposition_residual=max(errors)))
    diffs=[float(np.max(abs(states[i]-states[i+1]))) for i in (0,1)]
    out.append(dict(alpha=alpha,cases=cases,max_state_differences=diffs,
        convergence_ratio=diffs[0]/diffs[1],observed_order=float(np.log2(diffs[0]/diffs[1]))))
print(json.dumps(dict(scope=__doc__,results=out,
    limitations=['No gravity; moving slip boundary exchanges energy and angular momentum.',
    'Rotation changes include pressure/transport/inertia effects; added heat is not minus rotation change.',
    'Residual-transport closure remains unvalidated.']),indent=2))
