"""Read completed checkpoints; audit intermediate geometric masses only."""
import json
from pathlib import Path
import numpy as np
from spacetime_ale_mass import endpoint,velocity
rows=[]
for T,seed in [(.02,0),(.1,0),(.3,0),(.1,731),(.3,731),(.3,732)]:
    c=json.loads(Path(f'T{T:g}_seed{seed}.checkpoint.json').read_text())
    x=np.array(c['sites_initial']);w=np.array(c['weights_initial']);dw=np.array(c['weight_increment'])
    mass=np.array(c['initial'])[:,0];drift=[]
    for tau in (.25,.5,.75):
        actual=endpoint(x+tau*T*velocity(x),w+tau*dw,tau*T)[:,0]
        drift.append(float(np.max(abs(actual/mass-1))))
    rows.append(dict(duration=T,seed=seed,sample_tau=[.25,.5,.75],relative_mass_drift=drift))
print(json.dumps(rows,indent=2))
