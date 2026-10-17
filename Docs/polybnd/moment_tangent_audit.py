"""Static realizability tangent audit at zero-velocity pressure states."""
import json
import numpy as np
from adaptive_moment_batch import setup
import moment_ale_evolve as model

rows=[]
for amplitude in (.001,.1):
    mass,y=setup('wave',3,amplitude)
    dy,diag=model.rhs(0.,y,mass);d=dy[:-4].reshape(len(mass),20)
    # At v=B=Pi=0, Rdot = Pidot. Packing: x(2), w(1), raw[1:](17).
    pidot=d[:,15:19].reshape(-1,2,2)
    eig=np.array([np.linalg.eigvalsh((a+a.T)/2)/m for a,m in zip(pidot,mass)])
    rows.append(dict(amplitude=amplitude,min_Rdot=float(eig.min()),max_Rdot=float(eig.max()),
                     negative_cells=int(sum(np.min(eig,axis=1)<-1e-12)),cells=len(mass)))
print(json.dumps(rows,indent=2))
