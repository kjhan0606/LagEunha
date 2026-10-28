"""Differentiate covariance of stored moments including gravity and hydro.
Euler tangent positivity bound is diagnostic, not a production step prescription.
"""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');rows=[]
def covariance(r):
    m=r[0];q=r[1:3];Q=r[3:7].reshape(2,2);p=r[7:9];D=r[9:13].reshape(2,2);Pi=r[13:17].reshape(2,2)
    S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m;B=np.linalg.solve(S,C).T
    R=Pi-np.outer(p,p)/m-B@S@B.T
    return (R+R.T)/2,B
for mode,hydro_root,gravity_root in [('affine','Codex_continuous_initial_rhs_split_20261018','Codex_continuous_initial_gravity_v2_20261018'),('quadratic','Codex_continuous_initial_rhs_quadratic_split_20261018','Codex_continuous_quadratic_gravity_20261018')]:
    h=np.load(base/hydro_root/'rhs_r0.55_f8.npz');g=np.load(base/gravity_root/'gravity_q24.npz');mass=h['mass'];y=h['state'][:-4].reshape(-1,20)
    raw=np.column_stack([mass,y[:,3:]]);drate=g['combined_rates'][:-4].reshape(-1,20);d=np.column_stack([np.zeros(len(mass)),drate[:,3:]])
    active=np.flatnonzero(~h['frozen']);derivatives=[];bounds=[];fdmax=np.zeros(2)
    for i in active:
        r=raw[i];v=d[i];m=mass[i];q=r[1:3];p=r[7:9];dq=v[1:3];dp=v[7:9];R,B=covariance(r)
        dS=v[3:7].reshape(2,2)-(np.outer(dq,q)+np.outer(q,dq))/m
        dC=v[9:13].reshape(2,2)-(np.outer(dq,p)+np.outer(q,dp))/m
        dT=v[13:17].reshape(2,2)-(np.outer(dp,p)+np.outer(p,dp))/m
        dR=dT-B@dC-dC.T@B.T+B@dS@B.T;dR=(dR+dR.T)/2;derivatives.append(dR)
        L=h['covariance_factors'][i];normalized=np.linalg.solve(L,np.linalg.solve(L,dR).T).T
        minimum=np.linalg.eigvalsh((normalized+normalized.T)/2).min()
        if minimum<0:bounds.append((float(-1/minimum),int(h['ids'][i])))
        for k,eps in enumerate((1e-5,5e-6)):
            fd=(covariance(r+eps*v)[0]-covariance(r-eps*v)[0])/(2*eps)
            fdmax[k]=max(fdmax[k],float(np.max(abs(fd-dR))))
    derivatives=np.array(derivatives)
    rows.append(dict(mode=mode,active_cells=len(active),maximum_covariance_rate=float(np.max(abs(derivatives))),
        finite_difference_max_errors=fdmax.tolist(),negative_normalized_rate_cells=len(bounds),
        smallest_euler_tangent_bounds=sorted(bounds)[:5]))
print(json.dumps(dict(scope=__doc__,results=rows,limits='Fixed-mass derivative, actual rates; no covariance clipping. A finite tangent bound does not prove loss of positivity under the exact nonlinear evolution, nor stability below the bound. Finite differences suffer raw-moment cancellation.'),indent=2))
