"""Required intrinsic angular-momentum balance from existing ALE packets.

Midpoint allocation below is a diagnostic convention, not an exact integral
of angular-momentum flux over the face and not an implemented spin equation.
The *sum* required to cancel the pair orbital torque is independent of that
common application point.
"""
import json
import numpy as np
from analyze_ale_shadow_v1 import ROOT,FOCUS,load

target=ROOT/'rk4_diagnostics/ale_shadow_real_v1/angular_balance.json'
if target.exists():raise SystemExit('Refuse to overwrite angular analysis')
result={}
for name in ('fine_captured','fine_fresh_final'):
    p,c,a,b,metas=load(target.parent/name);z=p['a'];n=len(c)
    dx=(c['x'][b]-c['x'][a]+12)%24-12
    dy=(c['y'][b]-c['y'][a]+12)%24-12
    ox=.5*(z[:,0]+z[:,2]);oy=.5*(z[:,1]+z[:,3])
    values={}
    for label,start in (('full',5),('pressure',9),('advection',13),('riemann_residual',17),('native_canonical',21)):
        fx=z[:,start+1];fy=z[:,start+2]
        orbit=dx*fy-dy*fx
        si=-(ox*fy-oy*fx);sj=(ox-dx)*fy-(oy-dy)*fx
        total=orbit+si+sj
        cells=np.zeros(n);np.add.at(cells,a,si);np.add.at(cells,b,sj)
        interior_pairs=(np.hypot(c['x'][a]-12,c['y'][a]-12)<.5)&(np.hypot(c['x'][b]-12,c['y'][b]-12)<.5)
        values[label]=dict(max_abs_pair_orbital_rate=float(abs(orbit).max()),
            max_pair_orbit_plus_midpoint_intrinsic_residual=float(abs(total).max()),
            core_both_endpoints_Rlt_point5_abs_pair_orbital_rate_sum=float(abs(orbit[interior_pairs]).sum()),
            focus_midpoint_intrinsic_rate={str(ident):float(cells[np.flatnonzero(c['id']==ident)[0]]) for ident in FOCUS})
    result[name]=dict(scope='instantaneous angular balance of hypothetical full-flux momentum update; no state/spin/energy changes; periodic seam excluded from core summary',
        identity='orbital_pair=(xj-xi) cross Fp; intrinsic_i=-(xf-xi) cross Fp; intrinsic_j=+(xf-xj) cross Fp; sum zero for any common xf',
        limitation='Per-cell face-midpoint allocation is approximate when flux varies along the face; pair total cancellation is exact algebra. A finite-energy spin/affine closure is still needed.',terms=values)
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(output=str(target),cases={name:case['terms'] for name,case in result.items()}),indent=2),flush=True)
