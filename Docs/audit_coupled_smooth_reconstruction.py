"""Local candidate only: smoothly taper LSQ weights at the existing radius.

No production dynamics, particle motion, or cell weights are modified.
"""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
r=Path('/gpfs/kjhan/LagForce/Codex_coupled_generator_stages_20261018')
out={}
for sign in ['plus','minus']:
 states={};selected={}
 for call in [2,3]:
  a=np.loadtxt(r/f'{sign}_5e-05'/f'rhs.{call}');a=a[np.argsort(a[:,0])];states[call]=a
  tree=cKDTree(a[:,1:3]);selected[call]={i:set(tree.query_ball_point(a[i,1:3],2.5*np.sqrt(a[i,6])))-{i} for i in [155,848]}
 out[sign]={}
 for i in [155,848]:
  union=selected[2][i]|selected[3][i];fits={}
  for call,a in states.items():
   fits[str(call)]={}
   for method in ['hard','fixed_union','smooth']:
    ids=np.array(sorted(union if method=='fixed_union' else selected[call][i]));hs=np.sqrt(a[i,6])
    xy=(a[ids,1:3]-a[i,1:3])/hs;X,Y=xy.T;radius=np.linalg.norm(xy,axis=1)
    basis=np.column_stack((X,Y,X*X/2,X*Y/np.sqrt(2),Y*Y/2,X**3/6,X*X*Y/(2*np.sqrt(3)),X*Y*Y/(2*np.sqrt(3)),Y**3/6))
    sw=1/radius
    if method=='smooth':sw=sw*np.maximum(1-(radius/2.5)**2,0)**2
    design=sw[:,None]*basis
    coef,_,rank,singular=np.linalg.lstsq(design,sw[:,None]*(a[ids,3:5]-a[i,3:5]),rcond=None)
    assert rank==9
    identity,_,_,_=np.linalg.lstsq(design,design,rcond=None)
    reproduction=float(abs(identity-np.eye(9)).max());assert reproduction<1e-10
    fits[str(call)][method]={'gradient':(coef[:2].T/hs).reshape(-1).tolist(),'rank':int(rank),
     'condition_number':float(singular[0]/singular[-1]),'polynomial_coefficient_reproduction_error':reproduction}
  changes={k:float(np.linalg.norm(np.array(fits['3'][k]['gradient'])-fits['2'][k]['gradient'])) for k in ['hard','fixed_union','smooth']}
  out[sign][str(i)]={'fits':fits,'gradient_change_norm':changes,'removed':sorted(selected[2][i]-selected[3][i])}
result={'cells':out,'candidate':'LSQ row multiplier (1-(d/(2.5sqrtV))^2)^2 inside support; zero outside. Thus squared least-squares weight has fourth-power taper.',
 'scope':'Local independent reconstruction diagnostic, not an adopted solver. Fixed-union is a diagnostic control only. Smooth fit preserves cubic polynomial reproduction in these full-rank stencils. Boundary/fallback/rank transitions, conditioning elsewhere, shock response, conservation and long-time stability remain untested. Reconstruction weights are distinct from Laguerre cell weights.'}
with (r/'smooth_reconstruction_analysis.json').open('x') as f:json.dump(result,f,indent=2)
print(json.dumps({s:{i:{'gradient_change_norm':v['gradient_change_norm'],'removed':v['removed']} for i,v in cells.items()} for s,cells in out.items()},indent=2))
