"""Check the hard reconstruction-radius transition against saved RK states.

Independent full point-cloud fit; agreement is checked, never assumed to
prove that production candidate search used the same stencil.
"""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
r=Path('/gpfs/kjhan/LagForce/Codex_coupled_generator_stages_20261018')
out={}
for sign in ['plus','minus']:
 out[sign]={}
 for call in [2,3]:
  a=np.loadtxt(r/f'{sign}_5e-05'/f'rhs.{call}');a=a[np.argsort(a[:,0])]
  assert np.array_equal(a[:,0],np.arange(len(a)))
  tree=cKDTree(a[:,1:3]);rows={}
  for i in [103,154,155,156,312,734,848,849,863]:
   hs=np.sqrt(a[i,6]);ids=np.array(tree.query_ball_point(a[i,1:3],2.5*hs));ids=ids[ids!=i]
   dx=a[ids,1:3]-a[i,1:3];dist=np.linalg.norm(dx,axis=1);order=np.argsort(dist);ids=ids[order];dx=dx[order];dist=dist[order]
   X,Y=(dx/hs).T;sw=hs/dist
   design=sw[:,None]*np.column_stack((X,Y,.5*X*X,X*Y/np.sqrt(2),.5*Y*Y,X**3/6,X*X*Y/(2*np.sqrt(3)),X*Y*Y/(2*np.sqrt(3)),Y**3/6))
   coeff,residual,rank,singular=np.linalg.lstsq(design,sw[:,None]*(a[ids,3:5]-a[i,3:5]),rcond=None)
   fitted=(coeff[:2].T/hs).reshape(-1);actual=a[i,21:25]
   near=np.array(tree.query_ball_point(a[i,1:3],2.51*hs));near=near[near!=i]
   ratio=np.linalg.norm(a[near,1:3]-a[i,1:3],axis=1)/hs
   edge=np.argsort(abs(ratio-2.5))[:4]
   rows[str(i)]={'selected_ids':ids.tolist(),'rank':int(rank),'condition_number':float(singular[0]/singular[-1]),
    'fit_gradient':fitted.tolist(),'production_gradient':actual.tolist(),'max_gradient_abs_difference':float(abs(fitted-actual).max()),
    'nearest_cutoff_candidates':[{'id':int(near[j]),'distance_over_sqrtV':float(ratio[j])} for j in edge]}
  out[sign][str(call)]=rows
transitions={}
for sign in out:
 transitions[sign]={}
 for i in out[sign]['2']:
  before=out[sign]['2'][i];after=out[sign]['3'][i]
  transitions[sign][i]={'added':sorted(set(after['selected_ids'])-set(before['selected_ids'])),
   'removed':sorted(set(before['selected_ids'])-set(after['selected_ids'])),
   'max_fit_error':max(before['max_gradient_abs_difference'],after['max_gradient_abs_difference']),
   'ranks':[before['rank'],after['rank']]}
result={'states':out,'transitions':transitions,'scope':'Independent cubic weighted least-squares fits using saved stage positions and velocities, radius2.5sqrt(V). Interior cells only; compares to rounded production gradients. No production candidate-list export or trajectory modification.'}
with (r/'reconstruction_shell_analysis.json').open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(transitions,indent=2))
