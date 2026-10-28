from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_production_moment_geometry_20261018');r=json.loads((root/'exact_contact_all_registry_diagnostic.json').read_text());a=np.load(root/'stage0_inputs.npz');witness=[]
def unpack(x):return F(int(x[0]),int(x[1]))
def pack(x):return [str(x.numerator),str(x.denominator)]
for round in r['rounds']:
 for c in round['cells']:
  i=c['owner'];xi,yi=map(F,a['xy'][i]);wi=F(a['W'][i])
  for j in c['violating_global_ids']:
   xj,yj=map(F,a['xy'][j]);dx=xj-xi;dy=yj-yi;rhs=dx*dx+dy*dy+wi-F(a['W'][j])
   for e in c['edges']:
    pa=list(map(unpack,e['A_global_exact']));residual=2*((pa[0]-xi)*dx+(pa[1]-yi)*dy)-rhs
    if residual>0:witness.append(dict(iteration=round['iteration'],owner=i,missing_neighbor=j,vertex_global_exact=e['A_global_exact'],outgoing_tag=e['neighbor'],plane_dx_exact=pack(dx),plane_dy_exact=pack(dy),plane_rhs_relative_exact=pack(rhs),positive_violation_exact=pack(residual),violation_float=float(residual)))
ea=next(e for e in r['cells'][0]['edges'] if e['neighbor']==5);eb=next(e for e in r['cells'][1]['edges'] if e['neighbor']==3)
assert ea['A_global_exact']==eb['B_global_exact'] and ea['B_global_exact']==eb['A_global_exact']
report=dict(all_registry_count=r['all_registry_count'],candidate_completion_iterations=len(r['rounds'])-1,added_global_ids={'3':[],'5':[2]},witnesses=witness,final_edge_3_5_exact_positive=True,final_edge_length_approx=ea['length_approx'],global_endpoints_exact_reversed_equal=True,endpoint_A=ea['A_global_exact'],endpoint_B=ea['B_global_exact'],scope='Exact original binary64 interpretation; two stored stage0 cells, not full geometry validation.')
(root/'exact_contact_witness.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
