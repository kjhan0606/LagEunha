"""Exact rational topology of original binary64 generator and W inputs."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_production_moment_geometry_20261018');a=np.load(root/'stage0_inputs.npz');lookup={int(g):i for i,g in enumerate(a['selected'])}
def pair(p):return [str(p.numerator),str(p.denominator)]
def run(gid):
 row=lookup[gid];xi,yi=map(F,a['xy'][gid]);wi=F(a['W'][gid]);p=[(F(-3)-xi,F(-3)-yi),(F(3)-xi,F(-3)-yi),(F(3)-xi,F(3)-yi),(F(-3)-xi,F(3)-yi)];tags=[-1,-2,-3,-4]
 for j in a['candidates'][a['offsets'][row]:a['offsets'][row+1]]:
  j=int(j);xj,yj=map(F,a['xy'][j]);wj=F(a['W'][j]);dx=xj-xi;dy=yj-yi;threshold=(dx*dx+dy*dy+wi-wj)/2
  q=[];qt=[]
  for e,pa in enumerate(p):
   pb=p[(e+1)%len(p)];da=pa[0]*dx+pa[1]*dy;db=pb[0]*dx+pb[1]*dy;ai=da<=threshold;bi=db<=threshold
   if ai:q.append(pa);qt.append(tags[e])
   if ai!=bi:
    t=(threshold-da)/(db-da);q.append((pa[0]+t*(pb[0]-pa[0]),pa[1]+t*(pb[1]-pa[1])));qt.append(j if ai else tags[e])
  clean=[];ct=[]
  for pp,tt in zip(q,qt):
   if clean and pp==clean[-1]:ct[-1]=tt
   else:clean.append(pp);ct.append(tt)
  if len(clean)>1 and clean[-1]==clean[0]:clean.pop();ct.pop()
  p=clean;tags=ct
 edges=[]
 for e,j in enumerate(tags):
  pa=p[e];pb=p[(e+1)%len(p)];length2=(pb[0]-pa[0])**2+(pb[1]-pa[1])**2
  edges.append(dict(neighbor=j,exact_positive_length=length2>0,length_approx=float(length2)**.5,length_squared_exact=pair(length2),A_global_exact=[pair(pa[0]+xi),pair(pa[1]+yi)],B_global_exact=[pair(pb[0]+xi),pair(pb[1]+yi)],A_relative_float=list(map(float,pa)),B_relative_float=list(map(float,pb))))
 return dict(owner=gid,candidates=a['candidates'][a['offsets'][row]:a['offsets'][row+1]].tolist(),edges=edges)
report=dict(scope=__doc__,cells=[run(3),run(5)],limits='Exact clipping against saved candidate lists, not exact all-registry completeness or production implementation. No RHS or evolution.')
(root/'exact_contact_diagnostic.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
