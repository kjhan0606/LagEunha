import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_production_moment_geometry_20261018');data=np.load(root/'stage0_inputs.npz');row=3;i=data['selected'][row];xy=data['xy'];W=data['W'];bounds=data['bounds'];p=np.array([[bounds[0],bounds[2]],[bounds[1],bounds[2]],[bounds[1],bounds[3]],[bounds[0],bounds[3]]])-xy[i];tags=[-1,-2,-3,-4]
for j in data['candidates'][data['offsets'][row]:data['offsets'][row+1]]:
 d=xy[j]-xy[i];threshold=.5*(d@d+W[i]-W[j]);points=[];newtags=[]
 for e,a in enumerate(p):
  b=p[(e+1)%len(p)];da=a@d;db=b@d;ai=da<=threshold;bi=db<=threshold
  if ai:points.append(a.copy());newtags.append(tags[e])
  if ai!=bi:
   t=(threshold-da)/(db-da);points.append(a+t*(b-a));newtags.append(int(j) if ai else tags[e])
 q=[];qt=[]
 for a,t in zip(points,newtags):
  if len(q) and np.array_equal(a,q[-1]):qt[-1]=t
  else:q.append(a);qt.append(t)
 if len(q)>1 and np.array_equal(q[-1],q[0]):q.pop();qt.pop()
 p=np.array(q);tags=qt
 print(json.dumps(dict(cut=int(j),polygon=p.tolist(),tags=tags)),flush=True)
for e,j in enumerate(tags):
 a=p[e];b=p[(e+1)%len(p)];d=xy[j]-xy[i];edge=b-a
 print(json.dumps(dict(edge=e,neighbor=j,length=float(np.linalg.norm(edge)),outward_dot=float(edge[1]*d[0]-edge[0]*d[1]),plane_residual=(2*p[[e,(e+1)%len(p)]]@d-(d@d+W[i]-W[j])).tolist())),flush=True)
