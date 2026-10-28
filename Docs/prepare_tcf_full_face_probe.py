"""Recover finite interior Voronoi faces independently from saved t6 positions."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi

base=Path('/gpfs/kjhan/LagForce')
source=base/'Codex_TCF_rhs_identity_20261017/t6'
s=np.loadtxt(source/'rhs.dat'); s=s[np.argsort(s[:,0])]
assert np.array_equal(s[:,0],np.arange(len(s)))
f=np.loadtxt(source/'rhs.dat.faces')
x=s[:,1:3]; v=s[:,3:5]; P=s[:,5]
mesh=Voronoi(x)
ridges={tuple(sorted(pair)): vertices for pair,vertices in zip(mesh.ridge_points,mesh.ridge_vertices)}
active=np.linalg.norm(x-3,axis=1)<1.9
rows=[]; owners=[]; area_err=[]; missing=[]
F=np.zeros((len(s),2)); E=np.zeros(len(s)); V=np.zeros(len(s)); H=np.zeros(len(s))
for rec in f:
    i,j=rec[:2].astype(int)
    if not active[i]: continue
    verts=ridges.get(tuple(sorted((i,j))))
    if verts is None or len(verts)!=2 or min(verts)<0:
        missing.append([int(i),int(j)]); continue
    ends=mesh.vertices[verts]; A=np.linalg.norm(ends[1]-ends[0])
    dx=x[j]-x[i]; d=np.linalg.norm(dx); e=dx/d
    t=ends.mean(axis=0)-.5*(x[i]+x[j]); t-=np.dot(t,e)*e
    Q=rec[14]; ui=v[i]@e; uj=v[j]@e
    force=-A*((.5*(P[i]+P[j])+Q)*e+(P[j]-P[i])*t/d)
    energy=-A*(.5*(P[j]*ui+P[i]*uj)+(P[j]*(v[i]@t)-P[i]*(v[j]@t))/d+.5*Q*(ui+uj))
    vol=A*(.5*(ui+uj)+(v[i]-v[j])@t/d)
    heat=.5*A*Q*(ui-uj)
    F[i]+=force; E[i]+=energy; V[i]+=vol; H[i]+=heat
    rows.append([A,d,*e,*t,P[i],P[j],*v[i],*v[j],Q]); owners.append(i)
    area_err.append(abs(A-rec[4])/max(A,1e-300))
def compare(a,b):
    a=a[active]; b=b[active]
    return {'max_abs':float(np.max(abs(a-b))),
            'relative_l2':float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))}
result={'source':str(source),'center':[3,3],'active_cells':int(active.sum()),
    'faces':len(rows),'missing_face_count':len(missing),'missing_first10':missing[:10],
    'area_relative_max':float(max(area_err)),
    'force':compare(F,s[:,8:10]),'energy':compare(E,s[:,10]),
    'volume':compare(V,s[:,12]),'heat':compare(H,s[:,11]),
    'scope':'Independent SciPy Voronoi on saved full-precision positions; r<1.9. Recorded rounded Q. No QBP replay, weights or production geometry kernel.'}
root=base/'Codex_tcf_full_face_probe_20261018'; root.mkdir(exist_ok=False)
np.asarray(rows).tofile(root/'inputs.bin')
np.asarray(owners,dtype=np.int64).tofile(root/'owners.bin')
(root/'geometry_reference.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
