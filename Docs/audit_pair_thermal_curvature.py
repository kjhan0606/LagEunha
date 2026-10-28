"""Finite directional curvature of fixed-entropy Voronoi internal energy."""
import json,argparse
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi,cKDTree
ap=argparse.ArgumentParser();ap.add_argument('--tangent',action='store_true');ap.add_argument('--diagonal',action='store_true');args=ap.parse_args()
base=Path('/gpfs/kjhan/LagForce');p=base/'Codex_qbp_t10_decomposition_20261018'
s=np.loadtxt(p/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
x=s[:,1:3];r=np.linalg.norm(x-3,axis=1);gamma=float(np.float32(1.666667))
mask=r<2.4;ids=np.where(mask)[0];inside=np.where(r<.6)[0]
dist,near=cKDTree(x).query(x[inside],k=2);k=np.argmin(dist[:,1]);i=int(inside[k]);j=int(near[k,1])
direction=(x[j]-x[i])/np.linalg.norm(x[j]-x[i]);motion=np.zeros_like(x)
if args.tangent:direction=np.array([-direction[1],direction[0]])
if args.diagonal:direction=(direction+np.array([-direction[1],direction[0]]))/np.sqrt(2)
motion[i]=-s[j,7]/(s[i,7]+s[j,7])*direction;motion[j]=s[i,7]/(s[i,7]+s[j,7])*direction
def volumes(z):
 mesh=Voronoi(z);v=np.empty(len(ids))
 for k,a in enumerate(ids):
  reg=mesh.regions[mesh.point_region[a]];assert reg and min(reg)>=0
  q=mesh.vertices[reg]-z[a];v[k]=abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2
 return v
V0=volumes(x);P=s[ids,5];U=P*V0/(gamma-1)
gap=float(dist[k,1]);h=float(np.sqrt(.5*(s[i,6]+s[j,6])))
out={'scope':'Fixed-entropy internal-energy curvature only; geometry rebuilt for a center-of-mass-preserving pair separation. No gravity/rotation/AV/HG dynamics. Independent SciPy volumes; all cells r<2.4 included.',
 'direction':'diagonal' if args.diagonal else ('tangent' if args.tangent else 'radial'), 'pair':[i,j],'gap':gap,'gap_over_h':gap/h,'initial_volume_relative_max_error':float(np.max(abs(V0/s[ids,6]-1))), 'samples':[]}
for fraction in [1e-3,5e-4,2.5e-4]:
 eps=fraction*gap;vp=volumes(x+eps*motion);vm=volumes(x-eps*motion)
 # Cellwise expm1 avoids cancellation against the large unperturbed sum.
 duplus=U*np.expm1(-(gamma-1)*np.log(vp/V0));duminus=U*np.expm1(-(gamma-1)*np.log(vm/V0))
 vprime=(vp-vm)/(2*eps);vsecond=(vp+vm-2*V0)/eps**2
 stiff=float(np.sum(gamma*P/V0*vprime**2));geom=float(-np.sum(P*vsecond))
 out['samples'].append({'epsilon_over_gap':fraction,'energy_curvature':float(np.sum(duplus+duminus)/eps**2),
 'compressibility_curvature':stiff,'pressure_geometry_curvature':geom,
 'volume_first_sum':float(vprime.sum()),'volume_second_sum':float(vsecond.sum()),
 'energy_first_derivative':float(np.sum(duplus-duminus)/(2*eps))})
root=base/('Codex_pair_thermal_curvature_diagonal_20261018' if args.diagonal else ('Codex_pair_thermal_curvature_tangent_20261018' if args.tangent else 'Codex_pair_thermal_curvature_20261018'));root.mkdir(exist_ok=False)
(root/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
