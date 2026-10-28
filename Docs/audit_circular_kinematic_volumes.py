"""Prescribed exact orbits, independent Voronoi geometry; no force evolution."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.spatial import Voronoi
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
base=Path('/gpfs/kjhan/LagForce');a,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1);om=np.sqrt(gravity_coefficient(r));active=r<2.5
sel=r<1.9;ids=np.where(sel)[0];gamma=float(np.float32(1.666667));Vref=None
masks={'inner':r<.5,'edge':(r>=.5)&(r<.6),'core':(r>=.6)&(r<1.9)}
out={'scope':'Prescribed exact circular orbits for r<2.5, frozen exterior; independent unweighted Voronoi geometry. No hydrodynamics, no actual steering, no physical stability pass. Pressure proxy is fixed-mass isentropic P/P0=(V0/V)^gamma.', 'frames':[]}
for time in np.arange(0.,21.,2.):
 theta=om*time;z=x.copy();z[active,0]=(x[:,0]*np.cos(theta)-x[:,1]*np.sin(theta))[active];z[active,1]=(x[:,0]*np.sin(theta)+x[:,1]*np.cos(theta))[active]
 mesh=Voronoi(z);V=np.full(len(a),np.nan)
 for i in ids:
  region=mesh.regions[mesh.point_region[i]];assert region and min(region)>=0
  q=mesh.vertices[region]-z[i];V[i]=abs(np.dot(q[:,0],np.roll(q[:,1],-1))-np.dot(q[:,1],np.roll(q[:,0],-1)))/2
 assert np.all(V[sel]>0)
 if Vref is None:
  Vref=V.copy();out['initial_volume_relative_max_error']=float(np.max(abs(V[sel]/a['volume'][sel]-1)))
 ratio=V/Vref;pp=ratio**(-gamma)
 frame={'time':float(time),'regions':{}}
 for name,mask in masks.items():
  rms=lambda v:float(np.sqrt(np.average(v[mask]**2,weights=a['mass'][mask])))
  frame['regions'][name]={'volume_relative_RMS':rms(ratio-1),'isentropic_pressure_relative_RMS':rms(pp-1),'volume_ratio_min':float(ratio[mask].min()),'volume_ratio_max':float(ratio[mask].max())}
 out['frames'].append(frame)
p=base/'Codex_circular_kinematic_volumes_20261018';p.mkdir(exist_ok=False)
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'initial_volume_relative_max_error':out['initial_volume_relative_max_error'],'final':out['frames'][-1],
 'maximum_sampled_pressure_RMS':{name:max(f['regions'][name]['isentropic_pressure_relative_RMS'] for f in out['frames']) for name in masks}},indent=2))
