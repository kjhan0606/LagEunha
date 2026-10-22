from pathlib import Path
import json,sys,time
import numpy as np
from scipy.spatial import Voronoi
from audit_snapshot import read_snapshot
root=Path(__file__).resolve().parent
out={}
for label,left,right in [('t6','final_cpu_late','final_gpu_late'),('t8','final_cpu_t8','final_gpu_t8')]:
 def state(d):
  f=next(p for p in (root/d).glob('khout.*.dat') if not p.is_symlink());p,m=read_snapshot(f)
  ids=(p['u4if']&np.uint64((1<<56)-1)).astype(int);return p[np.argsort(ids)],m
 p,m=state(left);q,_=state(right);xy=np.column_stack([p['x']%m['Lx'],p['y']%m['Ly']]);pad=3.
 blocks=[xy]
 for dx in (-m['Lx'],0,m['Lx']):
  for dy in (-m['Ly'],0,m['Ly']):
   if dx==0 and dy==0:continue
   z=xy+np.array([dx,dy]);sel=(z[:,0]>=-pad)&(z[:,0]<=m['Lx']+pad)&(z[:,1]>=-pad)&(z[:,1]<=m['Ly']+pad);blocks.append(z[sel])
 start=time.monotonic();vor=Voronoi(np.vstack(blocks));area=np.empty(len(p));rad=np.empty(len(p))
 for i in range(len(p)):
  reg=vor.regions[vor.point_region[i]]
  if not reg or -1 in reg:raise ValueError('Unbounded base cell')
  v=vor.vertices[reg]-xy[i]
  area[i]=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))-np.dot(v[:,1],np.roll(v[:,0],-1)));rad[i]=np.sqrt(np.max(np.sum(v*v,axis=1)))
 clearance=pad+np.minimum.reduce([xy[:,0],xy[:,1],m['Lx']-xy[:,0],m['Ly']-xy[:,1]])
 if not np.all(2*rad<clearance):raise ValueError('Periodic tiling completeness not proved')
 changed=p['volume']!=q['volume'];r=np.hypot(p['x']-12,p['y']-12)
 out[label]={'time':m['t'],'n':len(p),'tiled_points':sum(len(b) for b in blocks),'all_cells_periodic_complete':True,'analysis_seconds':time.monotonic()-start,'area_sum':float(area.sum()),'expected_area':m['Lx']*m['Ly'],'cpu_max_rel_error':float(np.max(abs(p['volume']-area)/area)),'gpu_max_rel_error':float(np.max(abs(q['volume']-area)/area)),'cpu_cells_rel_error_gt_1e-6':int(np.sum(abs(p['volume']-area)/area>1e-6)),'gpu_cells_rel_error_gt_1e-6':int(np.sum(abs(q['volume']-area)/area>1e-6)),'changed_cells':int(changed.sum()),'changed_radius_min':float(r[changed].min()),'changed_radius_max':float(r[changed].max()),'min_area':float(area.min()),'max_area':float(area.max())}
 print(label,out[label])
(root/'final_volume_gpu/periodic_area_audit.json').write_text(json.dumps(out,indent=2)+'\n')
