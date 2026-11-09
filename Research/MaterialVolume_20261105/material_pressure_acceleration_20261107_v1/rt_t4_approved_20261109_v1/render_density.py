from pathlib import Path
import sys,json,time,hashlib,shutil
import numpy as np
R=Path('/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1')
D=R/'Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1'
O=D/'rt_t4_approved_20261109_v1';F=O/'figures'
r=json.loads((O/'result.json').read_text());assert r['passed'] and r['saved_time']==4.0
assert hashlib.sha256(Path(r['checkpoint']).read_bytes()).hexdigest()==r['checkpoint_sha256']
F.mkdir(exist_ok=False);shutil.copy2(__file__,F/'run.py')
sys.path.insert(0,str(D));from rt_power_polygons import reconstruct
a=np.fromfile(O/'state_ids.bin').reshape(-1,13);a=a[a[:,12].argsort()]
tick=time.monotonic();print('Reconstructing final saved power cells',flush=True)
polys,owners,g=reconstruct(a);print('Geometry certified',json.dumps(g),flush=True)
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
rho=a[:,6]/a[:,8]
plt.rcParams.update({'font.size':12})
fig=plt.figure(figsize=(4.8,7.0));ax=fig.add_axes([.16,.11,.63,.84])
pc=PolyCollection(polys,array=rho[owners],cmap='viridis',edgecolors='none',antialiaseds=False)
pc.set_clim(.8,2.8);ax.add_collection(pc)
ax.set(xlim=(0,.5),ylim=(0,1),xlabel='x',ylabel='y',aspect='equal')
ax.text(.04,.975,'RT 128 × 256\nt = 4',transform=ax.transAxes,va='top',fontsize=12,bbox=dict(facecolor='white',alpha=.88,edgecolor='none'))
cax=fig.add_axes([.84,.16,.035,.72]);fig.colorbar(pc,cax=cax,label=r'$\rho$',extend='both')
fig.savefig(F/'density_t4.png',dpi=240,bbox_inches='tight');fig.savefig(F/'density_t4.pdf',bbox_inches='tight');plt.close(fig)
g.update(time=4.0,shape=[128,256],source_checkpoint=r['checkpoint'],checkpoint_sha256=r['checkpoint_sha256'],state_sha256=hashlib.sha256((O/'state_ids.bin').read_bytes()).hexdigest(),geometry_helper_sha256=hashlib.sha256((D/'rt_power_polygons.py').read_bytes()).hexdigest(),density_min=float(rho.min()),density_max=float(rho.max()),display_limits=[.8,2.8],cells_below_display=int(np.sum(rho<.8)),cells_above_display=int(np.sum(rho>2.8)),drawing='Actual saved Laguerre power-cell polygons, periodic x pieces clipped to box, no smoothing or field interpolation.',wall_seconds=time.monotonic()-tick,reference_scope='Physical reflecting walls y=0,1; this is not a validated equal-BC comparison with the GIZMO diagnostic.')
(F/'report.json').write_text(json.dumps(g,indent=2)+'\n')
print('DENSITY_REPORT',json.dumps(g),flush=True)
