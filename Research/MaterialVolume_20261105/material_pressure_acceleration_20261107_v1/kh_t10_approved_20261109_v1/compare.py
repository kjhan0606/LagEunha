"""Saved-state KH density and conservation report; no evolution or filtering."""
from pathlib import Path
import ast,json,hashlib,sys
import numpy as np,h5py
from scipy.spatial import cKDTree
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;D=O.parent;nx=int(sys.argv[1]);target=float(sys.argv[2]);P=O/f'n{nx}';Q=P/'comparisons'/f't{target:g}';Q.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
res=json.loads((P/'native'/f'endpoint_{target:g}.json').read_text());assert res['passed'];npth=Path(res['checkpoint']).parent/'state_ids.bin';a=np.fromfile(npth).reshape(-1,13);a=a[a[:,12].argsort()]
gpath=None
for path in sorted((P/'gizmo/output').glob('snap_*.hdf5')):
 with h5py.File(path) as f:
  if abs(float(f['Header'].attrs['Time'])-target)<1e-10:gpath=path
assert gpath is not None
with h5py.File(gpath) as f:
 g=f['PartType0'];order=np.argsort(g['ParticleIDs'][:]);x=g['Coordinates'][:][order,:2];v=g['Velocities'][:][order,:2];m=g['Masses'][:][order];u=g['InternalEnergy'][:][order];rho=g['Density'][:][order];ids=g['ParticleIDs'][:][order];actual_time=float(f['Header'].attrs['Time'])
assert len(a)==2*nx*nx and np.array_equal(a[:,12]+1,ids) and np.array_equal(a[:,6],m)
assert np.isfinite(np.c_[x,v,m,u,rho]).all() and (u>0).all() and (rho>0).all()
a0=np.fromfile(P/'input/t0_ids.bin').reshape(-1,13);a0=a0[a0[:,12].argsort()]
def quantities(v,m,u):
 m=m.astype(np.longdouble);v=v.astype(np.longdouble);u=u.astype(np.longdouble);kin=np.sum(.5*m*np.sum(v*v,axis=1));internal=np.sum(m*u)
 return dict(mass=m.sum(),momentum_x=np.sum(m*v[:,0]),momentum_y=np.sum(m*v[:,1]),kinetic=kin,internal=internal,energy=kin+internal)
initial=quantities(a0[:,2:4],a0[:,6],a0[:,7]/a0[:,6]);nq=quantities(a[:,2:4],a[:,6],a[:,7]/a[:,6]);gq=quantities(v,m,u)
with h5py.File(P/'input/snap_001.hdf5') as f:
 g=f['PartType0'];restart=quantities(g['Velocities'][:,:2],g['Masses'][:],g['InternalEnergy'][:])
report=dict(time=target,gizmo_actual_snapshot_time=actual_time,cells=len(a),shape=[nx,2*nx],mass_equal_per_cell=True,native_checkpoint=res['checkpoint'],native_checkpoint_sha256=res['checkpoint_sha256'],gizmo_snapshot=str(gpath),gizmo_snapshot_sha256=sha(gpath),native=res,initial={k:float(v) for k,v in initial.items()},gizmo={k:float(v) for k,v in gq.items()},gizmo_energy_relative_from_t0=float((gq['energy']-initial['energy'])/initial['energy']),gizmo_energy_relative_from_restart=float((gq['energy']-restart['energy'])/restart['energy']),gizmo_momentum_change=[float(gq['momentum_x']-initial['momentum_x']),float(gq['momentum_y']-initial['momentum_y'])],gizmo_conservation_scope='HDF5 output-time predicted variables, not synchronized raw restart variables. Snapshot/time-bin convention is included in these diagnostics. No matched-accuracy speed claim.',initialization_caveat='128 GIZMO restarts from synchronized exported t=1.35 state to allow 4 MPI ranks; density is recomputed. Original t=0 GIZMO pressure offset about0.0581% remains in inherited history.')
# Reuse the established periodic Laguerre ownership certificate, unchanged.
source=D/'common_closed_volume_20261109_v1/full_kh128/assess_checkpoint.py';tree=ast.parse(source.read_text());func=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='raster'];assert len(func)==1;ns=dict(np=np,cKDTree=cKDTree);exec(compile(ast.Module(body=func,type_ignores=[]),str(source),'exec'),ns)
field,margin=ns['raster'](a,512);nfield=field[:,:,0]
xx,yy=np.meshgrid((np.arange(512)+.5)/512,2*(np.arange(1024)+.5)/1024);_,owner=cKDTree(x%[1,2],boxsize=[1,2]).query(np.c_[xx.ravel(),yy.ravel()]);gfield=rho[owner].reshape(1024,512)
plt.rcParams.update({'font.size':12});fig=plt.figure(figsize=(7.5,6.5));cmin,cmax=.98,2.03
for i,(f,label,note) in enumerate(zip([nfield,gfield],['LagForce','GIZMO MFM'],['Laguerre cell density','Nearest density samples'])):
 ax=fig.add_axes([.1+i*.35,.10,.35,.81]);im=ax.imshow(f,origin='lower',extent=(0,1,0,2),interpolation='nearest',cmap='viridis',vmin=cmin,vmax=cmax,aspect='equal');ax.set(xlabel='x',xticks=[.2,.5,.8]);ax.tick_params(direction='in',top=True,right=True)
 if i==0:ax.set_ylabel('y')
 else:ax.tick_params(labelleft=False)
 ax.text(.035,.975,f'{label}\n{nx} × {2*nx}, t = {target:g}',transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.85,edgecolor='none'));ax.text(.035,.02,note,fontsize=10,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
cax=fig.add_axes([.84,.13,.025,.75]);fig.colorbar(im,cax=cax,label=r'$\rho$');fig.savefig(Q/'density.png',dpi=180,bbox_inches='tight');fig.savefig(Q/'density.pdf',bbox_inches='tight');plt.close(fig)
report['display']=dict(native='Certified actual power-cell ownership, no smoothing',gizmo='Nearest stored density samples; not MFM volume boundaries',color_limits=[cmin,cmax],native_owner_margin=margin,native_density_range=[float(nfield.min()),float(nfield.max())],gizmo_density_range=[float(rho.min()),float(rho.max())],native_pixel_clipped_fraction=float(np.mean((nfield<cmin)|(nfield>cmax))),gizmo_pixel_clipped_fraction=float(np.mean((gfield<cmin)|(gfield>cmax))),note='Different density representations: do not interpret pixel differences as a physical error norm.')
(Q/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['time','cells','gizmo_energy_relative_from_t0','gizmo_momentum_change']}),flush=True)
