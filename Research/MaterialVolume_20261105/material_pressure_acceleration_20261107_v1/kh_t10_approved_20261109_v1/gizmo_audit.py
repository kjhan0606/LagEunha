from pathlib import Path
import json,hashlib,subprocess
import numpy as np,h5py
O=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def quantities(v,m,u):
 m=m.astype(np.longdouble);v=v.astype(np.longdouble);u=u.astype(np.longdouble);kin=(.5*m*(v*v).sum(axis=1)).sum();internal=(m*u).sum()
 return dict(mass=m.sum(),momentum_x=(m*v[:,0]).sum(),momentum_y=(m*v[:,1]).sum(),kinetic=kin,internal=internal,energy=kin+internal)
def audit(nx,snapshot):
 P=O/f'n{nx}';out=P/'gizmo/audits';out.mkdir(exist_ok=True);dest=out/(snapshot.stem+'.json')
 if dest.exists():return
 a=np.fromfile(P/'input/t0_ids.bin').reshape(-1,13);a=a[a[:,12].argsort()];initial=quantities(a[:,2:4],a[:,6],a[:,7]/a[:,6])
 with h5py.File(snapshot) as f:
  g=f['PartType0'];order=np.argsort(g['ParticleIDs'][:]);ids=g['ParticleIDs'][:][order];x=g['Coordinates'][:][order];v=g['Velocities'][:][order];m=g['Masses'][:][order];u=g['InternalEnergy'][:][order];rho=g['Density'][:][order];t=float(f['Header'].attrs['Time'])
 q=quantities(v,m,u);r=dict(time=t,cells=len(m),fixed_mass_exact=bool(np.array_equal(ids,a[:,12]+1) and np.array_equal(m,a[:,6])),finite=bool(np.isfinite(np.c_[x,v,m,u,rho]).all()),minimum_u=float(u.min()),density_range=[float(rho.min()),float(rho.max())],quantities={k:float(v) for k,v in q.items()},energy_relative_change=float((q['energy']-initial['energy'])/initial['energy']),momentum_change=[float(q['momentum_x']-initial['momentum_x']),float(q['momentum_y']-initial['momentum_y'])],source=str(snapshot),sha256=sha(snapshot),scope='Snapshot predicted quantities at output time; includes output convention and rank-change snapshot restart; not a synchronized energy certificate.')
 r['finite_mass_positive_pass']=r['fixed_mass_exact'] and r['finite'] and r['minimum_u']>0 and r['density_range'][0]>0
 tmp=dest.with_suffix('.tmp');tmp.write_text(json.dumps(r,indent=2));tmp.replace(dest)
 if not r['finite_mass_positive_pass']:raise RuntimeError(f'Invalid GIZMO snapshot {snapshot}')
def final(nx):
 P=O/f'n{nx}';out=P/'gizmo/audits';out.mkdir(exist_ok=True);dest=out/'synchronized_final.json'
 if dest.exists():return
 states=[];meta=[]
 for rank in range(4):
  src=P/f'gizmo/output/restartfiles/restart.{rank}';export=out/f'rank{rank}_final.bin'
  if not export.exists():
   r=subprocess.run([str(O/'bin/read_gizmo_rank'),str(src),str(export),'10'],capture_output=True,text=True);assert r.returncode==0,r.stderr;meta.append(json.loads(r.stdout))
  z=np.fromfile(export).reshape(-1,10);states.append(z)
 z=np.concatenate(states);z=z[z[:,7].argsort()];a=np.fromfile(P/'input/t0_ids.bin').reshape(-1,13);a=a[a[:,12].argsort()]
 assert np.array_equal(z[:,7],a[:,12]+1) and np.array_equal(z[:,4],a[:,6]) and np.isfinite(z).all() and (z[:,5:7]>0).all()
 initial=quantities(a[:,2:4],a[:,6],a[:,7]/a[:,6]);q=quantities(z[:,2:4],z[:,4],z[:,5]);report=dict(time=10.,rank_count=4,cells=len(z),fixed_mass_exact=True,energy_relative_change=float((q['energy']-initial['energy'])/initial['energy']),momentum_change=[float(q['momentum_x']-initial['momentum_x']),float(q['momentum_y']-initial['momentum_y'])],quantities={k:float(v) for k,v in q.items()},scope='Synchronized final restart at original integer timeline endpoint TimeMax=10; combined four rank exports.',rank_metadata=meta)
 dest.write_text(json.dumps(report,indent=2));z.tofile(out/'synchronized_final.bin')
