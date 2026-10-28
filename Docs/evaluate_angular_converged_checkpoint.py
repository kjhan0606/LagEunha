"""Profile one existing quadratic RHS evaluation; no time evolution or solver changes."""
import cProfile, hashlib, importlib.util, json, os, pstats, sys, time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
pair=False
shared=True
sys.path[:0]=[str(base/'Codex_quadratic_thermal_audit_modules_20261018'),
             str(base/'Codex_quadratic_split_modules_20261018'),
             str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
if shared:sys.path.insert(0,str(base/'Codex_angular_faces_modules_20261018'))
import moment_ale_evolve as model
import fixedmass_hllc
from kinetic_factor_rhs import rhs
from quadratic_upwind_face_nodes import QuadraticUpwindNodes
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
os.environ.update(LAGFORCE_DENSITY_QUADRATURE_ORDER='12',LAGFORCE_QUADRATIC_VELOCITY='1',
    LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',
    LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',
    LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',
    LAGFORCE_VECTOR_KINETIC_ROOT='0')
os.environ['LAGFORCE_SHARED_FACE_VELOCITY']='1' if shared else '0'
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
tag=sys.argv[1]
assert tag in ('coarse','fine')
source=base/('Codex_quadratic_extended_dt1e4_20261018/step0003.npz' if tag=='coarse' else 'Codex_quadratic_extended_dt5e5_20261018/step0006.npz')
saved=np.load(source);y=saved['state'];mass=saved['mass'];ids=saved['ids'];L=saved['factors'];t=float(saved['time'])
ic0=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
frozen=np.linalg.norm(ic0['state'][:-4].reshape(-1,20)[:,:2],axis=1)>=.60
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
model.geometry=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
model.boundary_velocity=lambda point:np.zeros(2)
provider=QuadraticUpwindNodes(y[:-4].reshape(-1,20)[:,:2],mass,L,saved['weight_rate_seed'],8,model.boundary_velocity)

root=base/f'Codex_angular_converged_{tag}_t4e4_20261018';root.mkdir(exist_ok=False)
assert abs(t-4e-4)<1e-14
os.environ['LAGFORCE_ANGULAR_FACE_CELLS']=','.join(str(i) for i,k in enumerate(ids) if k in (956,941,864,924))
start=time.monotonic()
prior=base/f'Codex_angular_faces_{tag}_t4e4_20261018'
prior_report=json.loads((prior/'report.json').read_text())
assert prior_report['sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
prior_rates=np.load(prior/'rates.npz')
assert np.array_equal(prior_rates['ids'],ids) and np.array_equal(prior_rates['mass'],mass)
z=prior_rates['rates'][:-4].reshape(-1,20)[:,2].copy()
changes=[float(abs(z-saved['weight_rate_seed']).max())]
for iteration in range(1,10):
    provider=QuadraticUpwindNodes(y[:-4].reshape(-1,20)[:,:2],mass,L,z,8,model.boundary_velocity)
    d,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=provider)
    newz=d[:-4].reshape(-1,20)[:,2].copy()
    changes.append(float(abs(newz-z).max()));z=newz
    event=dict(iteration=iteration,weight_rate_change=changes[-1])
    with (root/'events.jsonl').open('a') as f:f.write(json.dumps(event)+'\n')
    if changes[-1]<1e-11:break
else:raise RuntimeError('weight-rate outer iteration failed')
keys=['internal_angular_rate','pressure_internal_torque','transport_internal_angular_rate']
assert np.isfinite(d).all() and all(np.isfinite(diag[k]).all() for k in keys)
np.savez(root/'rates.npz',ids=ids,mass=mass,rates=d,**{k:diag[k] for k in keys})
report=dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),time=t,
 weight_rate_changes=changes,seconds=time.monotonic()-start,mass_flux=float(diag['mass_flux']),
 face_partition={int(ids[i]):{int(ids[j]) if j>=0 else -1:v for j,v in faces.items()} for i,faces in diag['angular_face_partition'].items()},
 selected=[dict(id=int(ids[i]),**{k:float(diag[k][i]) for k in keys}) for i in range(len(ids)) if ids[i] in (956,941,864,924)],
 limits='Hydro RHS after outer weight/split fixed-point change <1e-11; first evaluation reused with source SHA verification. No gravity or integrated attribution. No evolution.')
for i,faces in diag['angular_face_partition'].items():
    actual=sum(f['total'] for f in faces.values())
    assert abs(actual-diag['transport_internal_angular_rate'][i])<1e-14
    assert abs(sum(f['mean']+f['residual']-f['total'] for f in faces.values()))<1e-14
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
