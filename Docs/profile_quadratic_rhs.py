"""Profile one existing quadratic RHS evaluation; no time evolution or solver changes."""
import cProfile, hashlib, importlib.util, json, os, pstats, sys, time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
pair='--verify-pair' in sys.argv
shared='--shared' in sys.argv or pair
sys.path[:0]=[str(base/'Codex_quadratic_thermal_audit_modules_20261018'),
             str(base/'Codex_quadratic_split_modules_20261018'),
             str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
if shared:sys.path.insert(0,str(base/'Codex_shared_face_velocity_modules_20261018'))
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
source=base/'Codex_quadratic_r060_dt5e5_20261018/step0002.npz'
saved=np.load(source);y=saved['state'];mass=saved['mass'];ids=saved['ids'];L=saved['factors'];t=float(saved['time'])
ic0=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
frozen=np.linalg.norm(ic0['state'][:-4].reshape(-1,20)[:,:2],axis=1)>=.60
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
model.geometry=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
model.boundary_velocity=lambda point:np.zeros(2)
provider=QuadraticUpwindNodes(y[:-4].reshape(-1,20)[:,:2],mass,L,saved['weight_rate_seed'],8,model.boundary_velocity)
root=base/('Codex_shared_face_velocity_full_audit_20261018' if pair else ('Codex_shared_face_velocity_profile_20261018' if shared else 'Codex_quadratic_rhs_profile_20261018'));root.mkdir(exist_ok=False)
if pair:
    from factor_moment_evolve import advance_factors_sylvester
    values=[];seconds=[]
    for enabled in ('0','1'):
        os.environ['LAGFORCE_SHARED_FACE_VELOCITY']=enabled
        nodes=QuadraticUpwindNodes(y[:-4].reshape(-1,20)[:,:2],mass,L,saved['weight_rate_seed'],8,model.boundary_velocity)
        start=time.monotonic()
        values.append(rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=nodes))
        seconds.append(time.monotonic()-start)
    checked=[]
    def compare(a,b,path):
        if isinstance(a,dict):
            assert a.keys()==b.keys(),path
            for k in a:compare(a[k],b[k],path+'.'+k)
        elif isinstance(a,(list,tuple)):
            assert len(a)==len(b),path
            for k,(aa,bb) in enumerate(zip(a,b)):compare(aa,bb,path+f'[{k}]')
        else:
            assert np.array_equal(a,b),path
            checked.append(path)
    compare(values[0],values[1],'rhs_and_diagnostics')
    for dt in (2.5e-5,5e-5,1e-4):
        updated=[np.asarray(advance_factors_sylvester(L,v[1],dt)) for v in values]
        assert np.array_equal(*updated)
    report=dict(scope='Both paths in one process on the actual checkpoint; all RHS and diagnostic leaves plus factor advances.',
        source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        bitwise_equal=True,checked_leaves=len(checked),diagnostic_keys=sorted(values[0][1]),
        baseline_seconds=seconds[0],shared_seconds=seconds[1],factor_step_sizes=[2.5e-5,5e-5,1e-4],
        limits='One hydro RHS state and frozen-coefficient factor updates; not complete midpoint evolution, gravity or production CUDA/MPI validation.')
    (root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));sys.exit(0)
prof=cProfile.Profile();start=time.monotonic();prof.enable()
d,diag=rhs(t,y,mass,covariance_factors=L,frozen_mask=frozen,face_order=8,face_node_provider=provider)
prof.disable();wall=time.monotonic()-start;prof.dump_stats(str(root/'rhs.prof'))
assert np.isfinite(d).all()
np.savez(root/'rhs.npz',rates=d,ids=ids)
stats=pstats.Stats(prof)
entries=[dict(file=k[0],line=k[1],function=k[2],primitive_calls=v[0],calls=v[1],self_seconds=v[2],cumulative_seconds=v[3]) for k,v in stats.stats.items()]
report=dict(scope=__doc__,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    shared_face_velocity=shared,
    wall_seconds=wall,mass_flux=float(diag['mass_flux']),min_pressure=float(diag['min_pressure']),
    top_cumulative=sorted(entries,key=lambda v:v['cumulative_seconds'],reverse=True)[:35],
    top_self=sorted(entries,key=lambda v:v['self_seconds'],reverse=True)[:25],
    limits='Profiler overhead included; one hydro RHS with supplied split seed, not the outer fixed point or gravity. No physical validation claim.')
(root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
