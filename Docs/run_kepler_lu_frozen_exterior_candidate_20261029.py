import hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Kepler_LU_FrozenBoundaryFix_20261029')
old=Path('/gpfs/kjhan/LagForce/Codex_LU_kepler_preflight_n64_20261017_093046')
sys.path.insert(0,str(old));sys.path.insert(0,str(root))
from physical_field_accuracy_audit import environment
cfg=json.loads((root/'configuration.json').read_text());os.environ.update(environment(cfg));os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING='1',LAGFORCE_PRESSURE_RIEMANN='exact_euler_star')
from check_kepler_moment_preflight import install
install()
from kepler_disk_reference import frozen,CENTER
from kepler_spline_gravity import acceleration
from gravity_coupled_2d import step
path=root/'initial_state.json';raw=path.read_bytes()
if hashlib.sha256(raw).hexdigest()!='4722b8054be55b9a03e6cddd7039034c5c1bd34d3dc4fd6b548e93af0b9d25ab':raise ValueError('changed initial')
state=json.loads(raw);y=np.asarray(state['state']);mass=np.asarray(state['mass']);L=[np.asarray(v) for v in state['factors']];mask=frozen(0.,y[:-4].reshape(-1,20)[:,:2]);dt=cfg['dt'];start=time.monotonic()
out,nextL,diag=step(0.,y,mass,L,dt,lambda p:acceleration(p,center=CENTER,rc=.5),frozen=frozen)
arr=out[:-4].reshape(-1,20);before=y[:-4].reshape(-1,20)
np.testing.assert_array_equal(arr[mask],before[mask])
for i in np.flatnonzero(mask):np.testing.assert_array_equal(nextL[i],L[i])
if not np.all(np.isfinite(out)) or diag['stage_mass_flux']>2e-13 or diag['generator_velocity_defect']!=0:raise ValueError('step contract')
checkpoint={'t':dt,'mass':mass.tolist(),'state':out.tolist(),'factors':[v.tolist() for v in nextL]}
with (root/'one_step_checkpoint.json').open('x') as f:json.dump(checkpoint,f)
report={'completed':True,'initial_sha256':hashlib.sha256(raw).hexdigest(),'cells':len(mass),'frozen_cells':int(mask.sum()),'dt':dt,'stage_mass_flux':float(diag['stage_mass_flux']),'generator_velocity_defect':float(diag['generator_velocity_defect']),'gravity_work':float(diag['gravity_work']),'max_thermal_kick_change':float(diag['max_thermal_kick_change']),'wall_seconds':time.monotonic()-start,'scope':'One historical n64 Python moment step with frozen outer cells; not the current native GPU solver, no long-time or full physical accuracy claim.'}
with (root/'one_step_report.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(report),flush=True)
