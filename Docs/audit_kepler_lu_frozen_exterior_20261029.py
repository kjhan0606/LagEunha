"""Audit one historical n64 moment step after the frozen-exterior fix."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

root=Path(sys.argv[1]);docs=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=root/'initial_state.json';endpoint=root/'one_step_checkpoint.json'
a=json.loads(initial.read_text());b=json.loads(endpoint.read_text());report=json.loads((root/'one_step_report.json').read_text())
n=len(a['mass']);y0=np.array(a['state']);y1=np.array(b['state']);p0=y0[:-4].reshape(n,20);p1=y1[:-4].reshape(n,20)
S=lambda x:np.sum(x,dtype=np.longdouble)
gas_delta=float(S(p1[:,19])-S(p0[:,19]));ledger=y1[-4:]-y0[-4:]
energy_residual=gas_delta-ledger[2]-report['gravity_work']
def Lcenter(p):return S(p[:,12]-p[:,13]-3*p[:,10]+3*p[:,9])
angular_delta=float(Lcenter(p1)-Lcenter(p0));angular_ledger=float(ledger[3]-3*ledger[1]+3*ledger[0])
thermal=p1[:,19]-.5*(p1[:,15]+p1[:,18])
checks={'source_hash':sha(initial)=='4722b8054be55b9a03e6cddd7039034c5c1bd34d3dc4fd6b548e93af0b9d25ab',
        'one_step_completed':report['completed'] and b['t']==1e-5,
        'fixed_mass':np.array_equal(a['mass'],b['mass']),
        'finite_state':bool(np.isfinite(y1).all()),'positive_thermal':bool(np.isfinite(thermal).all() and np.all(thermal>0)),
        'finite_stage_mass_residual':np.isfinite(report['stage_mass_flux']) and report['stage_mass_flux']<2e-13,
        'active_generator_velocity_no_steering':report['generator_velocity_defect']==0,
        'gas_energy_boundary_and_gravity_budget':abs(energy_residual)<1e-10,
        'central_angular_momentum_boundary_budget':abs(angular_delta-angular_ledger)<1e-10}
result={'passed':all(checks.values()),'checks':{k:bool(v) for k,v in checks.items()},
        'cells':n,'frozen_cells':report['frozen_cells'],'dt':b['t'],'min_thermal_energy':float(thermal.min()),
        'gas_energy_delta':gas_delta,'boundary_energy_exchange':float(ledger[2]),'gravity_work':report['gravity_work'],
        'gas_energy_budget_residual':energy_residual,'central_angular_momentum_delta':angular_delta,
        'central_boundary_angular_exchange':angular_ledger,'central_angular_budget_residual':angular_delta-angular_ledger,
        'wall_seconds_shared_host':report['wall_seconds'],'source_hashes':{str(p):sha(p) for p in (initial,endpoint,root/'kinetic_factor_rhs.py',root/'resume_one_step.py',Path(__file__))},
        'scope':'One historical 4084-cell Python moment step from saved initial state with prescribed frozen outer cells, not the current no-steering production wall setup or native CUDA/MPI. Energy and angular budgets are one-step checks, not temporal convergence, endpoint geometry realizability, physical accuracy, performance or t>6 stability.'}
out=docs/'KEPLER_LU_FROZEN_EXTERIOR_STEP_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
