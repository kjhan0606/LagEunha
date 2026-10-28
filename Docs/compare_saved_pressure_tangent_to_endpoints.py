"""Check whether initial spatial pressure drift predicts saved endpoint errors."""
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
base=Path('/gpfs/kjhan/LagForce');ld=np.longdouble
ic=np.load(base/'Codex_full_curved_kepler_IC_20261019/initial.npz')
physics=json.loads((base/'Codex_full_curved_kepler_IC_20261019/report.json').read_text())['physics']
a=ic['state'][:-4].reshape(-1,20).astype(ld);m=ic['mass'].astype(ld)
V=np.load(base/'Codex_full_curved_disk_mesh_20261019/mesh.npz')['moments'][:,0].astype(ld)
initial_error=(ld(physics['gamma'])-1)*(a[:,19]-.5*(a[:,15]+a[:,18]))/(ld(physics['P0'])*V)-1
tangent=np.load(base/'Codex_full_disk_saved_pressure_tangent_20261019/initial.npz')['factor_relative_pressure_rate']
r0=np.linalg.norm(a[:,:2],axis=1)
regions=dict(all=np.ones(len(m),bool),inner=r0<.4,disk_initial_radius_0p5_to_2=(r0>.5)&(r0<2))
reports={}
for label in ('coarse','fine_final'):
    actual=np.load(base/f'Codex_full_disk_saved_pressure_budget_20261019/{label}.npz')['actual_relative_pressure_error']
    predicted=initial_error+ld(.001)*tangent
    result={}
    for name,mask in regions.items():
        w=m[mask];av=actual[mask];pv=predicted[mask];dv=av-pv
        rms=lambda v:float(np.sqrt(np.average(v*v,weights=w)))
        result[name]=dict(actual_RMS=rms(av),predicted_RMS=rms(pv),remaining_RMS=rms(dv),
            remaining_to_actual_RMS_ratio=rms(dv)/rms(av),
            weighted_uncentered_cosine=float(np.sum(w*av*pv)/np.sqrt(np.sum(w*av*av)*np.sum(w*pv*pv))))
    reports[label]=result
report=dict(status='completed_saved_pressure_tangent_endpoint_alignment',time=.001,paths=reports,
    provider_calls=0,evolution_steps=0,limits='First-order local prediction from initial saved tangent; agreement diagnoses early spatial equilibrium bias but does not identify a unique closure term or predict long-time drift.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
