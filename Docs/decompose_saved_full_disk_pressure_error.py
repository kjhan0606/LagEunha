"""Saved-state pressure budget identity; no provider or evolution call.

Separates moving volume, raw-moment thermal transport, Pi factor rebuild,
and stored energy-update rounding at the common endpoint t=.001.
"""
import argparse,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);base=Path('/gpfs/kjhan/LagForce');ld=np.longdouble
initial=np.load(base/'Codex_full_curved_kepler_IC_20261019/initial.npz');mass=initial['mass'];N=len(mass);a0=initial['state'][:-4].reshape(N,20).astype(ld);V0=np.load(base/'Codex_full_curved_disk_mesh_20261019/mesh.npz')['moments'][:,0].astype(ld);physics=json.loads((base/'Codex_full_curved_kepler_IC_20261019/report.json').read_text())['physics'];gm=ld(physics['gamma'])-1;P0=ld(physics['P0']);Pi0=a0[:,15:19].reshape(N,2,2);U0=a0[:,19]-.5*np.trace(Pi0,axis1=1,axis2=2);r0=np.linalg.norm(initial['state'][:-4].reshape(N,20)[:,:2],axis=1);reports={}
def stats(v,mask):
 x=v[mask];w=mass[mask].astype(ld);return dict(RMS=float(np.sqrt(np.average(x*x,weights=w))),maximum_absolute=float(abs(x).max()),signed_mass_mean=float(np.average(x,weights=w)))
for name,rates in [('coarse',[(.001,'Codex_full_disk_coarse_midpoint_stage_20261019')]),('fine_final',[(.0005,'Codex_full_disk_fine_midpoint_stage_20261019'),(.0005,'Codex_full_disk_fine_second_midpoint_stage_20261019')])]:
 saved=np.load(base/f'Codex_full_disk_{name}_accepted_20261019/accepted.npz');a=saved['state'][:-4].reshape(N,20).astype(ld);Pi=a[:,15:19].reshape(N,2,2);V=np.load(base/f'Codex_full_disk_{name}_endpoint_geometry_20261019/mesh.npz')['moments'][:,0].astype(ld);sumE=np.zeros(N,dtype=ld);sumPi=np.zeros((N,2,2),dtype=ld);sumGthermal=np.zeros(N,dtype=ld)
 for h,dirname in rates:
  stage=np.load(base/dirname/'stage0.npz')['stage'].astype(ld);h=ld(h);sumE+=h*stage[:,17];sumPi+=h*stage[:,13:17].reshape(N,2,2);sumGthermal+=h*(stage[:,47]-.5*(stage[:,43]+stage[:,46]))
 scale=gm/(P0*V);deltaPi=Pi-Pi0-sumPi
 parts=dict(volume_and_initial_pressure=(gm*U0/(P0*V)-1),raw_thermal_transport=scale*(sumE-.5*np.trace(sumPi,axis1=1,axis2=2)),Pi_factor_rebuild=-.5*scale*np.trace(deltaPi,axis1=1,axis2=2),stored_energy_update_rounding=scale*(a[:,19]-a0[:,19]-sumE))
 actual=gm*(a[:,19]-.5*np.trace(Pi,axis1=1,axis2=2))/(P0*V)-1;reconstructed=sum(parts.values());residual=actual-reconstructed
 primitive=1+scale*(abs(a[:,19])+abs(a0[:,19])+abs(sumE)+.5*np.sum(abs(Pi)+abs(Pi0)+abs(sumPi),axis=(1,2)));bound=128*np.finfo(ld).eps*primitive;assert np.all(abs(residual)<=bound)
 regions={}
 for region,mask in dict(all=np.ones(N,bool),inner=r0<.4,disk_initial_radius_0p5_to_2=(r0>.5)&(r0<2)).items():
  regions[region]=dict(cells=int(mask.sum()),actual_pressure_error=stats(actual,mask),parts={k:stats(v,mask) for k,v in parts.items()},gravity_thermal_subset_of_raw_transport=stats(scale*sumGthermal,mask))
 np.savez(args.output/(name+'.npz'),ids=initial['ids'],actual_relative_pressure_error=actual,identity_residual=residual,identity_bound=bound,gravity_thermal_subset=scale*sumGthermal,**parts)
 reports[name]=dict(time=float(saved['time']),regions=regions,maximum_identity_residual=float(abs(residual).max()),maximum_identity_bound_fraction=float((abs(residual)/bound).max()))
result=dict(status='completed_saved_pressure_budget_identity',scope=__doc__,paths=reports,provider_calls=0,evolution_steps=0,limits='Algebraic decomposition on saved rates/states/actual volumes. Gravity thermal is a subset, not an extra summand. Component RMS values need not add because components cancel. No term is changed, removed, or reinjected; this diagnostic does not by itself identify a unique physical cause.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
