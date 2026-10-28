"""Physical diagnostics against the steady field on each actual moving cell.

Reference data are diagnostic only. No state reset, provider, or evolution.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
for name in ('state','geometry','reference','density','output'):parser.add_argument('--'+name,type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False);ld=np.longdouble
initial_path=Path('/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz');initial=np.load(initial_path);actual=np.load(args.state);ref=np.load(args.reference/'reference.npz');mesh=np.load(args.geometry/'mesh.npz');n=len(actual['ids'])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
state_sha=sha(args.state);assert str(ref['actual_state_sha256'])==state_sha
assert json.loads((args.reference/'report.json').read_text())['status']=='passed_steady_reference_on_moving_cells'
density_report=json.loads((args.density/'report.json').read_text());assert density_report['status']=='passed_full_stage_Psi_density_preparation' and density_report['state_sha256']==state_sha
for d in (initial,ref,mesh):assert np.array_equal(d['ids'],actual['ids'])
assert np.array_equal(actual['mass'],initial['mass'])
a=actual['state'][:-4].reshape(n,20);assert np.array_equal(a[:,:2],mesh['sites']) and np.array_equal(a[:,2],mesh['weights'])
area=mesh['moments'][:,0];gc=mesh['moments'][:,1:3]/area[:,None];radius=np.linalg.norm(gc,axis=1);er=np.zeros_like(gc);defined=radius>0;er[defined]=gc[defined]/radius[defined,None];et=np.column_stack([-er[:,1],er[:,0]])
physics=json.loads((args.reference/'report.json').read_text())['physics'];gamma=physics['gamma'];P0=physics['P0']
def decode(saved,reference=False):
 m=saved['mass'];a=saved['state'][:-4].reshape(n,20);r=np.column_stack([m,a[:,3:]]);q=r[:,1:3];p=r[:,7:9];D=r[:,9:13].reshape(n,2,2);Q=r[:,3:7].reshape(n,2,2)
 S=Q-q[:,:,None]*q[:,None,:]/m[:,None,None];C=D-q[:,:,None]*p[:,None,:]/m[:,None,None];B=np.linalg.solve(S,C).transpose(0,2,1);L=saved['factors'];R=L@L.transpose(0,2,1);Pi=r[:,13:17].reshape(n,2,2);thermal=r[:,17]-.5*np.trace(Pi,axis1=1,axis2=2)
 rebuilt=p[:,:,None]*p[:,None,:]/m[:,None,None]+B@S@B.transpose(0,2,1)+R
 assert np.isfinite(r).all() and np.all(thermal>0) and np.all(np.linalg.eigvalsh(S)>0) and np.all(np.linalg.eigvalsh(R)>0)
 if reference:thermal=P0*area/(gamma-1)
 v=p/m[:,None];pressure=(gamma-1)*thermal/area
 result=dict(raw=r,mass=m,density=m/area,centroid=q/m[:,None],velocity=v,radial_velocity=np.sum(v*er,axis=1),azimuthal_velocity=np.sum(v*et,axis=1),intrinsic_ell=C[:,0,1]-C[:,1,0],origin_L=D[:,0,1]-D[:,1,0],R=R,thermal=thermal,specific_thermal=thermal/m,pressure=pressure,Psi=saved['Psi'].astype(ld))
 diagnostics=dict(minimum_thermal=float(thermal.min()),minimum_R_eigenvalue=float(np.linalg.eigvalsh(R).min()),maximum_Pi_factor_defect=float(abs(Pi-rebuilt).max()),minimum_pressure=float(pressure.min()),maximum_pressure=float(pressure.max()))
 return result,diagnostics
fields,diagnostics=decode(actual);target,reference_diagnostics=decode(ref,True)
phase_defined=(np.sum(fields['velocity']**2,axis=1)>0)&(np.sum(target['velocity']**2,axis=1)>0)
v=fields['velocity'];vr=target['velocity'];phase_error=np.degrees(np.arctan2(vr[:,0]*v[:,1]-vr[:,1]*v[:,0],np.sum(vr*v,axis=1)))
lo=ref['minimum_domain_radius'];hi=ref['maximum_domain_radius'];wall=ref['wall_touching'];padding=32*np.finfo(float).eps*np.maximum(1,hi)
contact_inner=(lo-padding<=.5)&(hi+padding>=.5);contact_outer=(lo-padding<=2)&(hi+padding>=2);initial_radius=np.linalg.norm(initial['state'][:-4].reshape(n,20)[:,:2],axis=1)
regions=dict(all=np.ones(n,bool),inner=initial_radius<.4,cavity=(hi<.5)&~contact_inner,disk_core=(lo>.5)&(hi<2)&~contact_inner&~contact_outer,inner_contact=contact_inner,outer_contact=contact_outer,dynamic_buffer=(lo>2)&~wall&~contact_outer,wall=wall)
comparison={};weights=actual['mass']
for name,mask in regions.items():
 row=dict(cells=int(mask.sum()),fields={})
 for key in fields:
  chosen=mask&defined if key in ('radial_velocity','azimuthal_velocity') else mask
  count=int(chosen.sum())
  if not count:continue
  delta=np.asarray(fields[key][chosen]-target[key][chosen],dtype=ld).reshape(count,-1);reference=np.asarray(target[key][chosen],dtype=ld).reshape(count,-1);w=weights[chosen].astype(ld)
  rms=np.sqrt(np.sum(w*np.sum(delta*delta,axis=1))/w.sum());reference_rms=np.sqrt(np.sum(w*np.sum(reference*reference,axis=1))/w.sum())
  row['fields'][key]=dict(maximum_absolute_error=float(abs(delta).max()),mass_weighted_RMS_error=float(rms),reference_mass_weighted_RMS=float(reference_rms),relative_RMS=float(rms/reference_rms) if reference_rms>0 else None)
 if mask.any():
  pe=(fields['pressure'][mask]-P0)/P0;me=(fields['mass'][mask]-target['mass'][mask])/target['mass'][mask]
  row['relative_pressure_error']=dict(RMS=float(np.sqrt(np.average(pe*pe,weights=weights[mask]))),maximum=float(abs(pe).max()))
  row['relative_mass_reference_error']=dict(RMS=float(np.sqrt(np.average(me*me,weights=weights[mask]))),maximum=float(abs(me).max()))
 phase_mask=mask&phase_defined
 if phase_mask.any():row['velocity_direction_error_degrees']=dict(cells=int(phase_mask.sum()),RMS=float(np.sqrt(np.average(phase_error[phase_mask]**2,weights=weights[phase_mask]))),maximum=float(abs(phase_error[phase_mask]).max()))
 comparison[name]=row
initial_raw=np.column_stack([initial['mass'],initial['state'][:-4].reshape(n,20)[:,3:]]).astype(ld);raw=fields['raw'].astype(ld);rawref=target['raw'].astype(ld)
global_values=dict(raw_totals_minus_reference=np.asarray((raw-rawref).sum(0),float).tolist(),origin_L_change=float(((raw[:,10]-raw[:,11])-(initial_raw[:,10]-initial_raw[:,11])).sum()),hydro_energy_change=float((raw[:,17]-initial_raw[:,17]).sum()),Psi_change=float((actual['Psi'].astype(ld)-initial['Psi'].astype(ld)).sum()),combined_energy_change=float((raw[:,17]-initial_raw[:,17]).sum()+(actual['Psi'].astype(ld)-initial['Psi'].astype(ld)).sum()),stored_P_E_L_ledger=actual['state'][-4:].tolist())
initial_density_report=json.loads(Path('/gpfs/kjhan/LagForce/Codex_full_disk_initial_density_20261019/report.json').read_text())
reconstruction_signed=density_report['signed_q48_Psi_error'];initial_reconstruction_signed=initial_density_report['signed_q48_Psi_error']
global_values.update(Psi_reconstruction_q48_defect_signed=reconstruction_signed,Psi_reconstruction_q48_defect_absolute_sum=density_report['absolute_q48_Psi_error'],initial_Psi_reconstruction_q48_defect_signed=initial_reconstruction_signed,combined_energy_change_using_q48_Psi_vs_IC_target=global_values['combined_energy_change']+reconstruction_signed,combined_energy_change_q48_vs_initial_q48=global_values['combined_energy_change']+reconstruction_signed-initial_reconstruction_signed,measured_reconstruction_defect_absolute_sums=density_report['absolute_q48_Psi_error']+initial_density_report['absolute_q48_Psi_error'])
raw_change=raw-initial_raw
update=np.r_[raw_change[:,7:9].sum(0),raw_change[:,17].sum(),(raw_change[:,10]-raw_change[:,11]).sum()]
global_values['P_E_L_update_ledger_residual']=np.asarray(update-(actual['state'][-4:]-initial['state'][-4:]),float).tolist()
np.savez(args.output/'fields.npz',ids=actual['ids'],velocity_direction_error_degrees=phase_error,phase_defined=phase_defined,area=area,radius=radius,relative_pressure_error=(fields['pressure']-P0)/P0,relative_density_reference_error=(fields['density']-target['density'])/target['density'],**{k:mask for k,mask in regions.items()},**{'actual_'+k:v for k,v in fields.items()},**{'reference_'+k:v for k,v in target.items()})
report=dict(status='completed_physical_state_comparison',time=float(actual['time']),scope=__doc__,state_sha256=state_sha,diagnostics=diagnostics,reference_diagnostics=reference_diagnostics,regions=comparison,global_values=global_values,mask_definition=dict(inner='Original stable registry radius<0.4, retained from previous inner diagnostics',contact='Working radial-bound envelope of each actual curved cell; density radii0.5/2.0',disk_core='Cell entirely between density contacts',buffer='Cell entirely beyond r2 and not touching wall',wall='Actual circular boundary arcs'),provider_calls=0,evolution_steps=0,limits='One state against its actual-cell steady reference. Numerical positivity/conservation alone does not prove physical accuracy; no time-order or long-term stability claim.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],time=report['time'],pressure={k:v.get('relative_pressure_error') for k,v in comparison.items()},global_values=global_values),indent=2))
