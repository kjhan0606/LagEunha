"""Locate per-branch and central-response acceleration discrepancies."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_generator_steps_20261018'
src=b/'Codex_coupled_response_forces_t3p4_20261018'
fields=np.load(src/'force_fields.npz');initial={};mass=None
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
r0=np.hypot(ic['x']-3,ic['y']-3);masks={'inner':r0<.5,'shear_band':(r0>=.3)&(r0<.5),'central':r0<.2,'edge':(r0>=.5)&(r0<.6)}
def vectors(a):
 x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']));k=gravity_coefficient(np.linalg.norm(x,axis=1))
 vc=np.sqrt(k)[:,None]*np.column_stack((-x[:,1],x[:,0]))
 return x,v,vc,-k[:,None]*x
for sign in ['plus','minus']:
 a,t=read(src/sign/'khout.000743.dat');assert np.array_equal(a['id'],ic['id']);initial[sign]=vectors(a)
mass=ic['mass'];eps=.001
du0=((initial['plus'][1]-initial['plus'][2])-(initial['minus'][1]-initial['minus'][2]))/(2*eps)
def norm(z,mask):return float(np.sqrt(np.sum(mass[mask,None]*z[mask]**2)))
out={}
for dt in [1e-4,5e-5]:
 measured={};predicted={};ref_derivative={};branch={}
 for sign in ['plus','minus']:
  a,te=read(root/f'{sign}_{dt:g}'/'khout.000744.dat');assert np.array_equal(a['id'],ic['id'])
  assert abs(te-t-dt)<1e-12
  x,v,vc,grav=initial[sign];xx,vv,vcc,_=vectors(a)
  force=fields['hydro_'+sign]+fields['hg_'+sign]
  predicted[sign]=grav+force/mass[:,None]
  measured[sign]=(vv-v)/dt
  ref_derivative[sign]=(vcc-vc)/dt
  branch[sign]={}
  for name,mask in masks.items():
   branch[sign][name]={'acceleration_error_over_total':norm(measured[sign]-predicted[sign],mask)/norm(predicted[sign],mask),
    'acceleration_error_over_hydro_HG':norm(measured[sign]-predicted[sign],mask)/norm(force/mass[:,None],mask),
    'position_forward_derivative_error_over_velocity':norm((xx-x)/dt-v,mask)/norm(v,mask)}
 da=(measured['plus']-measured['minus'])/(2*eps)
 dp=(predicted['plus']-predicted['minus'])/(2*eps)
 dr=(ref_derivative['plus']-ref_derivative['minus'])/(2*eps)
 response_rate=da-dr
 contribution=mass*np.sum(du0*(da-dp),axis=1)
 regions={}
 for name,mask in masks.items():
  candidates=np.flatnonzero(mask);top=candidates[np.argsort(abs(contribution[candidates]))[-10:][::-1]]
  regions[name]={'central_acceleration_error_relative_L2':norm(da-dp,mask)/norm(dp,mask),
   'central_acceleration_error_power':float(contribution[mask].sum()),
   'measured_linear_departure_power':float(np.sum(mass[mask,None]*du0[mask]*response_rate[mask])),
   'forward_energy_quadratic_term':float(.5*dt*np.sum(mass[mask,None]*response_rate[mask]**2)),
   'largest_error_power_cells':[{'id':int(i),'initial_radius':float(r0[i]),'power':float(contribution[i])} for i in top]}
 out[str(dt)]={'branches':branch,'central_response':regions}
result={'comparisons':out,'scope':'Forward acceleration from saved one-step states versus terminal hydro plus reconstructed HG generator and gravity. Does not separate finite-time truncation, topology/reconstruction changes, and roundoff.'}
with (root/'acceleration_analysis.json').open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
