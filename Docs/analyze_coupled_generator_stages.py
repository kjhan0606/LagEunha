"""Locate response-generator mismatch in actual split/RK force calls."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_generator_stages_20261018'
probe=b/'Codex_coupled_response_forces_t3p4_20261018';control=b/'Codex_coupled_generator_omp1_20261018'
fields=np.load(probe/'force_fields.npz');dt=5e-5;eps=.001
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat');mass=ic['mass'];r0=np.hypot(ic['x']-3,ic['y']-3)
masks={'inner':r0<.5,'shear_band':(r0>=.3)&(r0<.5),'central':r0<.2,'edge':(r0>=.5)&(r0<.6)}
effects={};departures={};checks={}
def norm(a,mask):return float(np.sqrt(np.sum(mass[mask,None]*a[mask]**2)))
for sign in ['plus','minus']:
 d=root/f'{sign}_5e-05';assert (d/'exit.code').read_text().strip()=='0'
 a,t=read(d/'khout.000743.dat');z,te=read(d/'khout.000744.dat');c,tc=read(control/d.name/'khout.000744.dat')
 assert te==tc and all(np.array_equal(z[k],c[k]) for k in z.dtype.names)
 rows=[]
 for call in range(1,7):
  s=np.loadtxt(d/f'rhs.{call}');s=s[np.argsort(s[:,0])]
  assert s.shape==(len(a),25) and np.array_equal(s[:,0],a['id'] & 72057594037927935)
  rows.append(s)
 x=np.column_stack((a['x'],a['y']));v=np.column_stack((a['vx'],a['vy']));vf=np.column_stack((z['vx'],z['vy']))
 assert np.array_equal(rows[0][:,1:3],x) and np.array_equal(rows[0][:,3:5],v)
 assert np.array_equal(rows[1][:,1:3],x)
 assert np.array_equal(rows[5][:,1:3],np.column_stack((z['x'],z['y'])))
 h1=(rows[1][:,3:5]-rows[0][:,3:5])/dt
 h2=(vf-rows[5][:,3:5])/dt
 rk=(rows[5][:,3:5]-rows[1][:,3:5])/dt
 weighted=np.zeros_like(v)
 for w,s in zip([1/6,1/3,1/3,1/6],rows[1:5]):
  xx=s[:,1:3]-3;g=-gravity_coefficient(np.linalg.norm(xx,axis=1))[:,None]*xx
  weighted+=w*(s[:,8:10]/mass[:,None]+g)
 xx=x-3;g0=-gravity_coefficient(np.linalg.norm(xx,axis=1))[:,None]*xx
 hydro0=fields['hydro_'+sign]/mass[:,None];hg0=fields['hg_'+sign]/mass[:,None]
 departures[sign]=v-np.sqrt(gravity_coefficient(np.linalg.norm(xx,axis=1)))[:,None]*np.column_stack((-xx[:,1],xx[:,0]))
 effects[sign]={'first_HG_minus_half_generator':h1-.5*hg0,
 'last_HG_minus_half_generator':h2-.5*hg0,
 'RK_minus_initial_hydro_gravity':rk-hydro0-g0,
 'RK_minus_stage_quadrature':rk-weighted,
 'first_force_minus_terminal_probe':rows[0][:,8:10]/mass[:,None]-hydro0}
 checks[sign]={'dump_passive_all_reader_fields':True,
  'initial_force_relative_L2':{k:norm(effects[sign]['first_force_minus_terminal_probe'],m)/norm(hydro0,m) for k,m in masks.items()},
  'RK_quadrature_absolute_L2':{k:norm(rk-weighted,m) for k,m in masks.items()}}
du=(departures['plus']-departures['minus'])/(2*eps);regions={}
for name,mask in masks.items():
 regions[name]={}
 for term in effects['plus']:
  err=(effects['plus'][term]-effects['minus'][term])/(2*eps)
  power=mass*np.sum(du*err,axis=1);ids=np.flatnonzero(mask);top=ids[np.argsort(abs(power[ids]))[-5:][::-1]]
  regions[name][term]={'power':float(power[mask].sum()),'largest_cells':[{'id':int(i),'power':float(power[i])} for i in top]}
out={'time':t,'dt':dt,'checks':checks,'regions':regions,
 'scope':'Actual first/last HG impulses and RK stage RHS, compared with initial generator. Central response power projections are discrepancies, not physical source powers. Finite-dt effects and floating-point differences remain possible.'}
with (root/'stage_analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
