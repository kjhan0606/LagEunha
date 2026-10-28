from pathlib import Path
import json,struct,re,hashlib,numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_early_growth_budget_20261017')
layout=json.loads(Path('/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json').read_text())
fields=['x','y','mass','vx','vy']; dtp=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*5,'offsets':[layout['u4if']]+[layout[k] for k in fields],'itemsize':layout['size']})
def read(p):
 raw=p.read_bytes(); magic,n,size,stress=struct.unpack_from('<iiii',raw); assert magic==-1 and size==layout['size'] and len(raw)==16+n*size+52
 a=np.frombuffer(raw,dtype=dtp,count=n,offset=16); ids=a['id']&layout['id_mask']; order=np.argsort(ids); return a[order],ids[order],struct.unpack_from('<d',raw,16+n*size)[0]
def ep(a):
 x,y=a['x']-3,a['y']-3; r=np.hypot(x,y); u=r/.5; m=np.ones_like(r); lo=u<=.5; mid=(u>.5)&(u<1)
 q=u[lo]; m[lo]=q**3*(32/3-192/5*q*q+32*q**3); q=u[mid]; m[mid]=-1/15+q**3*(64/3-48*q+192/5*q*q-32/3*q**3)
 om=np.sqrt(m/r**3); return .5*a['mass']*((a['vx']+y*om)**2+(a['vy']-x*om)**2),r

results={}
for entry in json.loads((root/'manifest.json').read_text()):
 case=entry['case'];d=root/case;assert (d/'exit.code').read_text().strip()=='0'
 step=entry['step'];start=d/f'khout.{step:06d}.dat';finish=d/f'khout.{step+1:06d}.dat'
 before,ib,t0=read(start);after,ia,t1=read(finish)
 assert np.array_equal(ib,ia) and np.array_equal(before['mass'],after['mass'])
 log=(d/'log').read_text();lines=[s for s in log.splitlines() if s.startswith('[ACCEPTED_STATE]')];assert len(lines)==1
 state=dict(re.findall(r'(\w+)=([^\s]+)',lines[0]));assert state['retries']=='0' and state['step']=='1'
 dt=float(state['dt']);assert abs((t1-t0)-dt)<2e-15
 e0,r0=ep(before);e1,r1=ep(after)
 def bins(r):return np.where(r>=2.5,3,np.where(r<.5,0,np.where(r<=2,1,2)))
 b0,b1=bins(r0),bins(r1);cross=int(np.count_nonzero(b0!=b1))
 hg={i:[] for i in range(4)}
 for line in log.splitlines():
  if line.startswith('[HG_PERTURBATION]'):
   v=dict(re.findall(r'(\w+)=([^\s]+)',line));hg[int(v['bin'])].append(float(v['dEpert']))
 assert all(len(v)==2 for v in hg.values())
 regions={}
 for i,name in enumerate(['inner','disk','outer_active','frozen']):
  mask=b0==i;membership_equal=bool(np.array_equal(b0==i,b1==i));delta=float((e1-e0)[mask].sum());h=sum(hg[i]);energy=float(e0[mask].sum())
  regions[name]={'n':int(mask.sum()),'Epert_start':energy,'full_step_delta':delta,'HG_delta':h,'HG_membership_matches_material_bin':membership_equal,'RK_delta_by_subtraction':delta-h if membership_equal else None,'net_rate':delta/dt,'HG_rate':h/dt,'RK_rate_by_subtraction':(delta-h)/dt if membership_equal else None,'instantaneous_amplitude_growth_estimate':delta/(2*energy*dt) if energy>0 else None}
 results[case]={'t_start':t0,'t_end':t1,'dt':dt,'mass_unchanged':True,'radial_bin_crossings':cross,'regions':regions,'binary':(d/'binary.sha256').read_text().strip(),'snapshot_sha256':hashlib.sha256(start.read_bytes()).hexdigest()}
comparison={}
for label in ['t2','t4','t6']:
 a=results[label+'_full']['regions']['inner'];b=results[label+'_half']['regions']['inner']
 comparison[label]={'full_dt_net_rate':a['net_rate'],'half_dt_net_rate':b['net_rate'],'relative_rate_change':(b['net_rate']-a['net_rate'])/abs(a['net_rate']),'half_dt_HG_rate':b['HG_rate'],'half_dt_RK_rate':b['RK_rate_by_subtraction'],'half_dt_amplitude_growth_estimate':b['instantaneous_amplitude_growth_estimate']}
comparison['t6']['quarter_dt_net_rate']=results['t6_quarter']['regions']['inner']['net_rate']
comparison['t6']['quarter_vs_half_relative_change']=(comparison['t6']['quarter_dt_net_rate']-comparison['t6']['half_dt_net_rate'])/abs(comparison['t6']['half_dt_net_rate'])
out={'scope':'One accepted full HG/RK4/HG step at each pre-existing evolved checkpoint, with half-step-size check. Endpoint bins follow start-time particle IDs; HG bins use instantaneous radius. RK subtraction is supplied only when membership is identical at both endpoints. One particle crosses r=2 in each t2 case, so those disk/outer RK subtractions are omitted. RK contribution is obtained by subtracting measured HG increments from endpoint increment, not an independent stage budget. Local amplitude-rate estimates are not long-time fitted exponents.','cases':results,'inner_comparison':comparison}
(root/'growth_budget.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(comparison,indent=2))
