from pathlib import Path
import struct,json,numpy as np,subprocess,re,csv
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=Path(__file__).resolve().parent
gamma=float(np.float32(1.4));exact=pkg/'general2d_1d_references_v1/exact_sod'
def total(a):return float(np.sum(a,dtype=np.longdouble))
def read(path):
 with path.open('rb') as f:h=f.read(64);f.seek(-52,2);footer=f.read(52)
 magic,version,n,stride,prefix,features=struct.unpack('<iIQIIQ',h[:32]);assert magic==-2147483607 and (version,stride,prefix)==(2,696,576)
 assert path.stat().st_size==64+n*stride+52
 raw=np.memmap(path,dtype='u1',mode='r',offset=64,shape=(n,stride))
 def c(o,t='f8'):return np.ndarray((n,),dtype=t,buffer=raw,offset=o,strides=(stride,)).copy()
 ids=c(8,'u8')&np.uint64(2**56-1);order=np.argsort(ids);assert len(np.unique(ids))==n
 a={k:c(o)[order] for k,o in [('x',24),('y',32),('m',48),('vx',56),('vy',64),('W',576),('pi',584),('U',672),('V',688)]};a['id']=ids[order]
 a['time']=struct.unpack('<d',footer[:8])[0];a['rho']=a['m']/a['V'];a['P']=(gamma-1)*a['U']/a['V'];a['E']=total(a['U']+.5*a['m']*(a['vx']**2+a['vy']**2));a['mom']=np.array([total(a['m']*a[k]) for k in ['vx','vy']])
 return a
rows=[]
for d in list(root.glob('*_N*')):
 if not (d/'result.json').exists():continue
 result=json.loads((d/'result.json').read_text())
 if result['exit']:rows.append(result);continue
 problem=result.get('problem','shear');n=result.get('Nside',result.get('N'));angle=result.get('angle','axis')
 initial=Path(json.loads((d/'process.json').read_text())['initial']) if (d/'process.json').exists() else (pkg/'general2d_trial_v1'/f'shear_N{n}_axis'/'initial.bin' if n!=128 else d/'initial.bin')
 a=read(initial);outputs=[f for f in d.glob('khout.*.dat') if not f.is_symlink()];assert len(outputs)==1
 b=read(outputs[0]);assert np.array_equal(a['id'],b['id'])
 row={**result,'problem':problem,'Nside':n,'angle':angle,'time':b['time'],'cells':len(b['id']),'max_mass_change':float(abs(a['m']-b['m']).max()),'relative_total_mass_change':(total(b['m'])-total(a['m']))/total(a['m']),'relative_energy_change':(b['E']-a['E'])/a['E'],'momentum_change':list(b['mom']-a['mom']),'volume_sum':total(b['V']),'max_W':float(abs(b['W']).max()),'max_pi':float(abs(b['pi']).max()),'min_U':float(b['U'].min()),'min_V':float(b['V'].min()),'min_entropy_ratio':float((b['U']/a['U']*(b['V']/a['V'])**(gamma-1)).min())}
 log=(d/'log').read_text();steps=re.findall(r'UNIFIED_KDK\] step=\d+ t=(\S+) dt=(\S+).*? total_s=(\S+)',log);row['accepted_steps']=len(steps);row['step_wall_seconds']=sum(float(r[2]) for r in steps)
 if problem in ['contact','moving_contact']:
  expected=(0,0) if problem=='contact' else (1,.5)
  row['max_periodic_position_error']=float(np.maximum(abs((b['x']-a['x']-expected[0]*b['time']+1)%2-1),abs((b['y']-a['y']-expected[1]*b['time']+1)%2-1)).max())
  row.update({'density_L1_against_initial_material_cell_averages':total(b['V']*abs(b['rho']-a['rho']))/4,'max_pressure_error':float(abs(b['P']-1).max()),'max_velocity_error':float(np.hypot(b['vx']-expected[0],b['vy']-expected[1]).max())})
 elif problem=='shear':
  ev=np.sin(np.pi*b['y']);row.update({'rho_L1':total(b['V']*abs(b['rho']-1))/4,'pressure_L1':total(b['V']*abs(b['P']-1))/4,'velocity_L1':total(b['V']*np.hypot(b['vx']-ev,b['vy']))/4,'max_abs_transverse_velocity':float(abs(b['vy']).max()),'relative_internal_energy_increase':(total(b['U'])-total(a['U']))/total(a['U'])})
 elif problem=='sod':
  phase=(b['x']/2+(b['y']/2 if angle=='diagonal' else 0))%1
  mask=(phase>=.5)&(phase<1)
  xn=.5+(phase[mask]-.75)*2/(np.sqrt(2) if angle=='diagonal' else 1)
  arg=''.join(f'{x:.17g}\n' for x in xn)
  r=subprocess.run([str(exact),str(b['time']),str(gamma)],input=arg,text=True,capture_output=True,check=True)
  reference=np.array([[float(x) for x in line.split(',')] for line in r.stdout.splitlines()])
  vn=(b['vx'][mask]+b['vy'][mask])/np.sqrt(2) if angle=='diagonal' else b['vx'][mask]
  vt=(b['vy'][mask]-b['vx'][mask])/np.sqrt(2) if angle=='diagonal' else b['vy'][mask]
  norm=total(b['V'][mask]);row['exact_point_sample_L1_volume_mean']={key:total(b['V'][mask]*abs(values-reference[:,k]))/norm for key,values,k in [('rho',b['rho'][mask],1),('normal_velocity',vn,2),('pressure',b['P'][mask],3)]}
  row['transverse_velocity_volume_L1']=total(b['V'][mask]*abs(vt))/norm;row['max_abs_transverse_velocity']=float(abs(vt).max())
  if angle=='axis':
   # Compare material columns with the matching one-dimensional solution.
   ids=b['id'];cols=(ids%n).astype(int);V=np.bincount(cols,weights=b['V'],minlength=n);M=np.bincount(cols,weights=b['m'],minlength=n);X=np.bincount(cols,weights=b['V']*b['x'],minlength=n)/V
   rho=M/V;vx=np.bincount(cols,weights=b['m']*b['vx'],minlength=n)/M;P=np.bincount(cols,weights=b['V']*b['P'],minlength=n)/V
   keep=(X>=1)&(X<=2);X=X[keep]-1;rho=rho[keep];vx=vx[keep];P=P[keep];V=V[keep]/2
   with (pkg/'general2d_1d_references_v1'/f'N{n//2}'/'final.csv').open() as f:dat=list(csv.DictReader(f))
   xx=np.array([float(r['x']) for r in dat]);faces=np.r_[0,.5*(xx[:-1]+xx[1:]),1];ix=np.clip(np.searchsorted(faces,X,side='right')-1,0,len(xx)-1)
   row['same_resolution_1d_point_comparison']={k:float(np.sum(V*abs(a2-np.array([float(r[k]) for r in dat])[ix]))) for k,a2 in [('rho',rho),('v',vx),('p',P)]}
   (d/'column_profile.csv').write_text('x,rho,v,p,width\n'+''.join(','.join(f'{z:.17g}' for z in r)+'\n' for r in zip(X,rho,vx,P,V)))
 rows.append(row)
(root/'assessment.json').write_text(json.dumps({'scope':'actual full periodic native 2D evolution; point-sampled analytic norms are distinguished from exact cell-average quadrature','gamma':gamma,'cases':rows},indent=2)+'\n')
for r in rows:
 print(r.get('name',f"shear_N{r.get('N')}"),'exit',r['exit'],'E',r.get('relative_energy_change'),'steps',r.get('accepted_steps'))
 if r['exit']==0:
  print({k:v for k,v in r.items() if k in ['density_L1_against_initial_material_cell_averages','max_pressure_error','max_velocity_error','rho_L1','pressure_L1','velocity_L1','exact_point_sample_L1_volume_mean','transverse_velocity_volume_L1','same_resolution_1d_point_comparison','momentum_change','min_entropy_ratio']})
