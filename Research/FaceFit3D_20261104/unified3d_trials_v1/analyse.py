from pathlib import Path
import json,struct,subprocess,numpy as np,re
root=Path(__file__).resolve().parent;pkg=root.parent
exact=pkg/'general2d_1d_references_v1/exact_sod'
def total(x):return float(np.sum(x,dtype=np.longdouble))
def read(path):
 with path.open('rb') as f:h=f.read(128)
 assert h[:8]==b'LAG3DV1\0'
 version,flags,n,t,gamma=struct.unpack_from('<iiqdd',h,8);assert version==1 and flags==2
 assert path.stat().st_size==128+88*n
 a=np.fromfile(path,dtype='<f8',count=10*n,offset=128).reshape(10,n)
 ids=np.fromfile(path,dtype='<i8',count=n,offset=128+80*n);order=np.argsort(ids);assert len(np.unique(ids))==n
 a=a[:,order];m=a[6];U=m*a[7];V=a[9];v=a[3:6].T
 return {'id':ids[order],'time':t,'gamma':gamma,'x':a[:3].T,'v':v,'m':m,'U':U,'V':V,'rho':m/V,'P':(gamma-1)*U/V,'E':total(U+.5*m*np.sum(v*v,axis=1)),'momentum':np.array([total(m*v[:,a]) for a in range(3)])}
rows=[]
for d in sorted(root.glob('*_N*')):
 r=json.loads((d/'result.json').read_text())
 if r['exit']:
  log=(d/'log').read_text();r['failure_lines']=[x for x in log.splitlines() if 'failed' in x or 'invalid' in x];rows.append(r);continue
 a=read(d/'snap_000000.l3d');b=read(d/'snap_000001.l3d');assert np.array_equal(a['id'],b['id']);assert np.array_equal(a['m'],b['m'])
 r.update({'time':b['time'],'reached_target':b['time']==r['target'],'cells':len(b['id']),'max_mass_change':float(abs(a['m']-b['m']).max()),'min_U':float(b['U'].min()),'min_V':float(b['V'].min()),'min_entropy_ratio':float((b['U']/a['U']*(b['V']/a['V'])**(b['gamma']-1)).min()),'volume_sum':total(b['V']),'relative_energy_change':(b['E']-a['E'])/a['E'],'momentum_change':list(b['momentum']-a['momentum'])})
 if r['problem'] in ['contact','moving_contact']:
  velocity=np.array([1,.5,.25]) if r['problem']=='moving_contact' else np.zeros(3)
  r.update({'density_L1_material_cell_average':total(b['V']*abs(b['rho']-a['rho']))/8,'max_pressure_error':float(abs(b['P']-1).max()),'max_velocity_error':float(np.linalg.norm(b['v']-velocity,axis=1).max()),'max_periodic_position_error':float(abs((b['x']-a['x']-b['time']*velocity+1)%2-1).max())})
 elif r['problem']=='shear':
  expected=np.zeros_like(b['v']);expected[:,0]=np.sin(np.pi*(b['x'][:,1]+b['x'][:,2]))
  r.update({'rho_L1':total(b['V']*abs(b['rho']-1))/8,'pressure_L1':total(b['V']*abs(b['P']-1))/8,'velocity_L1':total(b['V']*np.linalg.norm(b['v']-expected,axis=1))/8})
 elif r['problem']=='sod':
  normal=np.ones(3)/np.sqrt(3) if r['angle']=='diagonal' else np.array([1,0,0]);divisor=np.sqrt(3) if r['angle']=='diagonal' else 1
  phase=((np.sum(b['x'],axis=1) if r['angle']=='diagonal' else b['x'][:,0])/2)%1
  keep=(phase>=.5)&(phase<1);coord=.5+(phase[keep]-.75)*2/divisor
  sample=subprocess.run([str(exact),str(b['time']),str(b['gamma'])],input=''.join(f'{x:.17g}\n' for x in coord),text=True,capture_output=True,check=True)
  ref=np.array([[float(x) for x in line.split(',')] for line in sample.stdout.splitlines()])
  vn=b['v'][keep]@normal;vt=b['v'][keep]-vn[:,None]*normal[None,:];norm=total(b['V'][keep]);r['exact_point_sample_L1_volume_mean']={k:total(b['V'][keep]*abs(v-ref[:,c]))/norm for k,v,c in [('rho',b['rho'][keep],1),('normal_velocity',vn,2),('pressure',b['P'][keep],3)]}
  r['transverse_velocity_L1']=total(b['V'][keep]*np.linalg.norm(vt,axis=1))/norm;r['max_transverse_velocity']=float(np.linalg.norm(vt,axis=1).max())
 rows.append(r)
(root/'assessment.json').write_text(json.dumps({'scope':'actual original native 3D CPU/OpenMP evolution, periodic, W0, same common closure; point-sampled shock norms, not exact-cell-average norms','cases':rows},indent=2)+'\n')
for r in rows:print(json.dumps(r))
