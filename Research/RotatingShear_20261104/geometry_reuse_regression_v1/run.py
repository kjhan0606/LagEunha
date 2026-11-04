from pathlib import Path
import os,subprocess,json,time,shutil,re,hashlib,struct,numpy as np
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=Path(__file__).resolve().parent;exe=pkg/'geometry_reuse_build_v1/unified_native_general.exe'
cases=[('moving_contact_N64_diagonal',.5),('shear_N64_axis',.25),('sod_N64_axis',.2)]
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
(root/'registration.json').write_text(json.dumps({'purpose':'same-binary A/B geometry reuse regression, sequential single-process profiling; external GPU load uncontrolled','cases':cases,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'outputs':'one final snapshot per run; original ICs linked','expected_new_bytes':20000000,'free_bytes':shutil.disk_usage(root).free},indent=2)+'\n')
def read(p):
 with p.open('rb') as f:h=f.read(64);f.seek(-52,2);t=struct.unpack('<d',f.read(8))[0]
 n=struct.unpack_from('<Q',h,8)[0];raw=np.memmap(p,dtype='u1',mode='r',offset=64,shape=(n,696))
 def c(offset,t='f8'):return np.ndarray((n,),dtype=t,buffer=raw,offset=offset,strides=(696,)).copy()
 ids=c(8,'u8')&np.uint64(2**56-1);order=np.argsort(ids)
 return {k:c(o)[order] for k,o in [('x',24),('y',32),('m',48),('vx',56),('vy',64),('W',576),('pi',584),('U',672),('V',688)]},t
rows=[];comparisons=[]
for name,target in cases:
 source=pkg/'solid_angle_periodic2d_v1'/name;initial=Path(json.loads((source/'process.json').read_text())['initial'])
 for reuse in [0,1]:
  d=root/f'{name}_reuse{reuse}';d.mkdir(exist_ok=False);shutil.copy2(source/'params.dat',d/'params.dat');(d/'khout.000001.dat').symlink_to(initial)
  env=json.loads((source/'environment.json').read_text());env.pop('GFS_UNIFIED_PLANAR_AUDIT_N',None);env['GFS_UNIFIED_REUSE_GEOMETRY']=str(reuse);(d/'environment.json').write_text(json.dumps(env,indent=2)+'\n');start=time.monotonic()
  with (d/'log').open('x') as f:
   p=subprocess.Popen([str(exe),'params.dat','1'],cwd=d,env={**base,**env},stdout=f,stderr=subprocess.STDOUT);(d/'process.json').write_text(json.dumps({'pid':p.pid,'initial':str(initial)},indent=2)+'\n');rc=p.wait()
  tlog=(d/'log').read_text();times=re.findall(r'geometry_s=(\S+) pressure_s=(\S+) viscosity_s=(\S+) drift_s=(\S+) total_s=(\S+)',tlog)
  metrics=np.array([[float(x) for x in t] for t in times]);reuse_count=tlog.count('[UNIFIED_GEOMETRY_REUSE] reused=1')
  row={'case':name,'reuse':reuse,'exit':rc,'wall_seconds':time.monotonic()-start,'steps':len(times),'reused_steps':reuse_count,'phase_seconds':dict(zip(['geometry','pressure','viscosity','drift','total'],metrics.sum(axis=0)))}
  assert rc==0;outputs=[x for x in d.glob('khout.*.dat') if not x.is_symlink()];assert len(outputs)==1
  state,t=read(outputs[0]);assert t==target;row['time']=t;(d/'result.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(row),flush=True)
  if reuse==0:reference=state
  else:
   differences={}
   for k,v in state.items():
    delta=v-reference[k]
    if k in ['x','y']:delta-=2*np.round(delta/2)
    scale=max(float(abs(v).max()),float(abs(reference[k]).max()),1e-30) if k in ['m','U','V'] else max(1,float(abs(v).max()),float(abs(reference[k]).max()))
    differences[k]={'max_abs':float(abs(delta).max()),'max_scaled':float(abs(delta).max())/scale}
   pass_=all(x['max_scaled']<2e-10 for x in differences.values()) and reuse_count==len(times)-1
   comparisons.append({'case':name,'fields':differences,'pass':pass_});print(json.dumps(comparisons[-1]),flush=True)
(root/'results.json').write_text(json.dumps(rows,indent=2)+'\n');(root/'comparisons.json').write_text(json.dumps(comparisons,indent=2)+'\n')
assert all(c['pass'] for c in comparisons)
