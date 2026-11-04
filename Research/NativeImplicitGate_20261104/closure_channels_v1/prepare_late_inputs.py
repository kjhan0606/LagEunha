from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,os,time,shutil,numpy as np,struct
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=pkg/'closure_channels_v1'
exe=pkg/'native_coupled_build_v5/unified_native_coupled.exe';template=pkg/'native_rhs_comparison/t0_unified'
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
def columns(p):
 raw=np.memmap(p,dtype='u1',mode='r',offset=64,shape=(90352,696))
 def col(o,t='f8'):return np.ndarray((90352,),dtype=t,buffer=raw,offset=o,strides=(696,)).copy()
 ids=col(8,'u8')&np.uint64(2**56-1);order=np.argsort(ids)
 return ids[order],np.column_stack([col(o)[order] for o in [48,24,32,56,64]])
ic=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_w0_initial_20261101_v1/initial696.bin');ids,icv=columns(ic);R0=np.hypot(icv[:,1]-3,icv[:,2]-3)
def run(target):
 parent=pkg/'native_finite_heat_t01_v1'/f'base_to{target:g}';source=next(f for f in parent.glob('khout.*.dat') if not f.is_symlink());step=int(source.name.split('.')[1])
 d=root/f'rhs_t{target:g}';d.mkdir(exist_ok=False);(d/'params.dat').write_text((template/'params.dat').read_text());(d/f'khout.{step:06d}.dat').symlink_to(source)
 env=json.loads((template/'environment.json').read_text());env.update({'HYDRO_TSTOP':repr(target+.001),'I_MPI_PIN':'0','I_MPI_HYDRA_BOOTSTRAP':'fork','LD_LIBRARY_PATH':'/opt/ohpc/pub/cuda/13.0.2/lib64:'+os.environ.get('LD_LIBRARY_PATH','')})
 (d/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
 with (d/'log').open('x') as out:r=subprocess.run([str(exe),'params.dat',str(step)],cwd=d,env={**base,**env},stdout=out,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError((target,r.returncode,(d/'log').read_text()[-3000:]))
 cell=np.fromfile(d/'weight_native_cells.bin',dtype=[('id','i8'),('v','f8',(10,))]);aux=np.fromfile(d/'unified_rhs_aux.bin',dtype=[('id','i8'),('v','f8',(8,))]);assert np.array_equal(cell['id'],aux['id'])
 sourceids,values=columns(source);at=np.searchsorted(sourceids,cell['id']);assert np.array_equal(sourceids[at],cell['id'])
 at0=np.searchsorted(ids,cell['id']);assert np.array_equal(ids[at0],cell['id'])
 output=np.zeros(len(cell),dtype=[('id','i8'),('v','f8',(12,))]);output['id']=cell['id'];v=output['v'];v[:,0]=values[at,0];v[:,1]=cell['v'][:,7];v[:,2]=cell['v'][:,9];v[:,3:7]=values[at,1:5];v[:,7:11]=aux['v'][:,2:6];v[:,11]=R0[at0]
 meta=json.loads((d/'weight_native_meta.json').read_text());dst=root/f't{target:g}.bin'
 with dst.open('xb') as f:f.write(struct.pack('<QQQd',0x554e494632444131,len(cell),meta['faces'],meta['Gamma']));output.tofile(f);f.write((d/'weight_native_faces.bin').read_bytes())
 print(target,'RHS only complete',flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[.03,.1]))
