"""Restart the unchanged material problem in an isolated MPI qualification directory."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,struct,subprocess,time
p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('label');p.add_argument('--ranks',type=int,required=True);p.add_argument('--window',type=float,default=.0005);p.add_argument('--dt',type=float,default=.000025);p.add_argument('--volume',default='endpoint');p.add_argument('--binary',type=Path);p.add_argument('--step',type=int);p.add_argument('--timeout',type=int,default=600);a=p.parse_args()
O=Path(__file__).resolve().parent;R=O.parents[3];S=a.source.resolve();T=O/a.label;T.mkdir()
prefix='cylmat' if list(S.glob('cylmat.*.dat')) else 'rtout' if list(S.glob('rtout.*.dat')) else 'khout'
f=S/(prefix+'.%06d.dat'%a.step) if a.step else sorted(S.glob(prefix+'.*.dat'))[-1]
step=int(f.name.split('.')[1]);binary=(a.binary or O/'build_a10/native.exe').resolve()
for x in S.glob(f.name+'*'):shutil.copy2(x,T/x.name)
shutil.copy2(S/'params.dat',T/'params.dat')
with f.open('rb') as q:q.seek(-52,2);start=struct.unpack('d',q.read(8))[0]
env=os.environ.copy();env.update(json.loads((S/'environment.json').read_text()));env.update(json.loads((R/'Research/MaterialVolume_20261105/material_finite_slack_20261107_v1/native_runtime_paths.json').read_text()))
env.update(HYDRO_TSTOP=repr(start+a.window),GFS_MATERIAL_DT=repr(a.dt),GFS_MATERIAL_VOLUME_COORDINATE=a.volume,OMP_NUM_THREADS='1',OMP_PROC_BIND='false',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',I_MPI_FABRICS='shm',I_MPI_PIN='0',CUDA_VISIBLE_DEVICES=os.environ.get('GFS_TEST_DEVICES','2,3'))
for k in ['GFS_MATERIAL_TRACE_FAILURE_DUMP','GFS_MATERIAL_AUDIT_FIRST_RHS','GFS_CUDA_PROJECT_AUDIT','GFS_MATERIAL_FIXED_IC']:env.pop(k,None)
cmd=['mpirun','-n',str(a.ranks),str(binary),'params.dat',str(step)]
(T/'environment.json').write_text(json.dumps(env,indent=2));(T/'command.json').write_text(json.dumps(cmd));begin=time.monotonic()
with (T/'log').open('w') as out:
 try:rc=subprocess.run(cmd,cwd=T,env=env,stdout=out,stderr=out,timeout=a.timeout).returncode
 except subprocess.TimeoutExpired:rc=124
status={'returncode':rc,'host':os.uname().nodename,'ranks':a.ranks,'source':str(f),'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'initial_t':start,'target':start+a.window,'dt_cap':a.dt,'volume':a.volume,'wall_seconds':time.monotonic()-begin}
(T/'status.json').write_text(json.dumps(status,indent=2));print(a.label,json.dumps(status),flush=True)
if rc:print((T/'log').read_text()[-2000:])
raise SystemExit(rc)
