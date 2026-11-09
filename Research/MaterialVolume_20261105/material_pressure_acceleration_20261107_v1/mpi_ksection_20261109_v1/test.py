from pathlib import Path
import os,json,shutil,subprocess,sys,struct,time
O=Path(__file__).resolve().parent;D=O.parent;R=O.parents[3];B=O/'build'
case=sys.argv[1];rank=int(sys.argv[2]);label=sys.argv[3] if len(sys.argv)>3 else 'a'
P=D/'common_noh_volume_20261109_v1/regression_sequential'/case
if case=='noh':P=D/'common_noh_volume_20261109_v1/short/n32_endpoint'
T=O/(case+'_r'+str(rank)+'_'+label);T.mkdir()
step=int(json.loads((P/'command.json').read_text())[-1]) if (P/'command.json').exists() else 1;prefix='rtout' if case=='rt' else 'cylmat' if case=='cylinder' else 'khout'
name=prefix+'.%06d.dat'%step
for f in P.glob(name+'*'):shutil.copy2(f,T/f.name)
shutil.copy2(P/'params.dat',T/'params.dat')
with (T/name).open('rb') as f:f.seek(-52,2);t=struct.unpack('d',f.read(8))[0]
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
env.update(json.loads((R/'Research/MaterialVolume_20261105/material_finite_slack_20261107_v1/native_runtime_paths.json').read_text()))
env.update(json.loads((P/'environment.json').read_text()))
for k in ['OMP_PLACES','GFS_MATERIAL_TRACE_FAILURE_DUMP','GFS_MATERIAL_AUDIT_FIRST_RHS','GFS_CUDA_PROJECT_AUDIT']:env.pop(k,None)
env.update(HYDRO_TSTOP=str(t+0.000025),GFS_MATERIAL_DT='0.000025',GFS_MATERIAL_VOLUME_COORDINATE='endpoint',OMP_NUM_THREADS='1',OMP_PROC_BIND='false',I_MPI_FABRICS='shm',I_MPI_PIN='0')
env['CUDA_VISIBLE_DEVICES']=os.environ.get('GFS_TEST_DEVICES','0')
cmd=['mpirun','-n',str(rank),os.environ.get('GFS_TEST_BINARY',str(B/'native.exe')),'params.dat',str(step)]
(T/'environment.json').write_text(json.dumps(env,indent=2));(T/'command.json').write_text(json.dumps(cmd))
a=time.monotonic()
with (T/'log').open('w') as f:
 try:q=subprocess.run(cmd,cwd=T,env=env,stdout=f,stderr=f,timeout=180);rc=q.returncode
 except subprocess.TimeoutExpired:rc=124
(T/'status.json').write_text(json.dumps({'returncode':rc,'host':os.uname().nodename,'binary':cmd[3],'devices':env['CUDA_VISIBLE_DEVICES'],'wall_seconds':time.monotonic()-a,'initial_t':t,'target':t+0.000025}))
print(T,rc,flush=True)
print((T/'log').read_text()[-1000:])

raise SystemExit(rc)
