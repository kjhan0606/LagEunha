"""Bounded, user-approved KH continuation. Existing native solvers are unchanged."""
from pathlib import Path
import os,sys,json,shutil,hashlib,struct,subprocess,time,math,re,traceback,socket
import numpy as np
import h5py
O=Path(__file__).resolve().parent;D=O.parent;B=D.parent;R=B.parents[1]
nx=int(sys.argv[1]);kind=sys.argv[2];P=O/f'n{nx}';Q=P/kind;plan=json.loads((P/'plan.json').read_text())
cores=list(map(int,sys.argv[3].split(',')))
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_','OMP_','I_MPI_'))};base.update(json.loads((B/'material_finite_slack_20261107_v1/native_runtime_paths.json').read_text()));base.update(OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
extract=B/'material_kh_stop_20261106_v1/extract_ids'
def write(p,v):
 f=p.with_suffix(p.suffix+'.tmp');f.write_text(json.dumps(v,indent=2)+'\n');f.replace(p)
def state(s,**kw):write(Q/'status.json',dict(state=s,pid=os.getpid(),host=socket.gethostname(),updated_epoch=time.time(),**kw))
def when(p):
 with p.open('rb') as f:f.seek(-52,2);return struct.unpack('d',f.read(8))[0]
def load(p):
 a=np.fromfile(p).reshape(-1,13);return a[a[:,12].argsort()]
def quantities(a):
 a=a.astype(np.longdouble);m=a[:,6];k=np.sum(.5*m*np.sum(a[:,2:4]**2,axis=1));u=a[:,7].sum()
 return {k:float(v) for k,v in dict(mass=m.sum(),momentum_x=np.sum(m*a[:,2]),momentum_y=np.sum(m*a[:,3]),kinetic=k,internal=u,total_energy=k+u,area=a[:,8].sum()).items()}
def native_segment(folder,source,target,binary,params,e,initial):
 folder.mkdir();shutil.copy2(source,folder/'khout.000001.dat');(folder/'params.dat').write_text(params)
 e=dict(e,HYDRO_TSTOP=str(target),GFS_MATERIAL_TRACE_FAILURE_DUMP=str(folder/'failed_edges.bin'),OMP_NUM_THREADS=str(len(cores)),OMP_PLACES=','.join('{'+str(c)+'}' for c in cores),OMP_PROC_BIND='true',I_MPI_PIN='0')
 e['CUDA_VISIBLE_DEVICES']=os.environ.get('CUDA_VISIBLE_DEVICES','0')
 cmd=['taskset','-c',','.join(map(str,cores)),str(binary),'params.dat','1'];write(folder/'command.json',cmd);write(folder/'environment.json',e);write(folder/'lineage.json',dict(source=str(source),source_sha256=sha(source),start=when(source),target=target))
 tick=time.monotonic()
 with (folder/'log').open('x') as f:
  p=subprocess.Popen(cmd,cwd=folder,env=dict(base,**e),stdout=f,stderr=f);state('running',child_pid=p.pid,start=when(source),segment_target=target,target=10.,directory=str(folder));code=p.wait()
 lines=(folder/'log').read_text().splitlines();accepted=[x for x in lines if x.startswith('[MATERIAL_SSP_ACCEPT]')];vol=[float(re.search(r'max_relative=([^ ]+)',x)[1]) for x in lines if x.startswith('[MATERIAL_VOLUME_STEP]')]
 snaps=sorted(folder.glob('khout.*.dat'));snap=snaps[-1]
 r=dict(exit=code,saved_time=when(snap),target=target,steps=len(accepted),retries=sum(int(re.search(r'retries=(\d+)',x)[1]) for x in accepted),wall_seconds=time.monotonic()-tick,checkpoint=str(snap),checkpoint_sha256=sha(snap),max_volume_error=max(vol,default=None),passed=False)
 if code==0 and abs(r['saved_time']-target)<1e-12:
  subprocess.run([str(extract),str(snap),str(folder/'state_ids.bin')],check=True,stdout=subprocess.DEVNULL);a=load(folder/'state_ids.bin');qi=quantities(initial);q=quantities(a)
  e0=np.sum(initial[:,7].astype(np.longdouble)+.5*initial[:,6].astype(np.longdouble)*np.sum(initial[:,2:4].astype(np.longdouble)**2,axis=1));et=np.sum(a[:,7].astype(np.longdouble)+.5*a[:,6].astype(np.longdouble)*np.sum(a[:,2:4].astype(np.longdouble)**2,axis=1))
  r.update(fixed_mass_exact=bool(np.array_equal(a[:,12],initial[:,12]) and np.array_equal(a[:,6],initial[:,6])),finite=bool(np.isfinite(a).all()),min_U=float(a[:,7].min()),min_V=float(a[:,8].min()),quantities=q,energy_relative_change=float((et-e0)/e0),momentum_change=[q['momentum_x']-qi['momentum_x'],q['momentum_y']-qi['momentum_y']])
  r['passed']=bool(r['fixed_mass_exact'] and r['finite'] and r['min_U']>0 and r['min_V']>0 and r['max_volume_error'] is not None and r['max_volume_error']<=1e-11 and abs(r['energy_relative_change'])<1e-8)
 write(folder/'result.json',r);print(json.dumps(r),flush=True)
 if not r['passed']:raise RuntimeError(f'Native evolution failed; no retuning: {folder}')
 return snap,r
try:
 Q.mkdir(exist_ok=False);write(Q/'resources.json',dict(cores=cores,host=socket.gethostname(),mpi_ranks=1 if kind=='native' else 4,omp_threads=8 if kind=='native' else 2))
 if kind=='native':
  binary=Path(plan['native_binary']);assert sha(binary)==plan['native_binary_sha256'];source=P/'input/native.dat';initial=load(P/'input/t0_ids.bin');params=(P/'input/params.dat').read_text();e=json.loads((P/'input/environment.json').read_text())
  if nx==256:
   # Same KH portability case already qualified on Ada, before high-resolution H200 evolution.
   ref=D/'common_noh_volume_20261109_v1/regression_sequential/kh';first=min(ref.glob('khout.*.dat'));last=max(ref.glob('khout.*.dat'));subprocess.run([str(extract),str(first),str(Q/'gate_initial.bin')],check=True,stdout=subprocess.DEVNULL)
   _,gr=native_segment(Q/'h200_gate',first,when(last),binary,(ref/'params.dat').read_text(),json.loads((ref/'environment.json').read_text()),load(Q/'gate_initial.bin'))
   aa=load(Q/'h200_gate/state_ids.bin');bb=load(ref/'state_ids.bin');ok=bool(np.array_equal(aa[:,12],bb[:,12]) and np.array_equal(aa[:,6],bb[:,6]) and np.allclose(aa,bb,rtol=1e-8,atol=1e-10));write(Q/'h200_gate/parity.json',dict(passed=ok,max_absolute_columns=np.max(abs(aa-bb),axis=0).tolist()));assert ok
  results=[];start=plan['start'];index=0
  for major in plan['outputs']:
   steps=math.ceil((major-start)/.39);origin=start
   for j in range(1,steps+1):
    target=major if j==steps else origin+(major-origin)*j/steps;index+=1
    source,result=native_segment(Q/f'segment_{index:03d}',source,target,binary,params,e,initial);results.append(result);write(Q/'results.json',results)
   write(Q/f'endpoint_{major:g}.json',result);start=major
  state('completed',saved_time=10.,target=10.,segments=index,await_next_stage_approval=True)
 elif kind=='gizmo':
  assert sha(O/'bin/GIZMO')==plan['gizmo_binary_sha256'];shutil.copy2(P/'input/snap_001.hdf5',Q/'snap_001.hdf5');(Q/'output').mkdir()
  params=(B/'material_kh_128x256_20261106_v1/gizmo_t135/params.txt').read_text();params=params.replace('% GIZMO parameters for McNally 2D KH test (Hopkins Fig.16)','% Matched Lecoanet-style KH; not the McNally paper IC.')
  changes=dict(InitCondFile='snap_001',TimeBegin=str(plan['start']),TimeMax='10',OutputListOn='1',TimeOfFirstSnapshot=str(plan['start']),TimeBetSnapshot='20',TimeLimitCPU='165000',CpuTimeBetRestartFile='1800',MaxSizeTimestep=str(.001*128/nx),MaxMemSize='8000')
  for key,value in changes.items():params,n=re.subn(r'(?m)^'+key+r'\s+[^\n]+',key+' '+value,params);assert n==1,key
  params+='\nOutputListFilename outputs.txt\n';(Q/'params.txt').write_text(params);(Q/'outputs.txt').write_text('\n'.join(map(str,[plan['start'],*plan['outputs']]))+'\n')
  env=dict(base,OMP_NUM_THREADS='2',OMP_WAIT_POLICY='PASSIVE',OMP_PROC_BIND='close',OMP_PLACES='cores',I_MPI_PIN='1',I_MPI_PIN_DOMAIN='omp',I_MPI_PIN_PROCESSOR_LIST=','.join(map(str,cores)),I_MPI_JOB_RESPECT_PROCESS_PLACEMENT='0',I_MPI_HYDRA_BOOTSTRAP='fork')
  cmd=['taskset','-c',','.join(map(str,cores)),'/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun','-bootstrap','fork','-np','4','-ppn','4',str(O/'bin/GIZMO'),'params.txt']+(['2'] if nx==128 else [])
  write(Q/'command.json',cmd);write(Q/'environment.json',{k:v for k,v in env.items() if k.startswith(('OMP_','I_MPI_','LD_LIBRARY'))});tick=time.monotonic()
  with (Q/'log').open('x') as f:
   p=subprocess.Popen(cmd,cwd=Q,env=env,stdout=f,stderr=f);state('running',child_pid=p.pid,start=plan['start'],target=10.);code=p.wait()
  text=(Q/'log').read_text();snaps=[]
  for f in sorted((Q/'output').glob('snap_*.hdf5')):
   with h5py.File(f) as h:snaps.append(dict(path=str(f),time=float(h['Header'].attrs['Time']),cells=int(h['Header'].attrs['NumPart_Total'][0])))
  passed=bool(code==0 and 'Final time=10 reached.' in text and any(abs(x['time']-10)<1e-10 for x in snaps))
  write(Q/'result.json',dict(passed=passed,exit=code,wall_seconds=time.monotonic()-tick,snapshots=snaps,rank_count=4,restart_scope=plan['restart_scope']));assert passed
  state('completed',saved_time=10.,target=10.,await_next_stage_approval=True)
 else:raise ValueError(kind)
except Exception as exc:
 if Q.exists():state('failed',error=repr(exc),no_automatic_physical_changes=True)
 traceback.print_exc();sys.exit(2)
