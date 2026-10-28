"""Preserved-state physical comparison across the historical failure interval."""
import json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');old=b/'Codex_coupled_onset_to_t4_20261018/reference';seed=b/'Codex_coupled_onset_seed_20261018/reference'
root=b/'Codex_smooth_qr_physical_to_t10_20261018';target=10.002905917763403
assert (old/'exit.code').read_text().strip()=='0'
assert (b/'Codex_smooth_qr_steps_20261018/generator_comparison.json').is_file()
spec=[('legacy_continue',old,984,b/'Codex_orbit_work_20261018/eunha2.exe',0),
      ('smooth_from_t2',seed,184,b/'Codex_smooth_qr_build_20261018/eunha2.exe',1)]
prepared=[];estimate=0
for name,src,step,binary,on in spec:
 p=src/f'khout.{step:06d}.dat';a,t=read(p);assert t<target
 c=json.loads((old/'environment.json').read_text())
 c.update(HYDRO_TSTOP=str(target),GFS_MAX_DT='.0025',GFS_QR_SMOOTH=str(on),OMP_NUM_THREADS='4',EUNHA_DUMP_DT='2',EUNHA_DUMP_FINAL='1')
 for key in ['GFS_TCF_DUMP','GFS_QR_OPERATOR_DUMP','GFS_SHF_DUMP_FILE','GFS_TEST_REJECT_FINAL_ONCE']:
  assert key not in c
 assert not int(c.get('GFS_MPI_GEOMETRY_ONLY','0'))
 dumps=int((target-t)/2)+3;estimate+=(dumps+1)*p.stat().st_size
 prepared.append((name,src,step,binary,on,p,t,c,dumps))
free=shutil.disk_usage(b).free;assert free>2*estimate
print(json.dumps({'target':target,'snapshot_interval':2,'final_dump':True,'estimated_snapshot_bytes':estimate,'free_bytes':free}),flush=True)
root.mkdir(exist_ok=False);jobs=[];meta={'target_time':target,'cases':{},'comparison_scope':'Smooth operator starts from preserved t2 reference; legacy resumes completed t4 reference from that same t2 trajectory. Adaptive final times can differ. Evaluate physical error at actual times, not short-probe norms alone.'}
for name,src,step,binary,on,p,t,c,dumps in prepared:
 d=root/name;d.mkdir();shutil.copyfile(p,d/p.name);shutil.copyfile(src/'params.dat',d/'params.dat')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');proc=subprocess.Popen(['timeout','86400','mpirun','-np','1',str(binary),'params.dat',str(step)],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((name,d,proc,f));meta['cases'][name]={'source':str(p),'source_sha256':sha(p),'source_step':step,'start_time':t,'pid':proc.pid,'binary':str(binary),'binary_sha256':sha(binary),'smooth_qr':bool(on),'threads':4,'dt_cap':.0025,'dump_upper_estimate':dumps}
 (root/'manifest.json').write_text(json.dumps(meta,indent=2));print(name,proc.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
