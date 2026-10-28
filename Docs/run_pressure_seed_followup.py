"""Continue pressure seeds; separate matched half-dt comparison."""
import json,os,subprocess,shutil
from pathlib import Path
b=Path('/gpfs/kjhan/LagForce');old=b/'Codex_pressure_seed_20261018';root=b/'Codex_pressure_followup_20261018';root.mkdir(exist_ok=False)
jobs=[];meta={}
for name,source,step,dt,stop in [
 ('resume_plus',old/'plus',2313,.0005,10.0218167693911),
 ('resume_minus',old/'minus',2313,.0005,10.0218167693911),
 ('fine_plus',old/'plus',2305,.00025,10.0059417693911),
 ('fine_minus',old/'minus',2305,.00025,10.0059417693911),
 ('fine_reference',b/'Codex_entropy_multistep_20261018/off',2305,.00025,10.0059417693911)]:
 assert (source/'exit.code').read_text().strip()=='0'
 d=root/name;d.mkdir();shutil.copyfile(source/'params.dat',d/'params.dat');snap=f'khout.{step:06d}.dat';(d/snap).symlink_to(source/snap)
 config=json.loads((source/'environment.json').read_text());config.update(GFS_MAX_DT=str(dt),HYDRO_TSTOP=str(stop));(d/'environment.json').write_text(json.dumps(config,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(config)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','1200','mpirun','-np','1',str(b/'Codex_entropy_split_20261018/eunha2.exe'),'params.dat',str(step)],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((d,p,f));meta[name]={'source':str(source/snap),'restart_step':step,'dt':dt,'pid':p.pid}
(root/'manifest.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(d.name,code,flush=True)
