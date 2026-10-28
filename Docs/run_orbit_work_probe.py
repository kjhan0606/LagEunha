"""Passivity and timestep closure probe from preserved t2+.1 endpoint."""
import os,json,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_orbit_work_20261018'
assert (root/'build.exit').read_text().strip()=='0'
src=b/'Codex_onset_velocity_tangent_20261018/reference'
_,t=read(src/'khout.000224.dat');jobs=[]
for name,dt,flag in [('off',.0025,0),('on',.0025,1),('half',.00125,1)]:
 d=root/name;d.mkdir(exist_ok=False)
 for f in ['params.dat','khout.000224.dat']:shutil.copyfile(src/f,d/f)
 config=json.loads((src/'environment.json').read_text())
 config.update(GFS_MAX_DT=str(dt),HYDRO_TSTOP=str(t+.001875),GFS_ORBIT_WORK_AUDIT=str(flag))
 (d/'environment.json').write_text(json.dumps(config,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(config)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(root/'eunha2.exe'),'params.dat','224'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((name,p,f,d));print(name,p.pid,flush=True)
for name,p,f,d in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
