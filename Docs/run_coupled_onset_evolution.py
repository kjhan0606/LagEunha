import os,json,subprocess
from pathlib import Path
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_onset_seed_20261018'
meta=json.loads((root/'manifest.json').read_text());check=json.loads((root/'initialization_check.json').read_text())
jobs=[];manifest={'start':meta['time'],'target':meta['time']+.1,'cases':{}}
# Claim launch before starting any case; never relaunch a partly completed batch.
with (root/'evolution_launch.json').open('x') as f:json.dump(manifest,f)
for name in ['reference','plus','minus']:
 assert check[name]['terminal_probe_marker'] and check[name]['positions_velocities_exact']
 d=root/name;c=json.loads((d/'environment.json').read_text())
 c.update(HYDRO_TSTOP=str(meta['time']+.09875),GFS_MAX_DT='.0025',OMP_NUM_THREADS='4')
 (d/'evolution_environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'evolution.log').open('x');p=subprocess.Popen(['timeout','1800','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','184'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((name,d,p,f));manifest['cases'][name]={'pid':p.pid,'threads':4}
 (root/'evolution_launch.json').write_text(json.dumps(manifest,indent=2));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'evolution_exit.code').write_text(str(code));print(name,code,flush=True)
