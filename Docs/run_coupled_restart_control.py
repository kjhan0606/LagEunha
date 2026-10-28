"""Matched full-amplitude uninterrupted cases to isolate restart effects."""
import os,json,shutil,subprocess
from pathlib import Path
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018';root=b/'Codex_coupled_restart_control_20261018';root.mkdir(exist_ok=False)
target=json.loads((b/'Codex_coupled_convergence_20261018/manifest.json').read_text())['target_time'];jobs=[];meta={'target':target,'cases':{}}
for name in ['plus','minus']:
 d=root/name;d.mkdir()
 for f in ['params.dat','khout.000184.dat']:shutil.copyfile(src/name/f,d/f)
 c=json.loads((src/name/'evolution_environment.json').read_text());c.update(HYDRO_TSTOP=str(target-.000625),GFS_MAX_DT='.0025',EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','1800','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','184'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((name,d,p,f));meta['cases'][name]={'pid':p.pid};(root/'manifest.json').write_text(json.dumps(meta,indent=2));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
