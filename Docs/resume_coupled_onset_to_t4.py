"""Continue existing coupled branches through early growth, with checkpoints."""
import os,json,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_coupled_onset_seed_20261018';root=b/'Codex_coupled_onset_to_t4_20261018';root.mkdir(exist_ok=False)
meta=json.loads((src/'manifest.json').read_text());target=meta['time']+2
jobs=[];manifest={'source':str(src),'target':target,'checkpoint_interval':.2,'approximate_snapshot_bytes_per_case':11*23000000,'cases':{}}
for name in ['reference','plus','minus']:
 d=root/name;d.mkdir();a,t=read(src/name/'khout.000224.dat')
 for f in ['params.dat','khout.000224.dat']:shutil.copyfile(src/name/f,d/f)
 c=json.loads((src/name/'evolution_environment.json').read_text());c.update(HYDRO_TSTOP=str(target-.00125),EUNHA_DUMP_DT='.2',EUNHA_DUMP_FINAL='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','10800','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','224'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((name,d,p,f))
 manifest['cases'][name]={'pid':p.pid,'start_time':t,'threads':4};(root/'manifest.json').write_text(json.dumps(manifest,indent=2));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
