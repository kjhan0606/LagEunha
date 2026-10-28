"""Continue the instrumented checkpoint; preserve every prior run."""
import os,json,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_orbit_work_20261018/on'
root=b/'Codex_orbit_work_interval_20261018';root.mkdir(exist_ok=False)
_,t=read(src/'khout.000225.dat')
for f in ['params.dat','khout.000225.dat']:shutil.copyfile(src/f,root/f)
c=json.loads((src/'environment.json').read_text());c.update(HYDRO_TSTOP=str(t+.09875))
(root/'environment.json').write_text(json.dumps(c,indent=2))
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
with (root/'log').open('x') as f:
 p=subprocess.Popen(['timeout','1200','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','225'],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
 meta={'pid':p.pid,'source':str(src/'khout.000225.dat'),'start_time':t,'expected_steps':40,'expected_end':t+.1,'binary':str(b/'Codex_orbit_work_20261018/eunha2.exe')}
 (root/'manifest.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
 code=p.wait()
(root/'exit.code').write_text(str(code));print('exit',code,flush=True)
