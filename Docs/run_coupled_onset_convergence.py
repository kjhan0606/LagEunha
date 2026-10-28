"""Independent amplitude and timestep checks at the completed .1975 horizon."""
import os,json,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');full=b/'Codex_coupled_onset_seed_20261018';half=b/'Codex_coupled_halfseed_20261018';root=b/'Codex_coupled_convergence_20261018'
assert (half/'manifest.json').exists();root.mkdir(exist_ok=False)
_,target=read(b/'Codex_coupled_onset_to_t4_20261018/reference/khout.000263.dat')
jobs=[];manifest={'target_time':target,'cases':{},'reference_reuse':str(b/'Codex_coupled_onset_to_t4_20261018/reference/khout.000263.dat')}
spec=[('amplitude_plus',half/'plus',.0025,.0005),('amplitude_minus',half/'minus',.0025,-.0005),
 ('fine_reference',full/'reference',.00125,0.),('fine_plus',full/'plus',.00125,.001),('fine_minus',full/'minus',.00125,-.001)]
for name,src,dt,eps in spec:
 d=root/name;d.mkdir();a,t0=read(src/'khout.000184.dat')
 for f in ['params.dat','khout.000184.dat']:shutil.copyfile(src/f,d/f)
 c=json.loads((src/'environment.json').read_text());c.update(HYDRO_TSTOP=str(target-.000625),GFS_MAX_DT=str(dt),OMP_NUM_THREADS='4',EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','3600','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','184'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((name,d,p,f))
 manifest['cases'][name]={'pid':p.pid,'epsilon':eps,'dt':dt,'threads':4,'expected_steps':round((target-t0)/dt),'source':str(src/'khout.000184.dat')}
 (root/'manifest.json').write_text(json.dumps(manifest,indent=2));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(name,code,flush=True)
