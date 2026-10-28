import os,json,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_orbit_work_20261018/on'
root=b/'Codex_orbit_work_robustness_20261018';root.mkdir(exist_ok=False)
_,t=read(src/'khout.000224.dat');jobs=[]
for name in ['omp4','retry']:
 d=root/name;d.mkdir()
 for f in ['params.dat','khout.000224.dat']:shutil.copyfile(src/f,d/f)
 c=json.loads((src/'environment.json').read_text())
 if name=='omp4':c['OMP_NUM_THREADS']='4'
 else:c.update(GFS_MAX_DT='.005',GFS_TEST_REJECT_FINAL_ONCE='1')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','224'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((d,p,f));print(name,p.pid,flush=True)
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(d.name,code,flush=True)
