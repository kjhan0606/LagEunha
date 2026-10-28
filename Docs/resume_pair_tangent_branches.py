"""Continue completed branches without reapplying the initial perturbation."""
import os,json,shutil,subprocess
from pathlib import Path
base=Path('/gpfs/kjhan/LagForce');old=base/'Codex_pair_tangent_branches_20261018'
root=base/'Codex_pair_tangent_resume_20261018';root.mkdir(exist_ok=False)
jobs=[];manifest={'restart_step':2313,'additional_steps_requested':32,'sources':{}}
for case in ['reference','plus','minus','half_plus','half_minus']:
 src=base/'Codex_entropy_multistep_20261018/off' if case=='reference' else old/case
 assert (src/'exit.code').read_text().strip()=='0'
 d=root/case;d.mkdir();shutil.copyfile(src/'params.dat',d/'params.dat');(d/'khout.002313.dat').symlink_to(src/'khout.002313.dat')
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
 config=json.loads((src/'environment.json').read_text());config['HYDRO_TSTOP']='10.0218167693911';env.update(config)
 (d/'environment.json').write_text(json.dumps(config,indent=2));f=(d/'log').open('x')
 p=subprocess.Popen(['timeout','1200','mpirun','-np','1',str(base/'Codex_entropy_split_20261018/eunha2.exe'),'params.dat','2313'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((d,p,f));manifest['sources'][case]={'path':str(src/'khout.002313.dat'),'pid':p.pid}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest),flush=True)
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code)+'\n');print(d.name,code,flush=True)
