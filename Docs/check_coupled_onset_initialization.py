"""Actual solver initialization probes for prepared coupled branches."""
import os,json,shutil,subprocess,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_coupled_onset_seed_20261018';jobs=[]
template=json.loads((b/'Codex_pair_hg_modes_t2_20261018/environment.json').read_text())
for name in ['reference','plus','minus']:
 d=root/(name+'_probe');d.mkdir(exist_ok=False)
 for f in ['params.dat','khout.000184.dat']:shutil.copyfile(root/name/f,d/f)
 c=template.copy();c.update(OMP_NUM_THREADS='4',GFS_SHF_DUMP_FILE=str(d/'raw_faces.dat'))
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat','184'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((name,d,p,f));print(name,p.pid,flush=True)
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));assert code==0
out={}
for name,d,p,f in jobs:
 a,t=read(d/'khout.000184.dat');s=np.loadtxt(d/'raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
 assert np.array_equal(s[:,0],a['id'] & 72057594037927935)
 mask=np.hypot(a['x']-3,a['y']-3)<1.9
 errs={field:float(np.max(abs(s[mask,col]/a[field][mask]-1))) for field,col in [('volume',6),('pressure',5),('mass',7)]}
 for v in errs.values():assert v<1e-9
 for field,col in [('x',1),('y',2),('vx',3),('vy',4)]:assert np.array_equal(s[:,col],a[field])
 out[name]={'max_relative_initialization_errors_r_lt_1p9':errs,'positions_velocities_exact':True,'exit':code}
with (root/'initialization_check.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
