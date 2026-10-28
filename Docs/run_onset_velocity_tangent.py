"""Perturb the naturally evolved early velocity departure at t2."""
import sys,json,os,subprocess,shutil
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');src=b/'Codex_early_growth_budget_20261017/t2_half';root=b/'Codex_onset_velocity_tangent_20261018';root.mkdir(exist_ok=False)
layout=json.loads((b/'Codex_DPQ_backend_parity_20261017/layout.json').read_text());raw=(src/'khout.000184.dat').read_bytes()
fields=['id','x','y','vx','vy','mass'];dtype=np.dtype({'names':fields,'formats':['<u8']+['<f8']*5,'offsets':[layout['u4if']]+[layout[k] for k in fields[1:]],'itemsize':layout['size']})
a=np.frombuffer(raw,dtype=dtype,count=36864,offset=16);_,time=read(src/'khout.000184.dat')
x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1);vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]));seed=np.column_stack((a['vx'],a['vy']))-vc
mask=r<.6;seed[~mask]=0;seed[mask]-=np.average(seed[mask],axis=0,weights=a['mass'][mask]);np.save(root/'seed.npy',seed)
jobs=[];meta={'start_time':time,'seed':'Actual v-v_circular inside r<.6, mass-weighted translation removed. Velocities alone perturbed at restart, no steering later.','cases':{}}
for name,eps in [('reference',0.),('plus',1e-3),('minus',-1e-3)]:
 d=root/name;d.mkdir();buf=bytearray(raw);q=np.frombuffer(buf,dtype=dtype,count=36864,offset=16);q['vx']+=eps*seed[:,0];q['vy']+=eps*seed[:,1]
 (d/'khout.000184.dat').write_bytes(buf);shutil.copyfile(src/'params.dat',d/'params.dat')
 config=json.loads((src/'environment.json').read_text());config.update(GFS_MAX_DT='0.0025',HYDRO_TSTOP=str(time+.09875));(d/'environment.json').write_text(json.dumps(config,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(config)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','1200','mpirun','-np','1',str(b/'Codex_entropy_split_20261018/eunha2.exe'),'params.dat','184'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((d,p,f));meta['cases'][name]={'epsilon':eps,'pid':p.pid}
(root/'manifest.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(d.name,code,flush=True)
