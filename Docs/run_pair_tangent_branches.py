"""Isolated velocity perturbations, no steering during evolution."""
import json,os,shutil,subprocess,hashlib
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_pair_tangent_branches_20261018';root.mkdir(exist_ok=False)
src=base/'Codex_entropy_multistep_20261018/off';binary=base/'Codex_entropy_split_20261018/eunha2.exe'
layout=json.loads((base/'Codex_DPQ_backend_parity_20261017/layout.json').read_text())
raw=(src/'khout.002305.dat').read_bytes();n=36864
dtype=np.dtype({'names':['id','x','y','vx','vy','mass'],'formats':['<u8']+['<f8']*5,
 'offsets':[layout['u4if']]+[layout[k] for k in ['x','y','vx','vy','mass']],'itemsize':layout['size']})
original=np.frombuffer(raw,dtype=dtype,count=n,offset=16);ids=original['id']&layout['id_mask'];assert len(np.unique(ids))==n
ii,jj=[int(np.where(ids==v)[0][0]) for v in [156,157]]
normal=np.array([original['x'][jj]-original['x'][ii],original['y'][jj]-original['y'][ii]]);normal/=np.linalg.norm(normal)
B=np.array([normal,[-normal[1],normal[0]]]);h=json.loads((base/'Codex_pair_force_hessian_controls_20261018/result.json').read_text())
H=np.array(h['samples'][-1]['models']['actual']['matrix']);_,vectors=np.linalg.eigh(.5*(H+H.T));direction=B.T@vectors[:,0]
m=original['mass'][[ii,jj]];coef=np.array([-m[1],m[0]])/m.sum()
jobs=[];manifest={'pair':[156,157],'direction':direction.tolist(),'coefficients':coef.tolist(),'reference':str(src),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'cases':{}}
for name,eps in [('plus',1e-6),('minus',-1e-6),('half_plus',5e-7),('half_minus',-5e-7)]:
 d=root/name;d.mkdir();shutil.copyfile(src/'params.dat',d/'params.dat');buf=bytearray(raw);a=np.frombuffer(buf,dtype=dtype,count=n,offset=16)
 for k,c in zip([ii,jj],coef):
  a['vx'][k]+=eps*c*direction[0];a['vy'][k]+=eps*c*direction[1]
 (d/'khout.002305.dat').write_bytes(buf)
 config=json.loads((src/'environment.json').read_text());env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(config)
 (d/'environment.json').write_text(json.dumps(config,indent=2));log=(d/'log').open('x')
 command=['timeout','600','mpirun','-np','1',str(binary),'params.dat','2305']
 proc=subprocess.Popen(command,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT);jobs.append((d,proc,log))
 manifest['cases'][name]={'epsilon':eps,'pid':proc.pid,'snapshot_sha256':hashlib.sha256(buf).hexdigest()}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code)+'\n');print(d.name,code,flush=True)
