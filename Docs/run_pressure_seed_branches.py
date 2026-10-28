"""Energy-neutral distributed pressure seed, consistent restart fields."""
import os,json,subprocess,shutil
from pathlib import Path
import numpy as np
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_pressure_seed_20261018';src=b/'Codex_entropy_multistep_20261018/off'
layout=json.loads((root/'layout.json').read_text());assert layout['size']==624
fields=['id','x','y','ie','pressure','csound','stress.K'];dtype=np.dtype({'names':fields,'formats':['<u8']+['<f8']*6,'offsets':[layout['u4if']]+[layout[k] for k in fields[1:]],'itemsize':624})
raw=(src/'khout.002305.dat').read_bytes();a=np.frombuffer(raw,dtype=dtype,count=36864,offset=16);ids=(a['id']&layout['id_mask']).astype(int)
s=np.loadtxt(b/'Codex_qbp_t10_decomposition_20261018/raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
mask=np.hypot(a['x']-3,a['y']-3)<.6;z=-float(np.float32(1.666667))*s[ids,12]/s[ids,6]
z[~mask]=0;z[mask]-=np.average(z[mask],weights=a['ie'][mask]);z/=np.max(abs(z))
assert abs(np.dot(a['ie'],z))<1e-18
np.save(root/'seed.npy',z);jobs=[];meta={'seed':'-gamma Vdot/V inside r<.6, internal-energy-weighted mean removed, max absolute normalized to one','cases':{}}
for name,eps in [('plus',1e-5),('minus',-1e-5),('half_plus',5e-6),('half_minus',-5e-6)]:
 d=root/name;d.mkdir();buf=bytearray(raw);v=np.frombuffer(buf,dtype=dtype,count=36864,offset=16);ratio=1+eps*z
 for key in ['pressure','ie','stress.K']:v[key]*=ratio
 v['csound']*=np.sqrt(ratio)
 assert np.all(v['pressure']>0)
 (d/'khout.002305.dat').write_bytes(buf);shutil.copyfile(src/'params.dat',d/'params.dat')
 config=json.loads((src/'environment.json').read_text());(d/'environment.json').write_text(json.dumps(config,indent=2))
 env={k:w for k,w in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(config)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(b/'Codex_entropy_split_20261018/eunha2.exe'),'params.dat','2305'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);jobs.append((d,p,f))
 meta['cases'][name]={'epsilon':eps,'delta_total_U':float(np.sum(v['ie']-a['ie'])),'pid':p.pid}
(root/'manifest.json').write_text(json.dumps(meta,indent=2))
for d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));print(d.name,code,flush=True)
