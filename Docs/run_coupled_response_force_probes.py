"""Dump production hydro RHS and HG modes without advancing saved states."""
import json, os, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
b=Path('/gpfs/kjhan/LagForce')
src=b/'Codex_coupled_onset_to_t4_20261018'
root=b/'Codex_coupled_response_forces_t3p4_20261018'
step=743
template=json.loads((b/'Codex_pair_hg_modes_t2_20261018/environment.json').read_text())
sources={}
for name in ['plus','minus']:
 p=src/name/f'khout.{step:06d}.dat'
 a,t=read(p)
 sources[name]=(p,t,sha(p))
assert abs(sources['plus'][1]-sources['minus'][1])<1e-12
assert shutil.disk_usage(b).free>4*1024**3
root.mkdir(exist_ok=False)
jobs=[]
for name,(p,t,digest) in sources.items():
 d=root/name;d.mkdir()
 shutil.copyfile(p,d/p.name);shutil.copyfile(src/name/'params.dat',d/'params.dat')
 c=template.copy();c.update(OMP_NUM_THREADS='4',HYDRO_TSTOP=str(t+.001),GFS_SHF_DUMP_FILE=str(d/'raw_faces.dat'))
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 (d/'manifest.json').write_text(json.dumps({'source':str(p),'source_copy':p.name,'source_sha256':digest,'time':t,'purpose':'Terminal RHS/HG-mode diagnostic; no RK or HG application'},indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x')
 proc=subprocess.Popen(['timeout','600','mpirun','-np','1',str(b/'Codex_orbit_work_20261018/eunha2.exe'),'params.dat',str(step)],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((name,d,proc,f));print(name,proc.pid,flush=True)
codes=[]
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));codes.append(code);print(name,code,flush=True)
assert all(c==0 for c in codes)
