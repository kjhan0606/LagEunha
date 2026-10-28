"""Initial opt-in production reconstruction checks on preserved t3.4 states."""
import json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import sha
b=Path('/gpfs/kjhan/LagForce');build=b/'Codex_smooth_qr_build_20261018'
assert (build/'build.exit').read_text().strip()=='0'
src=b/'Codex_coupled_response_forces_t3p4_20261018';root=b/'Codex_smooth_qr_probes_20261018'
assert shutil.disk_usage(b).free>4*1024**3
root.mkdir(exist_ok=False);jobs=[];manifest={'binary_sha256':sha(build/'eunha2.exe'),'cases':{}}
for name,sign,enabled in [('off_plus','plus',0),('on_plus','plus',1),('on_minus','minus',1)]:
 d=root/name;d.mkdir()
 for f in ['params.dat','khout.000743.dat']:shutil.copyfile(src/sign/f,d/f)
 c=json.loads((src/sign/'environment.json').read_text());c.update(GFS_QR_SMOOTH=str(enabled),GFS_SHF_DUMP_FILE=str(d/'raw_faces.dat'),OMP_NUM_THREADS='1')
 if enabled:c['GFS_QR_OPERATOR_DUMP']=str(d/'qr_operator.bin')
 (d/'environment.json').write_text(json.dumps(c,indent=2))
 env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(c)
 f=(d/'log').open('x');p=subprocess.Popen(['timeout','600','mpirun','-np','1',str(build/'eunha2.exe'),'params.dat','743'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 jobs.append((name,d,p,f));manifest['cases'][name]={'sign':sign,'enabled':enabled,'pid':p.pid,'source_sha256':sha(d/'khout.000743.dat')};print(name,p.pid,flush=True)
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
codes=[]
for name,d,p,f in jobs:
 code=p.wait();f.close();(d/'exit.code').write_text(str(code));codes.append(code);print(name,code,flush=True)
assert all(c==0 for c in codes)
