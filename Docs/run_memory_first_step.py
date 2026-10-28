"""First actual hydro/gravity/HG/RK4 step with reversible memory enabled."""
import json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_reversible_memory_first_step_20261018'
source=b/'Codex_coupled_onset_seed_20261018/reference';snapshot=source/'khout.000184.dat'
binary=b/'Codex_reversible_memory_integrated_build_20261018/eunha2.exe'
assert sha(binary)=='8614d86c1714151a3ba22602fc8571a8378985ad2ca2a56a1b5e4e0aa061f5d6'
a,t=read(snapshot);assert shutil.disk_usage(b).free>1024**3
root.mkdir(exist_ok=False)
for filename in ['params.dat',snapshot.name]:shutil.copyfile(source/filename,root/filename)
case=json.loads((b/'Codex_adaptive_qr_physical_to_t10_20261018/adaptive_from_t2/environment.json').read_text())
case.update(GFS_REVERSIBLE_MEMORY='1',GFS_RMEM_INIT_ZERO='1',GFS_STEP_RETRY='1',
            GFS_MAX_DT='.0001',HYDRO_TSTOP=str(t+.00005),EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1')
(root/'environment.json').write_text(json.dumps(case,indent=2))
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
with (root/'log').open('x') as stream:
    job=subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat','184'],cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT)
    (root/'manifest.json').write_text(json.dumps({'pid':job.pid,'binary_sha256':sha(binary),'source_sha256':sha(snapshot),'source_time':t,'scope':'One actual step, explicit zero auxiliary state branch; source unchanged.'},indent=2))
    print('first_step',job.pid,flush=True);code=job.wait()
(root/'exit.code').write_text(str(code));print('exit',code,flush=True)
assert code==0
log=(root/'log').read_text();assert log.count('[ACCEPTED_STATE]')==1 and '[RMEM_ACCEPTED]' in log
assert sha(snapshot)==json.loads((root/'manifest.json').read_text())['source_sha256']
final=root/'khout.000185.dat';a,tf=read(final)
assert abs(tf-(t+.0001))<1e-12
assert (root/(final.name+'.rmem')).is_file()
print('accepted_time',tf,'sidecar_bytes',(root/(final.name+'.rmem')).stat().st_size,flush=True)
