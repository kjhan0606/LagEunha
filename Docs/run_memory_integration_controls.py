"""Independent integration controls, up to two simultaneous MPI1 processes."""
import concurrent.futures,json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_reversible_memory_controls_20261018'
first=b/'Codex_reversible_memory_first_step_20261018'
assert (first/'exit.code').read_text().strip()=='0'
seed=b/'Codex_coupled_onset_seed_20261018/reference'
new=b/'Codex_reversible_memory_integrated_build_20261018/eunha2.exe'
old=b/'Codex_adaptive_qr_periodic_build_20261018/eunha2.exe'
assert sha(new)=='8614d86c1714151a3ba22602fc8571a8378985ad2ca2a56a1b5e4e0aa061f5d6'
template=json.loads((first/'environment.json').read_text())
assert shutil.disk_usage(b).free>2*1024**3
root.mkdir(exist_ok=False);jobs=[];manifest={'cases':{},'scope':__doc__}
spec=[('direct_two',seed,184,new,1,4,.0001,2,False),
      ('resume',first,185,new,1,4,.0001,1,False),
      ('omp1',seed,184,new,1,1,.0001,1,False),
      ('retry',seed,184,new,1,4,.0002,1,True),
      ('off_new',seed,184,new,0,4,.0001,1,False),
      ('off_old',seed,184,old,0,4,.0001,1,False)]
for name,source,step,binary,on,threads,cap,steps,reject in spec:
    d=root/name;d.mkdir();snapshot=source/f'khout.{step:06d}.dat';_,t=read(snapshot)
    shutil.copyfile(snapshot,d/snapshot.name);shutil.copyfile(source/'params.dat',d/'params.dat')
    if name=='resume':shutil.copyfile(source/(snapshot.name+'.rmem'),d/(snapshot.name+'.rmem'))
    envcase=dict(template);envcase.update(GFS_REVERSIBLE_MEMORY=str(on),GFS_MAX_DT=str(cap),
        HYDRO_TSTOP=str(t+(steps-.5)*.0001),OMP_NUM_THREADS=str(threads))
    if name=='resume':envcase.pop('GFS_RMEM_INIT_ZERO',None)
    if reject:envcase['GFS_TEST_REJECT_FINAL_ONCE']='1'
    (d/'environment.json').write_text(json.dumps(envcase,indent=2))
    env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(envcase)
    manifest['cases'][name]={'source':str(snapshot),'source_sha256':sha(snapshot),'binary_sha256':sha(binary),
        'source_time':t,'expected_steps':steps,'final_step':step+steps,'memory_on':bool(on)}
    jobs.append((name,d,binary,step,env))
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
def execute(job):
    name,d,binary,step,env=job
    with (d/'log').open('x') as stream:
        p=subprocess.Popen(['timeout','1200','mpirun','-np','1',str(binary),'params.dat',str(step)],cwd=d,env=env,stdout=stream,stderr=subprocess.STDOUT)
        (d/'pid').write_text(str(p.pid));print(name,p.pid,flush=True);code=p.wait()
    (d/'exit.code').write_text(str(code));print(name,'exit',code,flush=True)
    return name,code
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    outcomes=list(pool.map(execute,jobs))
assert all(code==0 for _,code in outcomes),outcomes
