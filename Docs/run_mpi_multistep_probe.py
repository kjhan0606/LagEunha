"""Eight coupled steps from preserved t10 checkpoint, serial vs MPI4."""
import hashlib,json,os,shutil,subprocess
from pathlib import Path
source=Path('/gpfs/kjhan/LagForce/Codex_mpi_coupled_step_20261017')
root=Path('/gpfs/kjhan/LagForce/Codex_mpi_multistep_20261018')
root.mkdir(exist_ok=False)
binary=source/'eunha2.exe'
expected='fc2d1b5168dee1e531c44c2bfd77960f31f4e11e8f10680faa8b68679a25334e'
assert hashlib.sha256(binary.read_bytes()).hexdigest()==expected
jobs=[]
for case in ['serial','np4']:
    run=root/case;run.mkdir()
    for name in ['params.dat','restart.step']:
        shutil.copyfile(source/case/name,run/name)
    (run/'khout.002305.dat').symlink_to(source/case/'khout.002305.dat')
    env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
    config=json.loads((source/case/'environment.json').read_text())
    config['HYDRO_TSTOP']='10.0058167693911'
    config['EUNHA_DUMP_DT']='1000000'
    config['EUNHA_DUMP_FINAL']='1'
    env.update(config)
    (run/'environment.json').write_text(json.dumps(config,indent=2)+'\n')
    (run/'binary.sha256').write_text(expected+'\n')
    log=(run/'log').open('x')
    command=['timeout','600','mpirun','-np',config['PROBE_RANKS'],str(binary),'params.dat','2305']
    proc=subprocess.Popen(command,cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT)
    (run/'pid').write_text(str(proc.pid)+'\n')
    jobs.append((case,run,proc,log))
    print(json.dumps({'case':case,'pid':proc.pid,'command':command}),flush=True)
for case,run,proc,log in jobs:
    code=proc.wait();log.close();(run/'exit.code').write_text(str(code)+'\n')
    print(json.dumps({'case':case,'exit_code':code}),flush=True)
