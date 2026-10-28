"""Eight t10 steps with entropy component diagnostics off/on."""
import hashlib,json,os,shutil,subprocess
from pathlib import Path
source=Path('/gpfs/kjhan/LagForce/Codex_mpi_coupled_step_20261017')
root=Path('/gpfs/kjhan/LagForce/Codex_entropy_multistep_20261018')
root.mkdir(exist_ok=False)
binary=Path('/gpfs/kjhan/LagForce/Codex_entropy_split_20261018/eunha2.exe')
expected=Path('/gpfs/kjhan/LagForce/Codex_entropy_split_20261018/binary.sha256').read_text().strip()
assert hashlib.sha256(binary.read_bytes()).hexdigest()==expected
jobs=[]
for case in ['off','on']:
    run=root/case;run.mkdir()
    for name in ['params.dat','restart.step']:
        shutil.copyfile(source/'serial'/name,run/name)
    (run/'khout.002305.dat').symlink_to(source/'serial'/'khout.002305.dat')
    env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
    config=json.loads((source/'serial'/'environment.json').read_text())
    config['GFS_TCF_SACC']='entropy.dat'
    if case=='on':config['GFS_ENTROPY_SPLIT']='1'
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
