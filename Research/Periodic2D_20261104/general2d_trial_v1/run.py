from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,os,time,re,hashlib
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=Path(__file__).resolve().parent
preflight=pkg/'general2d_contact_preflight_v1'
exe=pkg/'general2d_build_v1/unified_native_general.exe';maker=pkg/'general2d_initial_v1/initial'
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
tasks=[('contact',32,'axis',.25),('contact',64,'axis',.25),
       ('moving_contact',32,'axis',.5),('moving_contact',64,'axis',.5),
       ('moving_contact',32,'diagonal',.5),('moving_contact',64,'diagonal',.5),
       ('shear',32,'axis',.25),('shear',64,'axis',.25),
       ('sod',64,'axis',.2),('sod',128,'axis',.2),('sod',64,'diagonal',.15)]
(root/'registration.json').write_text(json.dumps({'purpose':'Full native 2D periodic contacts, steady shear, and planar shocks using one common closure',
    'cases':tasks,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
    'new_outputs':'one exact final snapshot per case plus new IC when needed',
    'expected_new_storage_bytes':150000000,'matched_speed_benchmark':False},indent=2)+'\n')
def run(task):
    problem,n,angle,target=task;name=f'{problem}_N{n}_{angle}';d=root/name;d.mkdir(exist_ok=False)
    if task==tasks[0]:
        source=next(f for f in preflight.glob('khout.*.dat') if not f.is_symlink())
        step=int(source.name.split('.')[1]);(d/f'khout.{step:06d}.dat').symlink_to(source)
        initial=preflight/'initial.bin'
    else:
        step=1;initial=d/'initial.bin'
        r=subprocess.run([str(maker),str(n),'2',problem,angle,'1.4',str(initial)],text=True,capture_output=True)
        (d/'initial.json').write_text(r.stdout)
        if r.returncode:raise RuntimeError(r.stderr)
        (d/'khout.000001.dat').symlink_to(initial)
    params=(preflight/'params.dat').read_text()
    params=re.sub(r'(define Nx\s*=) 32',rf'\1 {n}',params)
    params=re.sub(r'(define Ny\s*=) 32',rf'\1 {n}',params)
    (d/'params.dat').write_text(params)
    env=json.loads((preflight/'environment.json').read_text());env.update({'EUNHA_IC':problem,'HYDRO_TSTOP':str(target)})
    (d/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    cmd=[str(exe),'params.dat',str(step)];start=time.monotonic()
    with (d/'log').open('x') as out:
        process=subprocess.Popen(cmd,cwd=d,env={**base,**env},stdout=out,stderr=subprocess.STDOUT)
        (d/'process.json').write_text(json.dumps({'pid':process.pid,'command':cmd,'initial':str(initial)},indent=2)+'\n')
        rc=process.wait()
    row={'name':name,'problem':problem,'Nside':n,'angle':angle,'target':target,'exit':rc,'wall_seconds':time.monotonic()-start}
    (d/'result.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True);return row
with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,tasks))
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
