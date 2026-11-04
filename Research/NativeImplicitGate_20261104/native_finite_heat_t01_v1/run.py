from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,os,json,time,hashlib
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1')
root=Path(__file__).resolve().parent
exe=pkg/'native_coupled_build_v5/unified_native_coupled.exe'
old=pkg/'native_coupled_t0_short_v3'
baseenv={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
env_template=json.loads((old/'environment.json').read_text())
def branch(name,factor,source):
    results=[]
    for target in [.01,.03,.1]:
        d=root/f'{name}_to{target:g}';d.mkdir(exist_ok=False)
        (d/'params.dat').write_text((old/'params.dat').read_text())
        # Restart step is named by the previous accepted writer, not guessed.
        step=1 if name=='half' and not results else int(source.name.split('.')[1])
        (d/f'khout.{step:06d}.dat').symlink_to(source)
        env={**env_template,'HYDRO_TSTOP':repr(target),'GFS_KDK_MAX_WINDOW':'.2',
             'GFS_UNIFIED_DT_FACTOR':str(factor),'GFS_UNIFIED_COHORT_FILE':str(root/'original_cohorts.bin')}
        (d/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
        cmd=[str(exe),'params.dat',str(step)]
        start=time.monotonic()
        with (d/'log').open('x') as out:
            process=subprocess.Popen(cmd,cwd=d,env={**baseenv,**env},stdout=out,stderr=subprocess.STDOUT)
            (d/'process.json').write_text(json.dumps({'pid':process.pid,'command':cmd,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()},indent=2)+'\n')
            rc=process.wait()
        row={'branch':name,'dt_factor':factor,'target':target,'exit':rc,'wall_seconds':time.monotonic()-start}
        (d/'result.json').write_text(json.dumps(row,indent=2)+'\n')
        results.append(row);print(json.dumps(row),flush=True)
        if rc:break
        outputs=[f for f in d.glob('khout.*.dat') if not f.is_symlink()]
        if len(outputs)!=1:raise RuntimeError(f'Expected one final snapshot in {d}')
        source=outputs[0]
    return results
if __name__=='__main__':
    initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_w0_initial_20261101_v1/initial696.bin')
    tasks=[('base',1.,old/'khout.000003.dat'),('half',.5,initial)]
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(branch,*task) for task in tasks]
        results=[future.result() for future in futures]
    (root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
