from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,json,subprocess,time,re,hashlib,shutil
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=Path(__file__).resolve().parent
exe=pkg/'solid_angle_build_v1/unified_native_general.exe';maker=pkg/'general2d_initial_v1/initial'
tasks=[('contact',n,'axis',.25) for n in [32,64]]+ [('moving_contact',n,angle,.5) for n in [32,64] for angle in ['axis','diagonal']]+[('shear',n,'axis',.25) for n in [32,64,128]]+[('sod',n,'axis',.2) for n in [64,128,256]]+[('sod',n,'diagonal',.15) for n in [64,128]]
free=shutil.disk_usage(root).free
(root/'registration.json').write_text(json.dumps({'purpose':'actual native 2D test of one solid-angle LSQ correction to all cases; preserve prior trajectories','cases':tasks,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'new_output_policy':'one final snapshot per case; reuse prior ICs; summaries for axis shocks','estimated_output_bytes':160000000,'free_bytes':free,'matched_speed_benchmark':False},indent=2)+'\n')
assert free>500000000
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))}
source=pkg/'general2d_trial_v1/sod_N64_axis'
def run(task):
 problem,n,angle,target=task;name=f'{problem}_N{n}_{angle}';d=root/name;d.mkdir(exist_ok=False)
 original=pkg/'general2d_trial_v1'/name/'initial.bin'
 if problem=='contact' and n==32:original=pkg/'general2d_contact_preflight_v1/initial.bin'
 if problem=='shear' and n==128:original=pkg/'general2d_shear_corrected_v1/N128/initial.bin'
 if original.exists():initial=original
 else:
  initial=d/'initial.bin';r=subprocess.run([str(maker),str(n),'2',problem,angle,'1.4',str(initial)],text=True,capture_output=True,check=True);(d/'initial.json').write_text(r.stdout)
 (d/'khout.000001.dat').symlink_to(initial)
 params=(source/'params.dat').read_text();params=re.sub(r'(define N[xy]\s*=) 64',rf'\1 {n}',params);(d/'params.dat').write_text(params)
 env=json.loads((source/'environment.json').read_text());env.update({'EUNHA_IC':problem,'HYDRO_TSTOP':str(target)})
 if problem=='sod' and angle=='axis':env['GFS_UNIFIED_PLANAR_AUDIT_N']=str(n)
 (d/'environment.json').write_text(json.dumps(env,indent=2)+'\n');t=time.monotonic()
 with (d/'log').open('x') as f:
  p=subprocess.Popen([str(exe),'params.dat','1'],cwd=d,env={**base,**env},stdout=f,stderr=subprocess.STDOUT)
  (d/'process.json').write_text(json.dumps({'pid':p.pid,'initial':str(initial)},indent=2)+'\n');rc=p.wait()
 row={'name':name,'problem':problem,'Nside':n,'angle':angle,'target':target,'exit':rc,'wall_seconds':time.monotonic()-t}
 (d/'result.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True);return row
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,tasks))
(root/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
