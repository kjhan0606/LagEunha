from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,os,time,hashlib,shutil
pkg=Path('/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1');root=Path(__file__).resolve().parent
assert json.loads((pkg/'unified3d_strict_geometry_fd_v1/result.json').read_text())['exit']==0
build=pkg/'unified3d_strict_build_v3';exe=build/'unified_native_general.exe'
cases=[('contact','axis',.25),('moving_contact','axis',.5),('moving_contact','diagonal',.5),('shear','diagonal',.25),('sod','axis',.2),('sod','diagonal',.15)]
tasks=[(32,*case) for case in cases if case[0] in ['shear','sod']]
(root/'registration.json').write_text(json.dumps({'purpose':'strict actual native 3D dimensional gate, same common pressure/solid-angle LSQ/BE/KDK','cases':tasks,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'outputs':'one initial and one final LAG3DV1 snapshot per case, or accepted-state snapshot with nonzero failure','expected_new_bytes':32000000,'free_bytes':shutil.disk_usage(root).free,'GPU_or_multi_rank_qualification':False},indent=2)+'\n')
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_','LAG3D_'))}
extra=json.loads((pkg/'unified3d_contact_preflight_v2/environment.json').read_text())
extra['GFS_UNIFIED_GEOMETRY_EACH_STAGE']='1'
def run(task):
 n,problem,angle,target=task;name=f'{problem}_N{n}_{angle}';d=root/name;d.mkdir(exist_ok=False)
 if False:
  old=pkg/'unified3d_contact_preflight_v2'
  for f in ['snap_000000.l3d','snap_000001.l3d','initial.l3d']:(d/f).symlink_to(old/f)
  row={'name':name,'problem':problem,'Nside':n,'angle':angle,'target':target,'exit':0,'reused_existing_preflight':str(old)}
 else:
  r=subprocess.run([str(build/'initial3d'),str(n),'2',problem,angle,'1.4',str(d/'initial.l3d')],text=True,capture_output=True,check=True);(d/'initial.json').write_text(r.stdout)
  params=(pkg/'unified3d_contact_preflight_v2/params.dat').read_text().replace('t_end = 0.25',f't_end = {target}');(d/'params.dat').write_text(params)
  (d/'environment.json').write_text(json.dumps(extra,indent=2)+'\n');start=time.monotonic()
  with (d/'log').open('x') as f:
   p=subprocess.Popen([str(exe),'params.dat'],cwd=d,env={**base,**extra},stdout=f,stderr=subprocess.STDOUT);(d/'process.json').write_text(json.dumps({'pid':p.pid})+'\n');rc=p.wait()
  row={'name':name,'problem':problem,'Nside':n,'angle':angle,'target':target,'exit':rc,'wall_seconds':time.monotonic()-start}
 (d/'result.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True);return row
with ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,tasks))
(root/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
