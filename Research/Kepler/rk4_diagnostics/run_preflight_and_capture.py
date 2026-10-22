"""Verify diagnostic noninterference, then capture first raw-Sedov failure."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
args=argparse.ArgumentParser()
args.add_argument('--run-dir',default='validation')
args.add_argument('--binary',default='v2/eunha2_rk4_diagnostics.exe')
args.add_argument('--preflight-only',action='store_true')
args=args.parse_args()
base=root/'rk4_diagnostics';dest=base/args.run_dir;dest.mkdir(exist_ok=False)
new=base/args.binary
old=root/'linear_velocity_limiter/local_box/eunha2_local_box.exe'
source=root/'linear_velocity_limiter/physics_v4/sod_x_vector_raw'
manifest={'new_binary_sha256':hashlib.sha256(new.read_bytes()).hexdigest(),'cases':[]}
def run(label,source,binary,extra):
    d=dest/label;d.mkdir()
    for name in ('params.dat','khout.000001.dat'):(d/name).symlink_to(source/name)
    (d/'eunha2.exe').symlink_to(binary)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',
        I_MPI_HYDRA_BOOTSTRAP='fork',EUNHA_IC='shock_probe',KEPLER_RESTART_STEP='1',
        HYDRO_TSTOP='1.e-12',EUNHA_DUMP_DT='1.e-9',GFS_GPU_FINAL_VOLUME='1',
        GFS_VELOCITY_LIMITER='1',GFS_GPU_RAW_FACES='1',GFS_GPU_FACE_DUMP_CALL='0')
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT',
                'GFS_FAIL_NONPOSITIVE_U','GFS_RK4_FAILURE_SNAPSHOT','GFS_EXACT_STOP',
                'GFS_RK4_AUDIT_STEP_MIN','GFS_RK4_AUDIT_STEP_MAX'):
        env.pop(key,None)
    env.update(extra)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','KEPLER_','EUNHA_'))},indent=2)+'\n')
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    manifest['cases'].append({'name':label,'exit_code':code})
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return d,code
reference,code=run('reference_default',source,old,{})
if code:raise RuntimeError('unchanged reference failed')
for label,env in (('new_default',{}),('new_strict_audit',{'GFS_FAIL_NONPOSITIVE_U':'1','GFS_RK4_AUDIT_STEP_MIN':'1'})):
    d,code=run(label,source,new,env)
    if code:raise RuntimeError(f'{label} failed')
    a,ma=read_snapshot(reference/'khout.000002.dat');b,mb=read_snapshot(d/'khout.000002.dat')
    errors={field:float(np.max(abs(a[field].astype(float)-b[field].astype(float)))) for field in a.dtype.names if field!='u4if'}
    if any(errors.values()) or not np.array_equal(a['u4if'],b['u4if']) or ma!=mb:
        raise RuntimeError(f'{label} changed native state: {errors}')
    manifest[label]={'all_recorded_fields_bitwise_equal':True,'meta':mb}
d,code=run('new_exact_endpoint',source,new,{'GFS_FAIL_NONPOSITIVE_U':'1','GFS_EXACT_STOP':'1','HYDRO_TSTOP':'1.e-6','EUNHA_DUMP_DT':'1.e-6'})
if code:raise RuntimeError('exact endpoint failed')
p,m=read_snapshot(d/'khout.000002.dat')
assert m['t']==1.e-6 and m['dt']==1.e-6,m
manifest['exact_endpoint']={'t':m['t'],'dt':m['dt'],'pass':True}
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('default diagnostic noninterference and exact endpoint passed',flush=True)
if args.preflight_only:raise SystemExit(0)
source=root/'linear_velocity_limiter/physics_v5/sedov_legacy_raw'
d,code=run('sedov_first_failure',source,new,{'GFS_FAIL_NONPOSITIVE_U':'1','GFS_RK4_FAILURE_SNAPSHOT':'1',
    'GFS_RK4_AUDIT_STEP_MIN':'440','GFS_RK4_AUDIT_STEP_MAX':'470',
    'GFS_VELOCITY_LIMITER':'0','HYDRO_TSTOP':'.05','EUNHA_DUMP_DT':'1000'})
snapshots=sorted(d.glob('khdiag.step*.start.dat'))
if code==0 or len(snapshots)!=1:raise RuntimeError(f'first-failure diagnostic did not capture expected invalid state: {code}')
p,m=read_snapshot(snapshots[0]);assert np.all(p['ie']>0)
manifest['sedov_failure']={'classification':'expected diagnostic abort; physics gate remains FAILED',
    'exit_code':code,'valid_step_start_snapshot':str(snapshots[0]),'start_time':m['t'],
    'negative_cells_in_start':int(np.sum(p['ie']<=0))}
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('captured first Sedov invalid stage',snapshots[0],m['t'],flush=True)
