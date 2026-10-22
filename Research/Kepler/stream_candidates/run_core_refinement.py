"""Controlled native evolution, sharing one integrated fine IC across limiters."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parent;sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
args=argparse.ArgumentParser()
args.add_argument('--run-dir',default='native')
args.add_argument('--binary',default='rk4_diagnostics/v3/eunha2_rk4_diagnostics.exe')
args.add_argument('--preflight-only',action='store_true')
args.add_argument('--face-dump',action='store_true')
args=args.parse_args()
base=root/'core_refinement';dest=base/args.run_dir;dest.mkdir(exist_ok=False)
binary=root/args.binary
assert (root/'rk4_diagnostics/validation_v3/manifest.json').exists(),'diagnostic default/exact endpoint gate must finish first'
manifest={'purpose':'fixed-eps integrated-IC core resolution test, native C/CUDA only',
          'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
          'stop':.2,'eps':.047,'Re':'infinity','steering':0,'weights':0,
          'output_policy':'one final native checkpoint per successful branch; on failure retain normal step-start and invalid stage, no timestep dump series',
          'cases':[]}
def prepare(label,ic,velocity,stop=.2,first_rhs=True):
    d=dest/label;d.mkdir()
    for name in ('params.dat','khout.000001.dat'):(d/name).symlink_to(base/ic/name)
    (d/'eunha2.exe').symlink_to(binary)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',
        I_MPI_HYDRA_BOOTSTRAP='fork',KEPLER_RESTART_STEP='1',EUNHA_IC='kepler',
        HYDRO_TSTOP=str(stop),EUNHA_DUMP_DT=str(stop),GFS_EXACT_STOP='1',
        GFS_FAIL_NONPOSITIVE_U='1',GFS_RK4_FAILURE_SNAPSHOT='1',
        GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',GFS_VELOCITY_LIMITER=str(velocity),
        GFS_GPU_FACE_DUMP_CALL='1' if args.face_dump else '0',GFS_GPU_TESS_DUMP='0')
    for k in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT','GFS_DIAG_DT_MAX',
              'GFS_AFFINE_UNLIMITED_VELOCITY','GFS_AFFINE_UNMERGED_FACES','GFS_RK4_AUDIT_STAGE',
              'GFS_RK4_AUDIT_STEP_MIN','GFS_RK4_AUDIT_STEP_MAX'):env.pop(k,None)
    if first_rhs:env.update(GFS_RK4_AUDIT_STEP_MIN='1',GFS_RK4_AUDIT_STEP_MAX='1',GFS_RK4_AUDIT_STAGE='k1_rhs')
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','EUNHA_','HYDRO_','KEPLER_'))},indent=2)+'\n')
    return d,env,{'name':label,'initial_data':ic,'velocity_limiter':velocity,'stop':stop}
def execute(item):
    d,env,row=item
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    return dict(row,exit_code=code,classification='completed short gate' if code==0 else 'FAILED; preserved diagnostic state')
pre=prepare('preflight_fine','fine_integrated',0,1e-7)
row=execute(pre)
if row['exit_code']:raise RuntimeError('fine geometry/positivity preflight failed; do not launch longer branches')
d=pre[0];p,m=read_snapshot(d/'khout.000002.dat')
assert m['t']==1e-7 and len(p)==71080
initial,_=read_snapshot(d/'khout.000001.dat')
states=np.genfromtxt(d/'rk4_stage_r0.tsv',names=True,dtype=None,encoding=None,delimiter='\t')
order=np.argsort(states['id'])
expected=np.load(base/'fine_integrated/geometry.npz')['volume']
relative=abs(states['V'][order]-expected)/expected
assert relative.max()<2e-6,relative.max()
assert np.all(states['U_raw']>0) and np.all(states['mass']>0)
assert 'halo_search_clipped=1' not in (d/'log').read_text()
manifest['preflight']={'native_first_force_volume_relative_error_max':float(relative.max()),
                       'compared_cells':len(p),'minimum_native_first_volume':float(states['V'].min()),
                       'time':m['t'],'pass':True}
if args.face_dump:
    subprocess.run([sys.executable,str(root/'audit_face_reciprocity_v2.py'),'--call','1',
        '--output',str(d/'reciprocity.json'),str(d)],check=True)
    geometry=next(iter(json.loads((d/'reciprocity.json').read_text()).values()))
    roundoff=512*np.finfo(float).eps*24
    assert geometry['geometry_abs_error_max']<roundoff,geometry
    assert geometry['orphan_length_max']<roundoff,geometry
    assert geometry['max_cell_surface_closure']<roundoff,geometry
    manifest['preflight']['shared_geometry']=geometry
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('fine native first geometry and positive-state preflight passed',flush=True)
if args.preflight_only:raise SystemExit(0)
work=[prepare('coarse_old_raw','coarse_integrated',0),prepare('fine_old_raw','fine_integrated',0),prepare('fine_vector_raw','fine_integrated',1)]
with ThreadPoolExecutor(max_workers=2) as pool:
    for future in as_completed([pool.submit(execute,item) for item in work]):
        row=future.result();manifest['cases'].append(row)
        (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(row,flush=True)
# A failed branch is evidence, not permission to extend it or repair U with a floor.
print('completed controlled group; inspect every branch classification',flush=True)
