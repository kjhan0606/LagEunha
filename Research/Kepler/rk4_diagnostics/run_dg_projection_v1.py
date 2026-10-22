"""Fresh native Euler spatial projections. Python launches; no evolution."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from audit_snapshot import read_snapshot
parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,choices=(2,4,6),default=2)
parser.add_argument('--version',default='v1');parser.add_argument('--controls',action='store_true')
parser.add_argument('--mpi4',action='store_true');args=parser.parse_args()
dest=ROOT/f'rk4_diagnostics/dg_projection_{args.version}_p{args.order}';dest.mkdir(exist_ok=False)
binary=ROOT/f'rk4_diagnostics/dg_projection_{args.version}_build/eunha2_dg_projection.exe'
assert (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
idle=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid','--format=csv,noheader'],text=True).strip()
if idle:raise SystemExit('GPU occupied; no native launch')
manifest=dict(order=args.order,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=[],
    output_policy='read-only native projection and polygon packets, <=0.6GiB across degrees/controls/MPI check; no evolved/checkpoint output')
def save():(dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def execute(name,initial,params,affine=None,slope=None,ranks=1):
    d=dest/name;d.mkdir();_,meta=read_snapshot(initial);assert meta['t']==0
    for key,path in [('khout.000001.dat',initial),('params.dat',params),('eunha2.exe',binary)]:
        (d/key).symlink_to(path)
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('GFS_'):env.pop(key)
    env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS=str(ranks),I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        EUNHA_IC='affine_probe' if affine else 'kepler',KEPLER_RESTART_STEP='1',HYDRO_TSTOP='.001',
        EUNHA_DUMP_DT='100',GFS_DG_PROJECT_ONLY='1',GFS_DG_ORDER=str(args.order),
        GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',GFS_KEPLER_HALO_CELLS='8',
        GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_TESS_DUMP='0')
    if affine:env['GFS_DG_AFFINE']=','.join(map(repr,affine))
    if slope:env['GFS_DG_DENSITY_SLOPE']=','.join(map(repr,slope))
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','EUNHA_','HYDRO_','NATIVE_','KEPLER_'))},indent=2)+'\n')
    before=hashlib.sha256(initial.read_bytes()).hexdigest()
    code=subprocess.run([sys.executable,str(ROOT/'run_native.py')],cwd=d,env=env).returncode
    log=(d/'log').read_text()
    row=dict(name=name,exit_code=code,ranks=ranks,initial=str(initial),initial_sha256=before,
        initial_unchanged=before==hashlib.sha256(initial.read_bytes()).hexdigest(),
        effective_params=str(d/'params.dat'),stopped_without_evolution='[DG_PROJECT_ONLY]' in log,
        packets=len(list(d.glob('dg_project_r*.bin'))),
        elapsed_seconds=json.loads((d/'timing.json').read_text())['elapsed_seconds'])
    manifest['cases'].append(row);save();print(json.dumps(row),flush=True)
    if code or not row['initial_unchanged'] or not row['stopped_without_evolution'] or row['packets']!=ranks:
        raise SystemExit('Native projection execution failed; preserve all outputs')
if args.controls:
    patches=ROOT/'rk4_diagnostics/ale_shadow_patches_v1';params=patches/'ic/params.dat'
    execute('rigid',patches/'ic/rigid.dat',params,[1.,1.e-6,.5,.5,.2,.1,0.,-.7,.7,0.])
    execute('shear_boost',patches/'ic/shear_boost.dat',params,[1.,1.e-6,.5,.5,3.2,-1.9,0.,0.,-.7,0.])
    execute('varying_rho_isobaric',patches/'ic/uniform.dat',params,[1.,1.e-6,.5,.5,.2,.1,0.,0.,0.,0.],[.4,-.2])
fresh=ROOT/'core_refinement/fine_integrated/khout.000001.dat'
execute('fine_fresh',fresh,fresh.parent/'params.dat')
if args.mpi4:execute('fine_mpi4',fresh,fresh.parent/'params.dat',ranks=4)
manifest['classification']='projection records produced; analyze pressure/conservation/quadrature before any evolution'
save()
