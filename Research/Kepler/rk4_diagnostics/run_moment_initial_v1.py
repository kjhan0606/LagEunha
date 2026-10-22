"""Native CUDA moment initialization gates only; no Python evolution."""
from pathlib import Path
import hashlib
import argparse
import json
import os
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from audit_snapshot import read_snapshot

parser=argparse.ArgumentParser()
parser.add_argument('--version',default='v1')
parser.add_argument('--skip-patches',action='store_true')
args=parser.parse_args()
dest=ROOT/('rk4_diagnostics/moment_initial_'+args.version)
dest.mkdir(exist_ok=False)
binary=ROOT/('rk4_diagnostics/moment_initial_'+args.version+'_build/eunha2_moment_initial.exe')
assert (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
manifest=dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=[],
    output_policy='initialization/realizability only; native full-stride snapshots plus moment audits, <=0.4GiB; no time evolution')

def save():
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def execute(name,initial,params,affine=None,ranks=1,step=1,legacy_guard=False):
    d=dest/name;d.mkdir()
    _,meta=read_snapshot(initial)
    (d/f'khout.{step:06d}.dat').symlink_to(initial)
    (d/'params.dat').symlink_to(params)
    (d/'eunha2.exe').symlink_to(binary)
    (d/'OWNER.txt').write_text('Codex native moment initialization gate; no evolution; preserve data.\n')
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('GFS_'):env.pop(key)
    env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS=str(ranks),I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        EUNHA_IC='affine_probe' if affine is not None else 'kepler',KEPLER_RESTART_STEP=str(step),
        HYDRO_TSTOP=repr(meta['t']+.001),EUNHA_DUMP_DT='100',GFS_EXACT_STOP='1',
        GFS_MOMENT_INIT_ONLY='1',GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',
        GFS_KEPLER_HALO_CELLS='8',GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_TESS_DUMP='0')
    if affine is not None:env['GFS_MOMENT_AFFINE']=','.join(map(repr,affine))
    if legacy_guard:env.pop('GFS_MOMENT_INIT_ONLY')
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','EUNHA_','KEPLER_','NATIVE_'))},indent=2)+'\n')
    before=hashlib.sha256(initial.read_bytes()).hexdigest()
    code=subprocess.run([sys.executable,str(ROOT/'run_native.py')],cwd=d,env=env).returncode
    log=(d/'log').read_text()
    row=dict(name=name,exit_code=code,ranks=ranks,initial=str(initial),initial_sha256=before,
        effective_params=str(d/'params.dat'),initial_unchanged=before==hashlib.sha256(initial.read_bytes()).hexdigest(),
        stopped_without_evolution='[MOMENT_INIT_ONLY]' in log,
        legacy_guard_rejected='[MOMENT_STATE] legacy restart rejected' in log,
        audit_files=len(list(d.glob('moment_initial_r*.bin'))),
        elapsed_seconds=json.loads((d/'timing.json').read_text())['elapsed_seconds'])
    manifest['cases'].append(row);save();print(json.dumps(row),flush=True)
    if not row['initial_unchanged']:raise SystemExit('Input mutated')
    if legacy_guard:
        if code==0 or not row['legacy_guard_rejected'] or row['audit_files']:
            raise SystemExit('Legacy moment-state guard failed')
    elif code or row['audit_files']!=ranks or not row['stopped_without_evolution']:
        raise SystemExit('Native moment initialization gate failed')
    return d/f'khout.{step+1:06d}.dat'

patches=ROOT/'rk4_diagnostics/ale_shadow_patches_v1'
params=patches/'ic/params.dat'
for name,coeff in ([] if args.skip_patches else [('rigid',[.2,.1,0.,-.7,.7,0.]),
                   ('shear',[.2,.1,0.,0.,-.7,0.]),
                   ('shear_boost',[3.2,-1.9,0.,0.,-.7,0.])]):
    execute(name,patches/f'ic/{name}.dat',params,[1.,1.e-6,.5,.5]+coeff)

fresh=ROOT/'core_refinement/fine_integrated/khout.000001.dat'
checkpoint=execute('fine_fresh',fresh,fresh.parent/'params.dat')
execute('fine_restart_mpi4',checkpoint,fresh.parent/'params.dat',ranks=4,step=2)
execute('legacy_guard',checkpoint,fresh.parent/'params.dat',step=2,legacy_guard=True)
manifest['classification']='native IC/moment-density/restart gates completed; no moment transport or long-time stability claim'
save()
