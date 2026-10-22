"""Bounded native shadow captures; Python prepares data and never evolves it."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import struct
import subprocess
import sys

import numpy as np
from scipy.spatial import Voronoi

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from audit_snapshot import read_snapshot

parser=argparse.ArgumentParser()
parser.add_argument('--phase',choices=('patches','real'),required=True)
args=parser.parse_args()
dest=ROOT/'rk4_diagnostics'/('ale_shadow_'+args.phase+'_v1')
dest.mkdir(exist_ok=False)
binary=ROOT/'rk4_diagnostics/ale_shadow_v1_build/eunha2_ale_shadow.exe'
assert (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
manifest=dict(phase=args.phase,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    output_policy='read-only first RHS captures; about0.25GiB total across patches and two real states; no evolution or final hydro dump',
    cases=[])


def save():
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def execute(name,initial,params,affine=None,ranks=1):
    d=dest/name;d.mkdir()
    _,meta=read_snapshot(initial)
    for key,path in [('khout.000001.dat',initial),('params.dat',params),('eunha2.exe',binary)]:
        (d/key).symlink_to(path)
    (d/'OWNER.txt').write_text('Codex actual-face ALE shadow; no evolution; preserve data.\n')
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('GFS_'):env.pop(key)
    env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS=str(ranks),I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        EUNHA_IC='affine_probe' if affine is not None else 'kepler',KEPLER_RESTART_STEP='1',
        HYDRO_TSTOP=repr(meta['t']+.001),EUNHA_DUMP_DT='100',GFS_EXACT_STOP='1',
        GFS_SSPRK3='1',GFS_ALE_SHADOW_CALL='1',GFS_ALE_SHADOW_ONLY='1',
        GFS_FAIL_NONPOSITIVE_U='1',GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',
        GFS_VELOCITY_LIMITER='0',GFS_KEPLER_HALO_CELLS='8',GFS_E_AUDIT='1',GFS_FLOOR_LOG='1',
        GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_TESS_DUMP='0')
    if affine is not None:env['GFS_ALE_SHADOW_AFFINE']=','.join(map(repr,affine))
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','EUNHA_','KEPLER_','NATIVE_'))},indent=2)+'\n')
    before=hashlib.sha256(initial.read_bytes()).hexdigest()
    code=subprocess.run([sys.executable,str(ROOT/'run_native.py')],cwd=d,env=env).returncode
    files=list(d.glob('ale_meta_r*_call1.json'))
    log=(d/'log').read_text()
    row=dict(name=name,exit_code=code,ranks=ranks,initial=str(initial),initial_sha256=before,
        effective_params=str(d/'params.dat'),meta_files=len(files),
        initial_unchanged=before==hashlib.sha256(initial.read_bytes()).hexdigest(),
        stopped_without_evolution='[ALE_SHADOW_ONLY]' in log,
        elapsed_seconds=json.loads((d/'timing.json').read_text())['elapsed_seconds'])
    manifest['cases'].append(row);save();print(json.dumps(row),flush=True)
    if code or len(files)!=ranks or not row['initial_unchanged'] or not row['stopped_without_evolution']:
        raise SystemExit('Stop: shadow capture/preparation gate failed')


if args.phase=='patches':
    icdir=dest/'ic';icdir.mkdir()
    template,_=read_snapshot(ROOT/'core_refinement/fine_integrated/khout.000001.dat')
    n=64;h=1/n;rng=np.random.default_rng(20261020)
    yy,xx=np.mgrid[:n,:n]
    xy=np.column_stack([(xx.ravel()+.5)*h,(yy.ravel()+.5)*h])
    xy+=rng.uniform(-.22*h,.22*h,xy.shape)
    cloud=np.concatenate([xy+[dx,dy] for dx,dy in [(0,0)]+[(a,b) for a in (-1,0,1) for b in (-1,0,1) if a or b]])
    vor=Voronoi(cloud);areas=np.empty(n*n)
    for i in range(n*n):
        reg=vor.regions[vor.point_region[i]];assert reg and -1 not in reg
        v=vor.vertices[reg]-xy[i]
        areas[i]=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))-np.dot(v[:,1],np.roll(v[:,0],-1)))
    assert abs(areas.sum()-1)<1e-12
    gamma=float(np.float32(1.666667))
    params=icdir/'params.dat'
    text=(ROOT/'rk4_diagnostics/validation/sedov_first_failure/params.dat').read_text()
    params.write_text(text)
    variants=[('uniform',[.2,.1,0.,0.,0.,0.]),
              ('rigid',[.2,.1,0.,-.7,.7,0.]),
              ('shear',[.2,.1,0.,0.,-.7,0.]),
              ('shear_boost',[3.2,-1.9,0.,0.,-.7,0.])]
    paths={}
    for name,coeff in variants:
        p=np.zeros(n*n,dtype=template.dtype)
        p['u4if']=np.arange(n*n,dtype=np.uint64)
        p['x'],p['y']=xy.T;p['mass']=areas;p['volume']=areas
        p['den']=1;p['pressure']=1.e-6;p['ie']=1.e-6*areas/(gamma-1)
        p['vx']=coeff[0]+coeff[2]*(xy[:,0]-.5)+coeff[3]*(xy[:,1]-.5)
        p['vy']=coeff[1]+coeff[4]*(xy[:,0]-.5)+coeff[5]*(xy[:,1]-.5)
        # Derived caches are prepared natively from this normal initial state.
        path=icdir/(name+'.dat')
        with path.open('wb') as fh:
            fh.write(struct.pack('<4i',-1,len(p),p.dtype.itemsize,1));fh.write(p.tobytes())
            fh.write(struct.pack('<ddfiiddff',0.,1.e-6,1.e-6,n,n,1.,1.,0.,0.))
        affine=[1.,1.e-6,.5,.5]+coeff
        paths[name]=(path,affine)
        execute(name,path,params,affine)
    path,affine=paths['shear'];execute('shear_mpi4',path,params,affine,4)
else:
    old=Path(json.loads((ROOT/'rk4_diagnostics/fine_cold_interval/manifest.json').read_text())['source'])
    fresh=ROOT/'rk4_diagnostics/ssprk3_fresh_v2/fine_fresh/khdiag.step000516.ssp_rejected_start.dat'
    execute('fine_captured',old,old.parent/'params.dat')
    execute('fine_fresh_final',fresh,fresh.parent/'params.dat')
manifest['classification']='native shadow captures completed; analyze flux/GCL/local thermal terms before evolution'
save()
