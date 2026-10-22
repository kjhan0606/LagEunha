"""One observer-equivalence replay; no new physical parameter variant."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import numpy as np

root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot

dest=root/'rk4_diagnostics/ale_shadow_observer_v1';dest.mkdir(exist_ok=False)
old=root/'rk4_diagnostics/ssprk3_intervals_v2/fine_base'
reference=next(c for c in json.loads((old.parent/'manifest.json').read_text())['cases'] if c['name']=='fine_base')
binary=root/'rk4_diagnostics/ale_shadow_v1_build/eunha2_ale_shadow.exe'
for key,path in [('khout.000001.dat',Path(reference['initial'])),('params.dat',old/'params.dat'),('eunha2.exe',binary)]:
    (dest/key).symlink_to(path)
env=os.environ.copy()
for key in list(env):
    if key.startswith('GFS_'):env.pop(key)
env.update(json.loads((old/'environment.json').read_text()))
env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
    GFS_ALE_SHADOW_CALL='1',GFS_SSP_MAX_RHS='16',GFS_SSP_MAX_STEPS='4')
for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
(dest/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','EUNHA_','KEPLER_','NATIVE_'))},indent=2)+'\n')
(dest/'OWNER.txt').write_text('Codex bounded observer check: repeat the existing three-step interval only to verify unchanged native evolution with read-only ALE instrumentation.\n')
code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=dest,env=env).returncode
row=dict(exit_code=code,binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    scope='instrumentation observer equivalence only; existing pressure-only evolution, shadow flux not applied',
    output_policy='one final native snapshot~39MiB plus first-RHS diagnostic packets~70MiB; no old data overwritten',
    reference=reference['final'])
files=sorted(f for f in dest.glob('khout.*.dat') if f.name!='khout.000001.dat')
if code==0 and files:
    p,m=read_snapshot(files[-1]);q,mq=read_snapshot(reference['final'])
    a=np.argsort(p['u4if']);b=np.argsort(q['u4if'])
    differing=[k for k in p.dtype.names if not np.array_equal(p[k][a],q[k][b])]
    row.update(final=str(files[-1]),all_recorded_fields_bitwise_equal=not differing,
        differing_fields=differing,metadata_equal=m==mq,
        accepted_steps=(dest/'log').read_text().count('[SSP_ACCEPT]'),
        elapsed_seconds=json.loads((dest/'timing.json').read_text())['elapsed_seconds'])
(dest/'manifest.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row,indent=2),flush=True)
if code or not row.get('all_recorded_fields_bitwise_equal') or not row.get('metadata_equal'):
    raise SystemExit('Observer equivalence failed; no physics promotion')
