"""One effective time-step control, reusing the graded positive preflight."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import numpy as np
from audit_snapshot import read_snapshot
from run_radial_resolution_v1 import audit_state, EXPECTED_BINARY

root = Path(sys.argv[1]).resolve()
base = root/'radial_resolution_20261021_v2'
nominal = base/'native_v2/graded/t0.05'
initial = base/'graded/khout.000001.dat'
source = base/'native_v2/graded/preflight_mpi1/khout.000002.dat'
dest = base/'half_step_t005_v1'
dest.mkdir(exist_ok=False)
read_dt = lambda text: [float(x) for x in re.findall(r'^Time=.*?\bdt=\s*([^\s]+)', text, re.M)]
nominal_dt = read_dt((nominal/'log').read_text())
if len(nominal_dt) < 10 or read_snapshot(source)[1]['t'] != 1e-7:
    raise ValueError('nominal accepted time series or common positive start missing')
ceiling = .5*float(np.median(nominal_dt))
binary = root/'raw_vertex_dedup/v4/eunha2_raw_vertex_dedup.exe'
if hashlib.sha256(binary.read_bytes()).hexdigest() != EXPECTED_BINARY:
    raise ValueError('frozen native binary changed')
(dest/'eunha2.exe').symlink_to(binary)
(dest/'khout.000002.dat').symlink_to(source)
(dest/'params.dat').symlink_to(base/'graded/params.dat')
env = os.environ.copy()
for key in list(env):
    if key.startswith(('GFS_', 'EUNHA_', 'HYDRO_', 'KEPLER_', 'SEDOV_', 'LAGMFM_')):
        env.pop(key)
env.update(json.loads((nominal/'environment.json').read_text()))
env['GFS_DIAG_DT_MAX'] = repr(ceiling)
for key in ('I_MPI_PMI_LIBRARY', 'PMI_RANK', 'PMI_SIZE', 'PMI_FD', 'PMI_PORT',
            'OMP_PLACES', 'GOMP_CPU_AFFINITY', 'KMP_AFFINITY'):
    env.pop(key, None)
(dest/'environment.json').write_text(json.dumps({k:v for k,v in env.items()
    if k.startswith(('GFS_', 'EUNHA_', 'HYDRO_', 'KEPLER_', 'SEDOV_', 'OMP_', 'I_MPI_', 'NATIVE_'))}, indent=2)+'\n')
row = dict(scope='ONE graded half-step control to matched t=.05; no full trajectory restart or dt scan',
           binary_sha256=EXPECTED_BINARY, positive_start=str(source), requested_stop=.05,
           nominal_dt_median=float(np.median(nominal_dt)), applied_dt_ceiling=ceiling,
           nominal_accepted_steps=len(nominal_dt), nominal_RHS=4*len(nominal_dt),
           long_goal_achieved=False)
(dest/'manifest.json').write_text(json.dumps(row, indent=2)+'\n')
code = subprocess.run([sys.executable, str(Path(__file__).with_name('run_native.py'))], cwd=dest, env=env).returncode
row['native_exit_code'] = code
dt = read_dt((dest/'log').read_text())
row['accepted_steps'] = len(dt)
row['completed_step_RHS'] = 4*len(dt)
row['actual_dt_max'] = max(dt) if dt else None
row['effective_half_step_control'] = bool(dt and max(dt) <= ceiling*(1+1e-5)
                                         and len(dt) >= 1.5*len(nominal_dt))
gamma = json.loads((base/'manifest.json').read_text())['gamma']
if code == 0:
    snapshots = [p for p in dest.glob('khout.*.dat') if not p.is_symlink()]
    final = max(snapshots, key=lambda p: int(p.name.split('.')[1]))
    if read_snapshot(final)[1]['t'] != .05:
        raise ValueError('half-step endpoint mismatch')
    row['physical_audit'] = audit_state(final, initial, gamma)
    native = json.loads((base/'native_v2/graded/manifest.json').read_text())
    row['nominal_physical_audit'] = next(s['physical_audit'] for s in native['stages'] if s['label']=='t0.05')
    row['status'] = 'MATCHED_HALF_STEP_COMPLETE' if row['effective_half_step_control'] else 'INEFFECTIVE_TIME_CONTROL'
else:
    starts = list(dest.glob('khdiag.step*.start.dat'))
    if starts:
        row['normal_start_physical_audit'] = audit_state(starts[-1], initial, gamma)
    row['status'] = 'NATIVE_FAILURE; INVALID_STAGE_PRESERVED'
(dest/'manifest.json').write_text(json.dumps(row, indent=2)+'\n')
print(json.dumps(row, indent=2))
raise SystemExit(code)
