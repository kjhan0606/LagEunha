"""Exercise ABBA auditor logic on explicitly mocked repeated real short logs."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

docs=Path(__file__).resolve().parent
base=Path('/gpfs/kjhan/LagForce/Kepler_WallGuardOMP_20261029')
audit=docs/'audit_wall_guard_omp_abba_20261029.py';initial=(base/'omp/khout.000022.dat').resolve()
with tempfile.TemporaryDirectory(prefix='mock_wall_guard_abba_') as tmp:
    root=Path(tmp)
    stage={'initial':str(initial),'initial_sha256':'0f60d80dec2668a7e53a6e99ec84760af824ed70fb7a868ddd605fe53ffe11bc',
           'binary_sha256':'da70844984a52fec77e2d6dec451e9b5695e761ba93a42d0f0d079c1091edf44',
           'target_time':0.0018340842946822263,'order':['A1','B1','B2','A2'],'minimum_speedup':0.01}
    for name in stage['order']:
        src=base/('omp' if name.startswith('B') else 'serial');d=root/name;d.mkdir()
        for p in ('khout.000026.dat','eunha2.exe','run.sh','exit_status.txt'):(d/p).symlink_to(src/p)
        # The real four-step logs are deliberately repeated to exercise the
        # minimum-sample path. This creates no performance or physics result.
        (d/'log').write_text((src/'log').read_text()*3)
        (d/'elapsed_seconds.txt').write_text('1.0\n')
    (root/'stage.json').write_text(json.dumps(stage))
    run=subprocess.run([sys.executable,str(audit),str(root)],capture_output=True,text=True)
    if run.returncode:raise RuntimeError(run.stdout[-1500:]+run.stderr[-1500:])
    print('PASS: mocked positive ABBA auditor control')
    (root/'abba_audit.json').unlink();stage['minimum_speedup']=1e9
    (root/'stage.json').write_text(json.dumps(stage))
    run=subprocess.run([sys.executable,str(audit),str(root)],capture_output=True,text=True)
    out=json.loads((root/'abba_audit.json').read_text())
    if run.returncode==0 or out['checks']['speedup_at_least_1p05']:
        raise RuntimeError('unmet speed gate was accepted')
    print('PASS: mocked insufficient-speed rejection')
