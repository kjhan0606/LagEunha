"""Stage independent H200 Q2 pressure comparison from checked t30 restarts."""
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

docs=Path(__file__).resolve().parent
root=Path('/gpfs/kjhan/LagForce/Kepler_Q2_PressureCompare_20261029')
pre=Path('/tmp/kepler_pressure_scale_preflight_20261029')
parent=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to30_a10_v1/khout.005247.dat')
binary=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/eunha2.exe')
params=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/params.dat')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if root.exists():raise FileExistsError(root)
audit=json.loads((pre/'audit.json').read_text())
if sha(parent)!=audit['source_sha256'] or sha(pre/'x1.dat')!=audit['x1_sha256'] or sha(pre/'x4b.dat')!=audit['x4_sha256']:
    raise ValueError('changed verified input')
if sha(binary)!='87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968':raise ValueError('changed native binary')
if audit['time']!=30 or audit['cells']!=55351 or audit['eos_relative_error_max']>2e-7:raise ValueError('preflight failed')
root.mkdir()
inputs=[]
for k,src in ((1,pre/'x1.dat'),(4,pre/'x4b.dat')):
    d=root/f'x{k}';d.mkdir()
    shutil.copyfile(src,d/'khout.000001.dat');shutil.copyfile(params,d/'params.dat')
    os.symlink(binary,d/'eunha2.exe')
    inputs.extend((d/'khout.000001.dat',d/'params.dat',d/'eunha2.exe'))
for name in ('kepler_q2_pressure_compare_20261029.slurm','audit_kepler_pressure_scaled_ic_20261029.py','kepler_scale_thermal_checkpoint_20261029.c'):
    shutil.copyfile(docs/name,root/name);inputs.append(root/name)
shutil.copyfile(pre/'audit.json',root/'ic_audit.json');inputs.append(root/'ic_audit.json')
contract={'source':str(parent),'source_sha256':sha(parent),'binary_sha256':sha(binary),
          'ic_audit_sha256':sha(root/'ic_audit.json'),'start_time':30.0,'target_time':32.501708696164926,
          'runs':[{'array_task':1,'pressure_scale':1},{'array_task':4,'pressure_scale':4}],
          'resources':'Each independent array element: 1 H200, 2 CPU, 32GB, 2h; no Slurm dependencies.',
          'outputs':'Each independent case writes one final ~32MB snapshot and one map plus logs; no periodic snapshots.',
          'scope':'Counterfactual short Q2 pressure sensitivity from a preserved t30 checkpoint, not evidence of t600 stability or the physical solution of the x4 state.'}
(root/'contract.json').write_text(json.dumps(contract,indent=2)+'\n');inputs.append(root/'contract.json')
(root/'inputs.sha256').write_text(''.join(sha(p)+'  '+str(p)+'\n' for p in inputs))
subprocess.run(['bash','-n',str(root/'kepler_q2_pressure_compare_20261029.slurm')],check=True)
print(json.dumps(contract,indent=2))
