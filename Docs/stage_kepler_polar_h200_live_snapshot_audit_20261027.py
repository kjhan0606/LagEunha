"""Stage read-only audit for the completed H200 t~1.5 snapshot."""
from pathlib import Path
import shutil
import subprocess
from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/live_t150_audit_cpu_v1')
FILES=('kepler_polar_h200_live_snapshot_audit_20261027.slurm',
       'audit_kepler_polar_long_gate_20261027.py',
       'audit_kepler_native_radial_20261026.py')

def main():
    if RUN.exists(): raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(DOCS/FILES[0])],check=True)
    RUN.mkdir()
    hashes={}
    for name in FILES:
        source=DOCS/name
        shutil.copy2(source,RUN/name)
        hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)

if __name__=='__main__': main()
