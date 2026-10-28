"""Stage complete moving-power-face kinematic mass-sweep audit."""
from pathlib import Path
import shutil
import subprocess
from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_weight_sweep_full_t6501_cpu_v1')
FILES=('kepler_fresh_weight_sweep_full_t6501_cpu_20261027.slurm',
       'audit_kepler_fresh_weight_sweep_full_speed_20261027.py',
       'audit_kepler_bounded_pressure_force_t10_20261027.py',
       'audit_kepler_sparse_area_t0825_20261027.py',
       'laguerre_newton_volume_t6_20261026.py',
       'laguerre_nonlinear_volume_probe_t6_20261026.py',
       'prepare_kepler_sparse_area_restart_20261027.py')

def main():
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(DOCS/FILES[0])],check=True)
    RUN.mkdir();hashes={}
    for name in FILES:
        source=DOCS/name;shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Changed {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)

if __name__=='__main__':main()
