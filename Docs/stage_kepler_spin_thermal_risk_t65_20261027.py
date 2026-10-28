"""Superset of the queued t=6.5 spin diagnostic, including thermal risk."""

from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT=ROOT/'affine_spin_energy_t65_v1'
RUN=ROOT/'affine_spin_thermal_risk_t65_v1'


def main():
    hashes_parent={}
    for line in (PARENT/'inputs.sha256').read_text().splitlines():
        value,name=line.split('  ',1);hashes_parent[name]=value
    template=(PARENT/'job.slurm').read_text()
    if digest(PARENT/'job.slurm')!=hashes_parent['job.slurm']:
        raise ValueError('Prior t=6.5 runner changed')
    for old,new in {'LagQ2AffineSpin65':'LagQ2SpinThermal65',
                    'w0_receipt.json --output w0.json':'w0_receipt.json --thermal-risk --output w0.json',
                    'weighted_receipt.json --output weighted.json':'weighted_receipt.json --thermal-risk --output weighted.json'}.items():
        if template.count(old)!=1:raise ValueError(f'Unexpected runner anchor: {old!r}')
        template=template.replace(old,new)
    if RUN.exists():raise FileExistsError(RUN)
    RUN.mkdir();(RUN/'job.slurm').write_text(template)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    hashes={'job.slurm':digest(RUN/'job.slurm')}
    for name,expected in hashes_parent.items():
        if name=='job.slurm':continue
        if digest(PARENT/name)!=expected:raise ValueError(f'Prior input changed: {name}')
        source=DOCS/name if name=='audit_kepler_affine_spin_energy_20261027.py' else PARENT/name
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
