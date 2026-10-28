"""Stage counterfactual spin-energy positivity analysis on saved states."""

from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT=ROOT/'affine_spin_energy_t1651_v1'
RUN=ROOT/'affine_spin_thermal_risk_t1651_v1'


def main():
    source_job=(PARENT/'job.slurm').read_text()
    substitutions={'LagQ2AffineSpin':'LagQ2SpinThermalRisk',
                   'w0_receipt.json --output w0.json':'w0_receipt.json --thermal-risk --output w0.json',
                   'weighted_receipt.json --output weighted.json':'weighted_receipt.json --thermal-risk --output weighted.json'}
    for old,new in substitutions.items():
        if source_job.count(old)!=1:raise ValueError(f'Unexpected job anchor: {old!r}')
        source_job=source_job.replace(old,new)
    if RUN.exists():raise FileExistsError(RUN)
    RUN.mkdir();(RUN/'job.slurm').write_text(source_job)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    hashes={'job.slurm':digest(RUN/'job.slurm')}
    for line in (PARENT/'inputs.sha256').read_text().splitlines():
        expected,name=line.split('  ',1)
        if name=='job.slurm':continue
        if digest(PARENT/name)!=expected:raise ValueError(f'Changed parent input: {name}')
        source=DOCS/name if name=='audit_kepler_affine_spin_energy_20261027.py' else PARENT/name
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
