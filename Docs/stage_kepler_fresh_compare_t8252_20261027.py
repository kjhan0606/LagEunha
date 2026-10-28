"""Stage same-initial-state comparison at t≈8.252."""
import json
from pathlib import Path
import shutil
import subprocess
from prepare_kepler_sparse_area_restart_20261027 import digest
import stage_kepler_weighted_fresh_snapshot_compare_20261027 as base

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_compare_step904_t8252_cpu_v1')

def main():
    if RUN.exists():raise FileExistsError(RUN)
    template=(DOCS/base.FILES[0]).read_text()
    replacements={'LagQ2SnapshotCompare':'LagQ2FreshCompare8252',
                  'jitter005_to65_a10_v1':'jitter005_to10_a10_v1',
                  'gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1':'gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1',
                  'step=002553':'step=000904','--parent-time .825':'--parent-time 6.5'}
    for old,new in replacements.items():
        if template.count(old)!=1:raise ValueError(f'Unexpected runner anchor: {old}')
        template=template.replace(old,new)
    RUN.mkdir();(RUN/'job.slurm').write_text(template)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    hashes={'job.slurm':digest(RUN/'job.slurm')}
    for name in base.FILES[1:]:
        source=DOCS/name;shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Changed {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'template_sha256':digest(DOCS/base.FILES[0]),'input_sha256':hashes,
                   'scope':'Same-IC, independently evolved W0 and weighted t~8.252 snapshots.'},out,indent=2)
        out.write('\n')
    print(RUN)

if __name__=='__main__':main()
