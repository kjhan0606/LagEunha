"""Compare wall-select/cache weighted candidate to v5 from matched states."""

import json
from pathlib import Path
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import DOCS, ROOT, digest

PARENT=ROOT/'gpu_energy_neutral_seedcert_w0_weighted_speed_v1'
BUILD=ROOT/'gpu_energy_neutral_wallselect_cache_capacity64_build_v10'
RUN=ROOT/'gpu_energy_neutral_wallselect_cache_pair_v1'


def main():
    old_manifest=json.loads((PARENT/'stage_manifest.json').read_text())
    old_candidate=old_manifest['parent_manifest']['candidate_binary_sha256']
    candidate=BUILD/'eunha2_wallselect_cache_capacity64_v10.exe'
    candidate_sha=(BUILD/'binary.sha256').read_text().split()[0]
    if digest(candidate)!=candidate_sha:raise ValueError('Candidate binary changed')
    template=(PARENT/'job.slurm').read_text()
    replacements={
        'LagQ2SeedCertSpeed':'LagQ2WallCachePair',
        old_candidate:candidate_sha,
        'export GFS_GPU_SEED_BIN_CERT=1':': # wall-select candidate has no seed-certificate option',
    }
    for old,new in replacements.items():
        if template.count(old)!=1:raise ValueError(f'Unexpected runner anchor: {old!r}')
        template=template.replace(old,new)
    if RUN.exists():raise FileExistsError(RUN)
    RUN.mkdir();(RUN/'job.slurm').write_text(template)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    root_names=('audit_kepler_seed_certificate_speed_20261027.py',
                'audit_kepler_sparse_area_t0825_20261027.py',
                'prepare_kepler_sparse_area_restart_20261027.py')
    root_hashes={'job.slurm':digest(RUN/'job.slurm')}
    for name in root_names:
        shutil.copy2(DOCS/name,RUN/name);root_hashes[name]=digest(RUN/name)
        if root_hashes[name]!=digest(DOCS/name):raise ValueError(f'Changed root input: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(root_hashes.items()):out.write(f'{value}  {name}\n')
    for case in ('w0','weighted'):
        (RUN/case).mkdir()
        for mode in ('baseline','enabled'):
            source=PARENT/case/mode;target=RUN/case/mode;target.mkdir();hashes={}
            for line in (source/'inputs.sha256').read_text().splitlines():
                expected,name=line.split('  ',1)
                if digest(source/name)!=expected:raise ValueError(f'Changed {case}/{mode}/{name}')
                actual_source=candidate if name=='eunha2.exe' and mode=='enabled' else source/name
                shutil.copy2(actual_source,target/name);hashes[name]=digest(target/name)
                if hashes[name]!=digest(actual_source):raise ValueError(f'Copy changed: {name}')
            if mode=='enabled' and hashes['eunha2.exe']!=candidate_sha:raise ValueError('Candidate mismatch')
            with (target/'inputs.sha256').open('x') as out:
                for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent':str(PARENT),'candidate_sha256':candidate_sha,
                   'scope':'v5 vs wall-select plus cache, matching W0 and weighted states; '
                           'a successful short pair does not establish long-time equivalence.'},
                  out,indent=2);out.write('\n')
    print(RUN)


if __name__=='__main__':main()
