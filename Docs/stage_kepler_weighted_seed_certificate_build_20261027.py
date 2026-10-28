"""Stage the opt-in seed-bin certificate on an isolated v5 source build."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS = Path(__file__).resolve().parent
V5 = ROOT/'gpu_energy_neutral_capacity64_build_v5'
RUN = ROOT/'gpu_energy_neutral_seed_certificate_capacity64_build_v9'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parent = json.loads((V5/'stage_manifest.json').read_text())
    template = (V5/'job.slurm').read_text()
    if digest(V5/'job.slurm') != parent['inputs']['job.slurm']['sha256']:
        raise ValueError('v5 job changed')
    replacements = {
        'LagQ2Cap64Build': 'LagQ2SeedCertBuild9',
        'lag_q2_capacity64_': 'lag_q2_seedcert9_capacity64_',
        'eunha2_capacity64_v5.exe': 'eunha2_seedcert_capacity64_v9.exe',
        'python3 patch_capacity.py "$build_dir" > patch_capacity_receipt.json':
            'python3 patch_capacity.py "$build_dir" > patch_capacity_receipt.json\n'
            'python3 patch_seedcert.py "$build_dir" > patch_seedcert_receipt.json\n'
            'g++ -std=c++14 -O2 -Wall -Wextra -Werror -I "$build_dir/Exam" '
            'test_seedcert.cpp -o "$build_dir/test_seedcert"\n'
            '"$build_dir/test_seedcert" > test_seedcert.log',
    }
    for old, new in replacements.items():
        count = 2 if old == 'eunha2_capacity64_v5.exe' else 1
        if template.count(old) != count:
            raise ValueError(f'Unexpected v5 job anchor: {old!r}')
        template = template.replace(old, new)
    names = ('source.tar', 'patch_guard.py', 'patch_iterations.py', 'patch_zero.py',
             'patch_capacity.py', 'gfs_pressure_neutral_projection_20261027.h',
             'test_projection.cpp')
    sources = {name: V5/name for name in names}
    sources.update({'patch_seedcert.py': DOCS/'apply_kepler_weighted_seed_certificate_20261027.py',
                    'gfs_seed_bin_certificate_20261027.h': DOCS/'gfs_seed_bin_certificate_20261027.h',
                    'test_seedcert.cpp': DOCS/'test_gfs_seed_bin_certificate_20261027.cpp'})
    for name in names:
        if digest(sources[name]) != parent['inputs'][name]['sha256']:
            raise ValueError(f'Changed v5 input: {name}')
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    (RUN/'job.slurm').write_text(template)
    subprocess.run(['bash', '-n', str(RUN/'job.slurm')], check=True)
    hashes = {'job.slurm': digest(RUN/'job.slurm')}
    for name, source in sources.items():
        if name == 'source.tar':
            (RUN/name).symlink_to(source.resolve())
        else:
            shutil.copy2(source, RUN/name)
        hashes[name] = digest(RUN/name)
        if hashes[name] != digest(source):
            raise ValueError(f'Input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent': str(V5), 'inputs': hashes,
                   'scope': 'Experimental opt-in seed-bin exclusion; native correctness, '
                            'fallback and performance validation required.'}, out, indent=2)
        out.write('\n')
    print(RUN)


if __name__ == '__main__':
    main()
