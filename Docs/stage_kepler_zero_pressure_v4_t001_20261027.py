"""Stage a v4 fresh native step against the completed v3 weighted pilot."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS=Path(__file__).resolve().parent
PILOT=ROOT/'gpu_energy_neutral_weighted_fresh_jitter64_t001_a10_v1'
BUILD=ROOT/'gpu_energy_neutral_zero_pressure_build_v4'
RUN=ROOT/'gpu_energy_neutral_zero_pressure_v4_t001_a10_v1'
V3='7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38'
V4='3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc'
INITIAL='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    if RUN.exists():raise FileExistsError(RUN)
    receipt=json.load(open(PILOT/'audit_t001.json'))
    if (not all(receipt['checks'].values()) or receipt['binary_sha256']!=V3
            or receipt['reference_sha256']!=INITIAL
            or digest(Path(receipt['endpoint']))!=receipt['endpoint_sha256']):
        raise ValueError('Reference v3 pilot is not accepted or changed')
    if digest(BUILD/'eunha2_zero_pressure_v4.exe')!=V4:
        raise ValueError('New v4 binary changed')
    names=('glass.bin','areas.bin','params.dat','accepted_initial.dat','ic_manifest.json',
           'disk_active_mask.bin','audit_kepler_q2_weighted_fresh_pilot_20261027.py',
           'audit_kepler_sparse_area_t0825_20261027.py',
           'audit_kepler_sparse_state_consistency_20261027.py',
           'laguerre_nonlinear_volume_probe_t6_20261026.py',
           'prepare_kepler_sparse_area_restart_20261027.py')
    sources={name:PILOT/name for name in names}
    sources.update({'eunha2.exe':BUILD/'eunha2_zero_pressure_v4.exe',
                    'job.slurm':DOCS/'kepler_q2_zero_pressure_v4_t001_20261027.slurm',
                    'audit_kepler_zero_pressure_v4_pair_20261027.py':DOCS/'audit_kepler_zero_pressure_v4_pair_20261027.py'})
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    hashes={}
    for name,source in sources.items():
        shutil.copy2(source,RUN/name)
        hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as f:
        for name,value in sorted(hashes.items()):f.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump({'run':str(RUN),'reference_receipt':str(PILOT/'audit_t001.json'),
                   'reference_receipt_sha256':digest(PILOT/'audit_t001.json'),
                   'binary_sha256':V4,'inputs':{name:{'source':str(sources[name]),'sha256':hashes[name]}
                                                 for name in sorted(hashes)},
                   'scope':'One native v4 nonzero-gradient regression step; no long-time promotion.'},f,indent=2)
        f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(hashes),'binary_sha256':V4},indent=2))


if __name__=='__main__':
    main()
