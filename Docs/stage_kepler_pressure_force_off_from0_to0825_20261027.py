"""Stage an isolated short t=0 pressure-force-off causal control."""

import hashlib
import json
from pathlib import Path
import shutil


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS=Path(__file__).resolve().parent
SOURCE=ROOT/'jitter005_t001_a10_v1'
DIAGNOSTIC=ROOT/'pressure_force_off_jitter_to10_a10_v1'
RUN=ROOT/'pressure_force_off_jitter_from0_to0825_a10_v1'
PARENT_AUDIT=ROOT/'jitter005_to65_a10_v1/audit_t65.json'
INITIAL_SHA='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'
BINARY_SHA='0a65c373613a5b31f8449b456387ebcf80d72136a81e46eec10c9692a8bed5e7'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    parent=json.loads(PARENT_AUDIT.read_text())
    if parent['time']!=6.5 or not all(parent['checks'].values()):
        raise ValueError('Parent native physics audit did not pass')
    if parent['initial_sha256']!=INITIAL_SHA or digest(Path(parent['initial']))!=INITIAL_SHA:
        raise ValueError('Accepted t=0 snapshot changed')
    origins={
        'eunha2.exe':DIAGNOSTIC/'eunha2.exe',
        'params.dat':SOURCE/'params.dat',
        'glass.bin':SOURCE/'glass.bin',
        'areas.bin':SOURCE/'areas.bin',
        'khout.000001.dat':Path(parent['initial']),
        'reference_native_t0825.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T0825_2026-10-27.json',
        'job.slurm':DOCS/'kepler_pressure_force_off_from0_to0825_a10_20261027.slurm',
    }
    for name in ('audit_kepler_pressure_force_off_initial_t0825_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'audit_kepler_sparse_native_ballistic_20261027.py',
                 'audit_kepler_sparse_state_consistency_20261027.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        origins[name]=DOCS/name
    if digest(origins['eunha2.exe'])!=BINARY_SHA:
        raise ValueError('Pressure-force-off binary changed')
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    hashes={}
    for name,source in origins.items():
        shutil.copy2(source,RUN/name)
        hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):
            raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as stream:
        for name,value in sorted(hashes.items()):
            stream.write(f'{value}  {name}\n')
    manifest={'run':str(RUN),'accepted_parent_audit':str(PARENT_AUDIT),
              'accepted_parent_audit_sha256':digest(PARENT_AUDIT),
              'initial_sha256':INITIAL_SHA,'diagnostic_binary_sha256':BINARY_SHA,
              'target_time':0.825,
              'inputs':{name:{'source':str(origins[name]),'sha256':value}
                        for name,value in hashes.items()},
              'scope':'Short t0→0.825 causal diagnostic only. Pressure force/work adjoint broken; not an accepted hydrodynamic solution.'}
    with (RUN/'stage_manifest.json').open('x') as stream:
        json.dump(manifest,stream,indent=2)
        stream.write('\n')
    print(json.dumps({'run':str(RUN),'inputs':len(hashes),'initial_sha256':INITIAL_SHA,
                      'binary_sha256':BINARY_SHA}))


if __name__=='__main__':
    main()
