"""Stage a one-step weighted conservation diagnostic from accepted t=0.825."""

import hashlib
import json
from pathlib import Path
import shutil


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT=ROOT/'gpu_energy_neutral_capacity64_v5_recover_to0825_a10_v1'
RUN=ROOT/'gpu_energy_neutral_capacity64_v5_fresh_localwork_t0826_ada_v1'
DOCS=Path(__file__).resolve().parent
EXPECTED_BINARY='feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'
EXPECTED_INITIAL='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    parent_receipt=PARENT/'audit_t0825.json'
    audit=json.loads(parent_receipt.read_text())
    if (audit['time']!=.825 or len(audit['checks'])!=21 or
            not all(audit['checks'].values()) or
            audit['initial_sha256']!=EXPECTED_INITIAL or
            audit['binary_sha256']!=EXPECTED_BINARY):
        raise ValueError('Fresh weighted parent is not accepted')
    endpoint=Path(audit['endpoint'])
    if (endpoint.parent!=PARENT or digest(endpoint)!=audit['endpoint_sha256'] or
            digest(Path(audit['initial']))!=EXPECTED_INITIAL or
            digest(PARENT/'eunha2.exe')!=EXPECTED_BINARY):
        raise ValueError('Weighted parent files changed')
    origins={
        'eunha2.exe':PARENT/'eunha2.exe',
        'params.dat':PARENT/'params.dat',
        'disk_active_mask.bin':PARENT/'disk_active_mask.bin',
        'initial.dat':Path(audit['initial']),
        'khout.000001.dat':endpoint,
        'job_ada.sh':DOCS/'kepler_q2_weighted_fresh_localwork_t0826_ada_20261027.sh',
    }
    for name in ('audit_kepler_q2_weighted_equilibrium_20261027.py',
                 'audit_kepler_bounded_native_evolution_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'audit_kepler_sparse_native_ballistic_20261027.py',
                 'audit_kepler_sparse_state_consistency_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        origins[name]=DOCS/name
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
    manifest={'parent':str(PARENT),'parent_receipt':str(parent_receipt),
              'parent_receipt_sha256':digest(parent_receipt),
              'parent_endpoint_sha256':audit['endpoint_sha256'],
              'initial_sha256':EXPECTED_INITIAL,'binary_sha256':EXPECTED_BINARY,
              'time_start':.825,'time_end':.826,
              'inputs':{name:{'source':str(origins[name]),'sha256':value}
                        for name,value in hashes.items()},
              'scope':'One-step fresh weighted local-work diagnostic on one Ada GPU; parent accepted 21/21. Separate from active long Slurm chain.'}
    with (RUN/'stage_manifest.json').open('x') as stream:
        json.dump(manifest,stream,indent=2)
        stream.write('\n')
    print(json.dumps({'run':str(RUN),'inputs':len(hashes),
                      'parent_sha256':audit['endpoint_sha256']}))


if __name__=='__main__':
    main()
