"""Stage a one-step local-weight-work probe from accepted fresh t=6.5."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT=ROOT/'gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1'
RUN=ROOT/'gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1'
DOCS=Path(__file__).resolve().parent
RECEIPT=DOCS/'KEPLER_FRESH_WEIGHTED_T65_2026-10-27/weighted_receipt.json'
BINARY_SHA='feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'
INITIAL_SHA='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    audit=json.loads(RECEIPT.read_text())
    endpoint=Path(audit['endpoint']);initial=Path(audit['initial'])
    if (audit['time']!=6.5 or len(audit['checks'])!=21 or
        not all(audit['checks'].values()) or audit['initial_sha256']!=INITIAL_SHA or
        audit['binary_sha256']!=BINARY_SHA or endpoint.parent!=PARENT or
        digest(endpoint)!=audit['endpoint_sha256'] or digest(initial)!=INITIAL_SHA or
        digest(PARENT/'eunha2.exe')!=BINARY_SHA):
        raise ValueError('Unaccepted or changed fresh t=6.5 parent')
    if RUN.exists():raise FileExistsError(RUN)
    template=DOCS/'kepler_q2_weighted_fresh_localwork_t0826_ada_20261027.sh'
    runner=template.read_text()
    for old,new in {'HYDRO_TSTOP=0.826':'HYDRO_TSTOP=6.501',
                    '--parent-time 0.825 --time 0.826':'--parent-time 6.5 --time 6.501',
                    '--output audit_t0826.json':'--output audit_t6501.json'}.items():
        if runner.count(old)!=1:raise ValueError(f'Unexpected runner anchor: {old}')
        runner=runner.replace(old,new)
    origins={'eunha2.exe':PARENT/'eunha2.exe','params.dat':PARENT/'params.dat',
             'disk_active_mask.bin':PARENT/'disk_active_mask.bin','initial.dat':initial,
             'khout.000001.dat':endpoint}
    for name in ('audit_kepler_q2_weighted_equilibrium_20261027.py',
                 'audit_kepler_bounded_native_evolution_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'audit_kepler_sparse_native_ballistic_20261027.py',
                 'audit_kepler_sparse_state_consistency_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        origins[name]=DOCS/name
    RUN.mkdir();hashes={}
    for name,source in origins.items():
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Changed while copying: {name}')
    (RUN/'job_ada.sh').write_text(runner);hashes['job_ada.sh']=digest(RUN/'job_ada.sh')
    with (RUN/'inputs.sha256').open('x') as stream:
        for name,value in sorted(hashes.items()):stream.write(f'{value}  {name}\n')
    manifest={'parent_receipt':str(RECEIPT),'parent_receipt_sha256':digest(RECEIPT),
              'parent_endpoint_sha256':audit['endpoint_sha256'],'initial_sha256':INITIAL_SHA,
              'binary_sha256':BINARY_SHA,'time_start':6.5,'time_end':6.501,
              'runner_template_sha256':digest(template),'inputs':hashes,
              'scope':'Isolated one-step Ada replay from accepted fresh-weighted t=6.5.'}
    with (RUN/'stage_manifest.json').open('x') as stream:
        json.dump(manifest,stream,indent=2);stream.write('\n')
    print(RUN)

if __name__=='__main__':main()
