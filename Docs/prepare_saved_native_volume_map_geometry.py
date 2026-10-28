"""Preserved native .015/actual .0175/.02 geometry only; exactly three calls.

No RHS, density fit, root solve, new accepted state or restart.
"""
import argparse,concurrent.futures,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    base=Path('/gpfs/kjhan/LagForce');native=base/'Codex_full_disk_native_samehost_original_and_coarse_20261019/mpi_coarse'
    build=base/'Codex_full_disk_native_evolution_build_v8_20261019'
    sys.path.insert(0,str(build));from full_moment_state_io import load,HEADER
    folders=dict(base=native/'accepted_0003',half=base/'Codex_full_disk_native_preserved_midpoint_e8_20261019',next=native/'accepted_0004')
    accepted=json.loads((native/'accepted_0004.json').read_text());h=accepted['h'];assert h==.005 and accepted['accepted_commit']
    commands=json.loads((build/'commands.json').read_text());assert '-ffp-contract=off' in commands['full_moment_update']
    assert not any('fast-math' in flag and flag!='-fno-fast-math' for flag in commands['full_moment_update'])
    data={name:load(path) for name,path in folders.items()}
    headers={name:HEADER.unpack((path/'header.bin').read_bytes()) for name,path in folders.items()}
    for name,epoch,t in [('base',7,3*h),('half',8,3.5*h),('next',9,4*h)]:
        assert headers[name][3]==epoch and float(data[name]['time'])==t
        assert headers[name][4]==headers['base'][4]==30070787499836263
        assert np.array_equal(data[name]['ids'],data['base']['ids']) and np.array_equal(data[name]['mass'],data['base']['mass'])
    n=len(data['base']['ids']);a=data['base']['state'][:-4].reshape(n,20);b=data['next']['state'][:-4].reshape(n,20)
    expected=a[:,2]+h*data['next']['weight_rate_seed']
    assert np.array_equal(expected.view(np.uint64),b[:,2].view(np.uint64)), 'Actual z_half seed/update lineage mismatch'
    preservation=json.loads((folders['half']/'preservation.json').read_text())
    assert all(sha(folders['half']/file)==digest for file,digest in preservation['sha256'].items())
    files=('header','ids','mass','state','factors','Psi','weight_rate_seed')
    hashes={str(path/(f+'.bin')):sha(path/(f+'.bin')) for path in folders.values() for f in files}
    estimate=2*1024**3;assert shutil.disk_usage(args.output).free>estimate
    manifest=dict(status='validated_saved_native_triplet_before_geometry',h=h,step=4,lineage=headers,source_folders={k:str(v) for k,v in folders.items()},
        source_sha256=hashes,update_contract='binary64 separate multiply/add; fp-contract=off',
        update_source_sha256=sha(build/'full_moment_update.c'),build_commands_sha256=sha(build/'commands.json'),
        expected_binary_sha256=(build/'binary.sha256').read_text().strip(),
        workers_per_geometry=4,concurrent_geometries=3,total_worker_cpus=12,estimated_output_upper_bytes=estimate,
        physical_providers=0,accepted_time_steps=0,density_fits=0,root_solves=0)
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    script=base/'Codex_full_disk_arc_roots_modules_v7_20261019/prepare_full_curved_stage_geometry_workers.py'
    for name,state in data.items():
        folder=args.output/name;folder.mkdir();np.savez(folder/'state.npz',**state)
    def geometry(name):
        folder=args.output/name
        with (folder/'geometry.log').open('x') as log:
            subprocess.run([sys.executable,str(script),'--state',str(folder/'state.npz'),'--workers','4','--output',str(folder/'geometry')],stdout=log,stderr=subprocess.STDOUT,check=True)
        return name
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        completed=list(pool.map(geometry,folders))
    assert all(sha(file)==digest for file,digest in hashes.items())
    manifest.update(status='completed_three_saved_native_geometries',completed=completed,source_geometry_script_sha256=sha(script),original_sources_unchanged=True)
    (args.output/'report.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(status=manifest['status'],completed=completed,physical_providers=0,accepted_time_steps=0)),flush=True)


if __name__=='__main__':main()
