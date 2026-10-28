"""Read-only identity capture on the existing allocation; no RHS or test.

Inspect only the explicitly identified process(es), their direct children,
and named immutable source files. Never scan process trees or directory trees.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import subprocess


def sha(file):
    h=hashlib.sha256()
    with Path(file).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--oracle',type=Path)
    p.add_argument('--native-log',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    base=Path('/gpfs/kjhan/LagForce');pids=[]
    if args.oracle:
        preflight=json.loads((args.oracle/'preflight.json').read_text())
        pids.append(int(re.fullmatch(r'eunha_original_oracle_(\d+)_\d+',Path(preflight['memory']).name).group(1)))
        assert preflight['host'].split('.')[0]==socket.gethostname().split('.')[0]
    if args.native_log:
        for line in args.native_log.read_text().splitlines():
            if '[FULL_MOMENT_NATIVE_RUNTIME]' in line:
                pids.append(int(re.search(r'pid=(\d+)',line).group(1)))
    assert pids
    processes={}
    # One direct-child level covers each persistent provider or current clean
    # original RHS process; it does not walk multiprocessing worker descendants.
    for pid in list(pids):
        children=Path(f'/proc/{pid}/task/{pid}/children')
        if children.exists():pids.extend(map(int,children.read_text().split()))
    for pid in sorted(set(pids)):
        root=Path(f'/proc/{pid}')
        if not root.exists():continue  # A clean RHS child may have exited.
        try:
            executable=str((root/'exe').resolve(strict=True))
            cmdline=(root/'cmdline').read_bytes().replace(b'\0',b' ').decode()
            mapped=set()
            for line in (root/'maps').read_text().splitlines():
                fields=line.split(maxsplit=5)
                if len(fields)==6 and fields[5].endswith('/liblagmoment_owned.so'):mapped.add(fields[5])
            processes[str(pid)]=dict(executable=executable,executable_sha256=sha(executable),
                cmdline=cmdline,actual_mapped_operator_libraries={path:sha(path) for path in mapped})
        except FileNotFoundError:continue
    assert processes and any(row['actual_mapped_operator_libraries'] for row in processes.values())
    devices=subprocess.check_output(['nvidia-smi','--query-gpu=name,uuid','--format=csv,noheader'],text=True).splitlines()
    gpu_processes=[]
    for line in subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid','--format=csv,noheader'],text=True).splitlines():
        fields=[s.strip() for s in line.split(',')]
        if fields[0] in processes:gpu_processes.append(dict(pid=int(fields[0]),GPU_UUID=fields[1]))
    fast=base/'Codex_full_disk_fast_modules_v1_20261019'
    modules=[fast/(name+'.py') for name in ('full_disk_density_memory','full_disk_velocity_memory','full_disk_stage_memory','full_disk_factor_memory')]
    if args.oracle:modules.append(args.oracle.parent/'check_full_disk_samehost_original_midpoint.py')
    snapshots={}
    for version in ('v7','v8','v10'):
        folder=base/f'Codex_full_disk_native_evolution_build_{version}_20261019'
        snapshots[version]=dict(binary_sha256=sha(folder/'eunha2.exe'),
            provider_sha256=sha(folder/'full_moment_provider_service.py'),
            snapshot_manifest_sha256=sha(folder/'provenance.json'))
    ic=base/'Codex_full_curved_kepler_IC_20261019/initial.npz'
    result=dict(status='captured_live_full_moment_runtime_identity',host=socket.gethostname(),
        slurm_job_id=os.environ.get('SLURM_JOB_ID'),CUDA_VISIBLE_DEVICES=os.environ.get('CUDA_VISIBLE_DEVICES'),
        allocation_node_gpu_inventory=devices,own_live_gpu_contexts=gpu_processes,processes=processes,
        initial_sha256=sha(ic),initial_seed_source='Immutable initial.npz field weight_rate_seed; same IC',
        original_provider_module_sha256={str(path):sha(path) for path in modules},
        native_snapshots=snapshots,physical_providers=0,accepted_time_steps=0)
    with args.output.open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(dict(status=result['status'],host=result['host'],processes=list(processes),own_gpu_contexts=gpu_processes)),flush=True)


if __name__=='__main__':main()
