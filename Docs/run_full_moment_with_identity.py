"""Launch unchanged native MPI executable, capture identity from its log events.

No polling loop: process headers and GPU selection messages trigger one read
of that exact PID's mappings. No state update, timestep control or retry.
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
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--provider',type=Path,required=True)
    p.add_argument('--params',type=Path,required=True)
    p.add_argument('--ranks',type=int,choices=(1,2),required=True)
    p.add_argument('--log',type=Path,required=True)
    p.add_argument('--identity',type=Path,required=True)
    args=p.parse_args()
    expected=sha(args.exe);operator_sha=os.environ.get('EUNHA_MOMENT_OPERATOR_SHA','51014e2fb7de34bb8eb1ed99c8ed0e22ecbbe5690a67327f8fec9892d5968471')
    assert expected==os.environ['EUNHA_MOMENT_BINARY_SHA']
    provenance=json.loads((args.provider.parent/'provenance.json').read_text())
    assert provenance['binary_sha256']==expected
    sources={name:sha(args.provider.parent/name) for name in provenance['sources']}
    assert sources==provenance['sources'],'Published provider/native snapshot modified'
    record=dict(host=socket.gethostname(),slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        binary_sha256=expected,params_sha256=sha(args.params),
        provider_script=str(args.provider.resolve()),provider_snapshot_sources_sha256=sources,
        ranks={},GPU_UUID_by_rank={},identity_errors=[],physical_providers=0)
    command=['mpirun','-np',str(args.ranks),'-ppn',str(args.ranks),str(args.exe),str(args.params)]
    with args.log.open('x') as log:
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        for line in process.stdout:
            log.write(line);log.flush()
            if '[FULL_MOMENT_NATIVE_RUNTIME]' in line:
                try:
                    rank=int(re.search(r'rank=(\d+)',line).group(1))
                    pid=int(re.search(r'pid=(\d+)',line).group(1))
                    assert re.search(r'binary_sha=([0-9a-f]{64})',line).group(1)==expected
                    root=Path(f'/proc/{pid}')
                    executable=(root/'exe').resolve(strict=True)
                    assert executable==args.exe.resolve() and sha(executable)==expected
                    mapped=set()
                    for row in (root/'maps').read_text().splitlines():
                        parts=row.split(maxsplit=5)
                        if len(parts)==6 and parts[5].endswith('/liblagmoment_owned.so'):mapped.add(parts[5])
                    assert len(mapped)==1
                    libraries={path:sha(path) for path in mapped}
                    assert list(libraries.values())==[operator_sha]
                    children=list(map(int,(root/'task'/str(pid)/'children').read_text().split()))
                    provider_processes={}
                    for child in children:
                        cmdline=Path(f'/proc/{child}/cmdline').read_bytes().split(b'\0')
                        if str(args.provider.resolve()).encode() in cmdline:
                            provider_processes[str(child)]=[part.decode() for part in cmdline if part]
                    assert len(provider_processes)==1
                    record['ranks'][str(rank)]=dict(pid=pid,actual_executable=str(executable),
                        actual_mapped_operator_libraries=libraries,actual_provider_process=provider_processes)
                except Exception as exc:
                    record['identity_errors'].append(repr(exc))
            if '[FULL_MOMENT_ACTUAL_GPU]' in line:
                rank=re.search(r'rank=(\d+)',line).group(1)
                uuid=re.search(r'UUID=([0-9a-f]{32})',line).group(1)
                previous=record['GPU_UUID_by_rank'].setdefault(rank,uuid)
                if previous!=uuid:record['identity_errors'].append('Rank changed GPU UUID: '+rank)
        status=process.wait()
    if status==0:
        if set(record['ranks'])!={str(rank) for rank in range(args.ranks)}:
            record['identity_errors'].append('Missing native process identity')
        if len(record['GPU_UUID_by_rank'])!=args.ranks or len(set(record['GPU_UUID_by_rank'].values()))!=args.ranks:
            record['identity_errors'].append('Missing or duplicate native GPU UUID')
    record.update(native_exit_code=status,status='passed_live_native_runtime_identity'
        if status==0 and not record['identity_errors'] else 'failed_native_run_or_identity_capture')
    with args.identity.open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(dict(status=record['status'],native_exit_code=status,identity_errors=record['identity_errors'])),flush=True)
    return status if status else (1 if record['identity_errors'] else 0)


if __name__=='__main__':raise SystemExit(main())
