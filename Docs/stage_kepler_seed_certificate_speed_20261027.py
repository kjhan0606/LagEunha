"""Stage same-GPU short timing with certificate diagnostic atomics disabled."""

from pathlib import Path
import json
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import DOCS, ROOT, digest

PARENT=ROOT/'gpu_energy_neutral_seedcert_w0_weighted_pair_v1'
RUN=ROOT/'gpu_energy_neutral_seedcert_w0_weighted_speed_v1'


def main():
    gate=json.loads((PARENT/'pair_audit.json').read_text())
    if not all(gate['checks'].values()):raise ValueError('Device equivalence gate failed')
    manifest=json.loads((PARENT/'stage_manifest.json').read_text())
    template=(PARENT/'job.slurm').read_text()
    replacements={
        'LagQ2SeedCertPair':'LagQ2SeedCertSpeed',
        'for mode in baseline enabled fallback; do':'for mode in baseline enabled; do',
        'export GFS_GPU_SEED_BIN_CERT=1 GFS_GPU_SEED_BIN_CERT_DIAG=1':'export GFS_GPU_SEED_BIN_CERT=1',
        'HYDRO_TSTOP=6.501 EUNHA_DUMP_DT=6.501':'HYDRO_TSTOP=6.525 EUNHA_DUMP_DT=6.525',
        'HYDRO_TSTOP=.826 EUNHA_DUMP_DT=.826':'HYDRO_TSTOP=.85 EUNHA_DUMP_DT=.85',
        '--time 6.501 --output audit.json':'--time 6.525 --output audit.json',
        '--parent-time .825 --time .826':'--parent-time .825 --time .85',
        'time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1':
            '{ time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1; } 2> runtime.log',
        'audit_kepler_seed_certificate_pair_20261027.py . --output pair_audit.json > pair_audit.log':
            'audit_kepler_seed_certificate_speed_20261027.py . --output speed_audit.json > speed_audit.log',
    }
    for old,new in replacements.items():
        count=2 if old=='HYDRO_TSTOP=.826 EUNHA_DUMP_DT=.826' else 1
        if template.count(old)!=count:raise ValueError(f'Unexpected runner anchor {old!r}')
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
        if root_hashes[name]!=digest(DOCS/name):raise ValueError(f'Changed {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(root_hashes.items()):out.write(f'{value}  {name}\n')
    for case in ('w0','weighted'):
        (RUN/case).mkdir()
        for mode in ('baseline','enabled'):
            source=PARENT/case/mode; target=RUN/case/mode; target.mkdir()
            hashes={}
            for line in (source/'inputs.sha256').read_text().splitlines():
                expected,name=line.split('  ',1)
                if digest(source/name)!=expected:raise ValueError(f'Changed {case}/{mode}/{name}')
                shutil.copy2(source/name,target/name);hashes[name]=digest(target/name)
                if hashes[name]!=expected:raise ValueError('Copied input changed')
            with (target/'inputs.sha256').open('x') as out:
                for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'device_gate':str(PARENT/'pair_audit.json'),
                   'device_gate_sha256':digest(PARENT/'pair_audit.json'),
                   'parent_manifest':manifest,'mode':'diagnostic atomics disabled'},out,indent=2)
        out.write('\n')
    print(RUN)


if __name__=='__main__':main()
