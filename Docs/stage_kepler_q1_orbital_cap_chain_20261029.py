"""Stage, but never submit, one native H200 orbital-cap checkpoint chain."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stage(root, mode):
    docs=Path(__file__).resolve().parent
    parent=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029')
    preflight=parent/'gravity_cap_preflight_v1';screen=preflight/'gravity_cap_screen.json'
    contract=preflight/'contract.json';initial=parent/'ic_v1/khout.000001.dat'
    startup=preflight/'native_eta006_ada';startup_audit=startup/'audit.json'
    startup_core=startup/'core_screen.json';startup_result=json.loads(startup_audit.read_text())
    starting=Path(startup_result['endpoint']);starting_step=int(starting.name.split('.')[1])
    if mode not in ('slurm','ada-baseline','ada-cache'):raise ValueError('Choose slurm, ada-baseline or ada-cache')
    cache_parity=preflight/'eta006_cache_parity.json'
    use_cache=mode=='ada-cache';use_ada=mode.startswith('ada-')
    if use_cache:
        exe=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_const_cache_v1/eunha2.exe')
        expected_exe='7364a05c09eeb393d4e8528f01a82efb7f6d9c6fba9b394d2d1e62e6bd78345f'
        parity=json.loads(cache_parity.read_text())
        if not parity['passed'] or parity['candidate_binary_sha256']!=expected_exe:
            raise ValueError('Eta-specific cache parity failed')
    else:
        exe=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v2/eunha2.exe')
        expected_exe='f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444'
    baseline=parent/'chain_h200_v2/to10p0';validator=parent/'post_t10_restart_dt_h200_v1'
    result=json.loads(screen.read_text())
    if not result['passed'] or result['selected_eta']!=.06 or result['restart_selection']['selected_source']!=str(initial):
        raise ValueError('Different gravity or restart preflight result')
    if sha(contract)!=result['contract_sha256'] or sha(initial)!=result['restart_selection']['selected_source_sha256']:
        raise ValueError('Preregistered contract or preserved original IC changed')
    if (not startup_result['passed'] or startup_result['orbital_eta']!=.06
            or not json.loads(startup_core.read_text())['passed']
            or sha(starting)!=startup_result['endpoint_sha256']):
        raise ValueError('Existing capped native startup is not an accepted core-checked restart')
    if sha(exe)!=expected_exe:raise ValueError('Native executable changed')
    if root.exists():raise FileExistsError('Never overwrite an existing experiment')
    root.mkdir(parents=True)
    sources=['run_kepler_q1_orbital_cap_chain_20261029.py','audit_kepler_q1_checkpoint_geometry_20261029.py',
             'audit_kepler_q1_restart_dt_bundle_20261029.py','audit_kepler_q1_gravity_cap_preflight_20261029.py',
             'kepler_native_global_ids_20261029.py','audit_kepler_native_radial_20261026.py',
             'make_kepler_bounded_polygon_ic_20261027.py','audit_kepler_q1_core_rk4_control_20261029.py',
             'audit_kepler_q1_bounded_gate_20261029.py','audit_kepler_const_config_cache_20261029.py',
             'audit_kepler_restart_bytes_20261029.py']
    for name in sources:shutil.copy2(docs/name,root/name)
    for source,name in [(contract,'contract.json'),(screen,'gravity_cap_screen.json'),
                        (baseline/'params.dat','params.dat'),(baseline/'run.sh','native_template.sh')]:
        shutil.copy2(source,root/name)
    validator_cfg=json.loads((validator/'bundle.json').read_text())
    for path,expected in validator_cfg['frozen_hashes'].items():
        if sha(path)!=expected:raise ValueError('Existing environment validator source changed: '+path)
    artifacts=[root/name for name in sources]+[root/name for name in ['contract.json','gravity_cap_screen.json','params.dat','native_template.sh']]
    external=[initial,exe,starting,startup_audit,startup_core,startup/'independent_geometry.json',
              startup/'checkpoint_gate.json',startup/'preflight.json',startup/'run.sh',startup/'log',
              Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py'),
              Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/snapshot_offsets.json')]
    if use_cache:external.append(cache_parity)
    artifacts += external+[Path(path) for path in validator_cfg['frozen_hashes']]
    cfg={'initial':str(initial),'binary':str(exe),'eta':.06,'gravity_screen':str(root/'gravity_cap_screen.json'),
         'execution_mode':'ada' if use_ada else 'slurm','cache_mode':int(use_cache),'cache_parity':str(cache_parity) if use_cache else None,
         'ada_gpu_uuid':'GPU-42b3367f-1df7-2b2d-7961-7781c87b022d',
         'starting_checkpoint':str(starting),'starting_step':starting_step,'startup_audit':str(startup_audit),
         'startup_core_screen':str(startup_core),
         'environment_validator':str(validator/'run_kepler_q1_restart_dt_bundle_20261029.py'),
         'template_directory':str(baseline),'stages':[{'directory':'to0p4125','time':.4125},
                 {'directory':'to0p825','time':.825}],
         'frozen_hashes':{str(p):sha(p) for p in artifacts},
         'output_policy':'Two final 52 MB checkpoints, no periodic, full-particle TSV, face or tessellation dumps. First is retained for recovery. All original and core screens gate each stage.',
         'scope':'One finite W0 native Q1 accuracy branch, starting from the accepted eta=.06 startup derived from the unchanged saved original IC. Original current H200 chain is independent and continues. No GIZMO long-time success claim.'}
    (root/'chain.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (root/'chain.sha256').write_text(sha(root/'chain.json')+'  chain.json\n')
    job=f'''#!/bin/bash
#SBATCH --job-name=KepQ1GravityCap
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --chdir={root}
#SBATCH --output={root}/slurm-%j.out
#SBATCH --error={root}/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
sha256sum -c chain.sha256
/home/kjhan/miniconda3/bin/python3 {root}/run_kepler_q1_orbital_cap_chain_20261029.py {root}
'''
    if mode=='slurm':
        job_file=root/'job.slurm'
    else:
        job_file=root/'job_ada.sh'
        job='\n'.join(line for line in job.splitlines() if not line.startswith('#SBATCH'))+'\n'
        job=job.replace('set -euo pipefail\n','set -euo pipefail\ncd '+str(root)+'\nunset SLURM_JOB_ID SLURM_CPUS_PER_TASK I_MPI_HYDRA_BOOTSTRAP\nexport CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d\n')
    job_file.write_text(job)
    subprocess.run(['bash','-n',str(job_file)],check=True)
    print(json.dumps({'root':str(root),'staged_only':True,'eta':.06,'starting_checkpoint':str(starting),'stages':cfg['stages'],
          'snapshots':2,'expected_total_snapshot_bytes':104085640,
          'hardware':'one H200, 8 OpenMP CPUs, 32 GB, maximum 6 hours' if mode=='slurm' else 'one RTX5000 Ada and 8 OpenMP CPUs',
          'chain_sha256':sha(root/'chain.json'),'job_sha256':sha(job_file)},indent=2))


if __name__=='__main__':stage(Path(sys.argv[1]).resolve(),sys.argv[2])
