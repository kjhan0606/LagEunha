"""Stage one H200 ABBA and wall-impact parity qualification for thermal-only cache."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_ThermalOnlyCacheQual_H200_20261029')
DOCS = Path(__file__).resolve().parent
BASE = Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/candidate_v2/eunha2.exe')
CAND = Path('/gpfs/kjhan/LagForce/Kepler_ThermalCacheOnly_20261029/project/eunha2.exe')
EARLY = Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/gravity_cap_preflight_v1/native_eta006_ada/khout.000022.dat')
IMPACT = Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/const_config_cache_bench_v1/gapped_late_crossing/khout.011580.dat')
PARAM = Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/params.dat')
REF = DOCS/'KEPLER_NATIVE_CONST_CONFIG_CACHE_2026-10-29/runners'
EXPECT = {BASE:'f26816122df101ea71f72b4f3fb9233e965d06784ed94e88cfb46245fa251444',
          CAND:'beca93e14e2e090b7f6c69926474f88070ec7ecbe19a3bafa147488643573a93',
          EARLY:'0f60d80dec2668a7e53a6e99ec84760af824ed70fb7a868ddd605fe53ffe11bc',
          IMPACT:'5946fad7b92ac439a51e01ad469ada3b9055d95d540be70f3c9c5c0b49ccee9c'}


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''): h.update(b)
    return h.hexdigest()


def replace_once(s, old, new):
    if s.count(old) != 1: raise ValueError('Unexpected runner: '+old)
    return s.replace(old, new)


def main():
    if ROOT.exists(): raise FileExistsError(ROOT)
    for p,h in EXPECT.items():
        if sha(p) != h: raise ValueError('Changed input '+str(p))
    ROOT.mkdir()
    stage={'initial':str(EARLY),'initial_sha256':EXPECT[EARLY], 'target_time':.003056807157803711,
           'order':['A1','B1','B2','A2'], 'baseline_binary_sha256':EXPECT[BASE],
           'candidate_manifest':{'eunha2.exe':EXPECT[CAND]},
           'predeclared_gates':{'ABBA_median_wall_per_step_speedup_min':1.05},
           'output_policy':'Six checkpoints ~272 MB total, existing initial snapshots reused; no periodic/full particle dumps.',
           'scope':'Short H200 performance and actual collision parity, not long-time or MPI qualification.'}
    (ROOT/'stage.json').write_text(json.dumps(stage,indent=2)+'\n')
    for name in ['A1','B1','B2','A2','impact_A','impact_B']:
        d=ROOT/name;d.mkdir();impact=name.startswith('impact');candidate=name.startswith('B') or name=='impact_B'
        initial=IMPACT if impact else EARLY
        (d/initial.name).symlink_to(initial.resolve());(d/'eunha2.exe').symlink_to(CAND if candidate else BASE)
        shutil.copy2(PARAM,d/'params.dat')
        template=(REF/('gapped_late_crossing' if impact else 'A1')/'run.sh').read_text()
        olddir='/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/const_config_cache_bench_v1/'+('gapped_late_crossing' if impact else 'A1')
        template=replace_once(template,olddir,str(d))
        if not impact:
            template=replace_once(template,'EUNHA_KEPLER_DTETA=0','EUNHA_KEPLER_DTETA=0.06')
            template=replace_once(template,'HYDRO_TSTOP=0.82845','HYDRO_TSTOP=0.003056807157803711')
            template=replace_once(template,'params.dat 1218','params.dat 22')
        template=replace_once(template,'GFS_NATIVE_CONST_CONFIG='+('1' if impact else '0'),
                              'GFS_NATIVE_CONST_CONFIG='+('1' if candidate else '0'))
        (d/'run.sh').write_text(template)
        if impact:
            impact_stage={'target_time':48.4, 'binary_sha256':EXPECT[CAND if candidate else BASE],
                          'input_sha256':EXPECT[IMPACT], 'predeclared_gates':{
                              'collision_ID_side':[110573,543],
                              'max_step_momentum_and_angular_residual':1e-10,
                              'restart_relative_total_energy_change_max':1e-7}}
            (d/'stage.json').write_text(json.dumps(impact_stage,indent=2)+'\n')
        subprocess.run(['bash','-n',str(d/'run.sh')],check=True)
    audits=[DOCS/n for n in ['audit_kepler_const_config_cache_20261029.py', 'audit_kepler_thermal_cache_impact_20261029.py',
                            'audit_kepler_wall_events_20261029.py', 'audit_kepler_restart_bytes_20261029.py',
                            'kepler_native_global_ids_20261029.py', 'audit_kepler_native_radial_20261026.py']]
    frozen={str(p):sha(p) for p in [*EXPECT,PARAM,ROOT/'stage.json',ROOT/'impact_A/stage.json',ROOT/'impact_B/stage.json',*audits,*[ROOT/n/'run.sh' for n in ['A1','B1','B2','A2','impact_A','impact_B']]]}
    (ROOT/'inputs.sha256').write_text(''.join(h+'  '+p+'\n' for p,h in frozen.items()))
    job='''#!/bin/bash
#SBATCH --job-name=KepThermCacheQual
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir='''+str(ROOT)+'''
#SBATCH --output='''+str(ROOT)+'''/slurm-%j.out
#SBATCH --error='''+str(ROOT)+'''/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
sha256sum -c inputs.sha256
for case_name in A1 B1 B2 A2 impact_A impact_B; do bash "$case_name/run.sh"; done
/home/kjhan/miniconda3/bin/python3 '''+str(DOCS)+'''/audit_kepler_thermal_cache_impact_20261029.py .
/home/kjhan/miniconda3/bin/python3 '''+str(DOCS)+'''/audit_kepler_const_config_cache_20261029.py . --output abba_audit.json
'''
    (ROOT/'job.slurm').write_text(job)
    subprocess.run(['bash','-n',str(ROOT/'job.slurm')],check=True)
    print(ROOT)


if __name__=='__main__': main()
