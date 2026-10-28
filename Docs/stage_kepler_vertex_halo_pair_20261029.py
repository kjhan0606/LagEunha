"""One native GPU job: executable qualifier and independent halo-geometry gate."""
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

DOCS=Path(__file__).resolve().parent
Q1=Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1')
ROOT=Q1/'vertex_halo_pair_t2654918_a100_20261029'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    build_root=Q1/'vertex_halo_native_build_20261029'
    build=json.loads((build_root/'build_result.json').read_text())
    candidate=Path(build['binary'])
    assert sha(candidate)==build['sha256']=='a22e0e8ebea0eca807ab3762d95bf037184b7784218a42fd2eb1ae7a24526bfc'
    source=Q1/'tess_padding_reuse_pair_20261028/np1_control'
    old=(Q1/'tess_padding_reuse_pair_20261028/qualifier_old/eunha2.exe').resolve()
    assert sha(old)=='6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'
    checkpoint=(source/'khout.023008.dat').resolve()
    assert sha(checkpoint)=='c34a846dedaf2e01a8e5548de2aa091f98db53f761424ce6cd931107ff0fc4a4'
    ROOT.mkdir()
    for name in ['audit_kepler_vertex_halo_pair_20261029.py','audit_kepler_geometry_force_dump_20261029.py',
                 'audit_kepler_padding_reuse_pair_20261028.py','audit_kepler_polar_long_gate_20261027.py']:
        shutil.copy2(DOCS/name,ROOT/name)
    start=2.6549183918680868
    cases=[('qualifier_old',old,0,40,start+1e-5,0),
           ('qualifier_new',candidate,0,40,start+1e-5,0),
           ('control40',candidate,0,40,2.66,4),
           ('vertex40',candidate,1,40,2.66,4),
           ('vertex20',candidate,1,20,2.66,4)]
    for name,binary,enabled,halo,target,dumps in cases:
        run=ROOT/name;run.mkdir()
        inputs={'eunha2.exe':binary,'params.dat':(source/'params.dat').resolve(),
                'khout.000000.dat':(source/'khout.000000.dat').resolve(),
                'khout.023008.dat':checkpoint}
        for helper in ['audit_kepler_polar_long_gate_20261027.py','audit_kepler_native_radial_20261026.py']:
            inputs[helper]=(source/helper).resolve()
        for dest,original in inputs.items(): (run/dest).symlink_to(original)
        script=(source/'run.sh').read_text()
        assert script.count('GFS_KEPLER_HALO_CELLS=40')==1
        script=script.replace('GFS_KEPLER_HALO_CELLS=40',f'GFS_KEPLER_HALO_CELLS={halo}')
        script=script.replace('HYDRO_TSTOP=2.6600000000000001',f'HYDRO_TSTOP={target:.17g}')
        script=script.replace('--time 2.6600000000000001',f'--time {target:.17g}')
        assert 'GFS_REUSE_TESS_PADDING=0' in script and 'GFS_GPU_LOCAL_PERIODIC_SYNC=0' in script
        script=script.replace("TIMEFORMAT=",f'export GFS_GPU_VERTEX_HALO_CERT={enabled} GFS_GEOM_FORCE_DUMP_MAXCALL={dumps}\nTIMEFORMAT=',1)
        (run/'run.sh').write_text(script)
        subprocess.run(['bash','-n',str(run/'run.sh')],check=True)
        (run/'inputs.sha256').write_text(''.join(f'{sha(path)}  {dest}\n' for dest,path in inputs.items())
                                       +f'{sha(run/"run.sh")}  run.sh\n')
    job=f'''#!/usr/bin/env bash
#SBATCH --job-name=KepVertexHaloPair
#SBATCH --partition=a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --chdir={ROOT}
#SBATCH --output={ROOT}/slurm-%j.out
#SBATCH --error={ROOT}/slurm-%j.err
set -euo pipefail
sha256sum -c audit_inputs.sha256
for name in qualifier_old qualifier_new control40 vertex40 vertex20; do
  cd {ROOT}/$name
  bash run.sh > runner.out 2> runner.err
  if [[ "$name" == qualifier_new ]]; then
    /home/kjhan/miniconda3/bin/python3 {ROOT}/audit_kepler_vertex_halo_pair_20261029.py --qualifier > {ROOT}/qualifier.log
  fi
done
/home/kjhan/miniconda3/bin/python3 {ROOT}/audit_kepler_vertex_halo_pair_20261029.py > {ROOT}/comparison.log
'''
    (ROOT/'job.slurm').write_text(job)
    subprocess.run(['bash','-n',str(ROOT/'job.slurm')],check=True)
    audit_files=[ROOT/p.name for p in ROOT.iterdir() if p.suffix=='.py']+[ROOT/'job.slurm']
    (ROOT/'audit_inputs.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in audit_files))
    receipt={'recorded_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
             'root':str(ROOT),'start':start,'target':2.66,'checkpoint':str(checkpoint),
             'cases':[{'name':n,'binary_sha256':sha(b),'vertex_certificate':e,'halo':h,'target':t,'read_only_force_dumps':d} for n,b,e,h,t,d in cases],
             'native_binary':build,'qualifier':'Canonical executable vs candidate with all new options off; exact fields and accepted-step trace',
             'gate':'Exact ten endpoint fields, accepted-step traces, owned/face counts, all four RK4 per-cell force/work/timestep arrays and complete face geometry plus normalized neighbor/image identities',
             'safety':'Existing 15 physics checks and every energy ledger injection/count zero; actual-padding domain banners required; reduced halo must exercise the new certificate',
             'relationship_to_job_439029':'Independent experiment: this changes geometric halo certification/width; 439029 measures redundant host padding and MPI/local alias refresh. New feature stays off in that job.',
             'storage_policy':'Five endpoint dumps about 425 MB plus twelve binary face/force audits about 650 MB; parent states linked. About 1.2 GB total including executables/logs.',
             'resources':{'partition':'a100_pcie','GPU':1,'CPUs':2,'memory_GB':48,'minutes':40},
             'production_modified':False,'production_promoted':False}
    (ROOT/'staging_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    result=subprocess.run(['ssh','syntax','sbatch','--parsable',str(ROOT/'job.slurm')],capture_output=True,text=True)
    (ROOT/'submission.stdout').write_text(result.stdout)
    (ROOT/'submission.stderr').write_text(result.stderr)
    result.check_returncode()
    job_id=result.stdout.strip().split(';')[0];assert job_id.isdecimal(),result.stdout
    receipt['job_id']=int(job_id)
    (ROOT/'staging_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (DOCS/'KEPLER_VERTEX_HALO_PAIR_STAGE_2026-10-29.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'job_id':job_id,'root':str(ROOT),'cases':len(cases),'GPU':1},indent=2))


if __name__=='__main__':main()
