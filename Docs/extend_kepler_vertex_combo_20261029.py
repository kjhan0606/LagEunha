"""Safely replace the unstarted five-case GPU request with seven cases."""
import ast
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import socket
import subprocess

from stage_kepler_vertex_halo_pair_20261029 import ROOT, job_script

DOCS=Path(__file__).resolve().parent
OLD='439065'
ORIGINAL=['qualifier_old','qualifier_new','control40','vertex40','vertex20']
ADDED=['vertex20_reuse','vertex20_reuse_local_sync']


def run(*args):
    p=subprocess.run(args,capture_output=True,text=True)
    if p.returncode: raise RuntimeError(f'{args}\n{p.stdout}\n{p.stderr}')
    return p.stdout


def sha(data): return hashlib.sha256(data).hexdigest()


def main():
    assert socket.gethostname().split('.')[0]=='syntax'
    receipt_path=ROOT/'combo_staging_receipt_20261029.json'
    assert not receipt_path.exists(), 'Inspect existing submission; do not duplicate it'
    before=run('scontrol','show','job',OLD,'-o')
    assert 'JobState=PENDING ' in before and 'RunTime=00:00:00 ' in before,before
    assert not any((ROOT/name/'log').exists() for name in ORIGINAL+ADDED)
    test=subprocess.run(['/home/kjhan/miniconda3/bin/python3',str(DOCS/'test_kepler_vertex_combo_gate_20261029.py')],capture_output=True,text=True)
    assert test.returncode==0,test.stdout+test.stderr
    source=(ROOT/'vertex20/run.sh').read_text()
    source_inputs=(ROOT/'vertex20/inputs.sha256').read_text()
    workers={}
    for name in ORIGINAL+ADDED:
        script=(ROOT/name/'run.sh').read_text() if name in ORIGINAL else source
        original_script=script
        if name in ADDED:
            assert script.count('GFS_REUSE_TESS_PADDING=0')==1
            script=script.replace('GFS_REUSE_TESS_PADDING=0','GFS_REUSE_TESS_PADDING=1')
            if name.endswith('local_sync'):
                assert script.count('GFS_GPU_LOCAL_PERIODIC_SYNC=0')==1
                script=script.replace('GFS_GPU_LOCAL_PERIODIC_SYNC=0','GFS_GPU_LOCAL_PERIODIC_SYNC=1')
        assert 'record_kepler_vertex_case_context' not in script
        script=script.replace('TIMEFORMAT=',f'"$PY" {ROOT}/record_kepler_vertex_case_context_20261029.py\nTIMEFORMAT=',1)
        checks=(ROOT/name/'inputs.sha256').read_text() if name in ORIGINAL else source_inputs
        old_hash=sha(original_script.encode())
        assert checks.count(old_hash+'  run.sh')==1
        checks=checks.replace(old_hash+'  run.sh',sha(script.encode())+'  run.sh')
        workers[name]=(script.encode(),checks.encode())
    new_job=ROOT/'job_with_combinations.slurm'
    with new_job.open('x') as stream:stream.write(job_script(ROOT))
    run('bash','-n',str(new_job))
    for name in ADDED:
        folder=ROOT/name;folder.mkdir()
        for filename in ['eunha2.exe','params.dat','khout.000000.dat','khout.023008.dat',
                         'audit_kepler_polar_long_gate_20261027.py','audit_kepler_native_radial_20261026.py']:
            (folder/filename).symlink_to((ROOT/'vertex20'/filename).resolve())
    # Stage syntax-checked worker text separately before holding the job.
    prepared=ROOT/'prepared_combo_workers_20261029';prepared.mkdir()
    for name,(script,checks) in workers.items():
        (prepared/(name+'.sh')).write_bytes(script)
        run('bash','-n',str(prepared/(name+'.sh')))
    for name in ['audit_kepler_vertex_halo_pair_20261029.py','record_kepler_vertex_case_context_20261029.py']:
        ast.parse((DOCS/name).read_text())
    archive=ROOT/'unstarted_five_case_20261029';archive.mkdir()
    backup={ROOT/'audit_inputs.sha256':(ROOT/'audit_inputs.sha256').read_bytes(),
            ROOT/'audit_kepler_vertex_halo_pair_20261029.py':(ROOT/'audit_kepler_vertex_halo_pair_20261029.py').read_bytes()}
    for name in ORIGINAL:
        for filename in ['run.sh','inputs.sha256']:
            path=ROOT/name/filename;backup[path]=path.read_bytes()
    for path,data in backup.items():
        dest=archive/path.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
    shutil.copy2(ROOT/'job.slurm',archive/'original_job.slurm')
    shutil.copy2(ROOT/'staging_receipt.json',archive/'staging_receipt.json')
    cancelled=False
    run('scontrol','hold',OLD)
    try:
        held=run('scontrol','show','job',OLD,'-o')
        assert 'JobState=PENDING ' in held and 'RunTime=00:00:00 ' in held,held
        run('scontrol','write','batch_script',OLD,str(archive/'submitted_batch.slurm'))
        for name,(script,checks) in workers.items():
            (ROOT/name/'run.sh').write_bytes(script)
            (ROOT/name/'inputs.sha256').write_bytes(checks)
        for name in ['audit_kepler_vertex_halo_pair_20261029.py','record_kepler_vertex_case_context_20261029.py']:
            shutil.copy2(DOCS/name,ROOT/name)
        audit_files=[ROOT/name for name in ['audit_kepler_vertex_halo_pair_20261029.py','audit_kepler_geometry_force_dump_20261029.py',
                     'audit_kepler_padding_reuse_pair_20261028.py','audit_kepler_polar_long_gate_20261027.py',
                     'record_kepler_vertex_case_context_20261029.py','job_with_combinations.slurm']]
        (ROOT/'audit_inputs.sha256').write_text(''.join(f'{sha(p.read_bytes())}  {p.name}\n' for p in audit_files))
        run('scancel',OLD);cancelled=True
        submitted=run('sbatch','--parsable',str(new_job)).strip().split(';')[0]
        assert submitted.isdecimal(),submitted
    except BaseException:
        if not cancelled:
            for path,data in backup.items():path.write_bytes(data)
            run('scontrol','release',OLD)
        raise
    receipt={'recorded_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
             'old_job':int(OLD),'new_job':int(submitted),'before':before.strip(),'held_state':held.strip(),
             'archive':str(archive),'job_script':str(new_job),'job_script_sha256':sha(new_job.read_bytes()),
             'cases':ORIGINAL+ADDED,'native_binary_sha256':sha((ROOT/'vertex20/eunha2.exe').read_bytes()),
             'native_build_repeated':False,'pending_GPU_job_count_change':0,
             'reason':'Astra recommends two sequential controlled combinations to validate single-GPU acceleration without waiting for the independent two-GPU allocation',
             'physics_or_time_window_changed':False,'checkpoint_start':2.6549183918680868,'full_case_target':2.66,
             'gate':'Exact first four RK4 complete face/support/image and force/work arrays plus strict endpoint and every zero-injection ledger; direct predecessor controls; all seven cases one allocation with exact requested flags',
             'time_guard':'40 minutes total; require 450 seconds for each full case or 90 for a qualifier plus 180 seconds comparison reserve; insufficient time records incomplete and exits without shortening evolution',
             'storage_policy':'About 1.8 GB: seven endpoints plus twenty full face/force dumps; parents linked',
             'tests':test.stdout+test.stderr,'long_run_and_its_continuation_modified':False,'production_promoted':False}
    with receipt_path.open('x') as stream:json.dump(receipt,stream,indent=2);stream.write('\n')
    (DOCS/'KEPLER_VERTEX_HALO_COMBO_STAGE_2026-10-29.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['before','held_state']},indent=2))


if __name__=='__main__':main()
