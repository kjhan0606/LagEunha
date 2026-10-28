"""Install a post-gate-only hook; do not modify a running native evolution."""
import ast
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import socket
import subprocess

DOCS=Path(__file__).resolve().parent
MAIN=Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo40_from_t2654918_a100_20261028')
BENCH=MAIN.parent/'vertex_halo_pair_t2654918_a100_20261029'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def command(*args):
    p=subprocess.run(args,text=True,capture_output=True)
    if p.returncode:raise RuntimeError(p.stdout+p.stderr)
    return p.stdout


def main():
    assert socket.gethostname().split('.')[0]=='syntax'
    receipt_path=MAIN/'checkpoint_validation_install_20261029.json'
    assert not receipt_path.exists(),'Hook already installed; inspect its state instead of reinstalling'
    parent=command('scontrol','show','job','438957','-o')
    assert 'JobState=RUNNING ' in parent and f'WorkDir={MAIN}/to3p0 ' in parent,parent
    stage=MAIN/'to3p0'
    assert not (stage/'audit_gate.json').exists() and not (stage/'continuation_gate.json').exists()
    with (stage/'log').open('rb') as stream:
        stream.seek(max(0,(stage/'log').stat().st_size-12000));tail=stream.read().decode(errors='replace')
    times=re.findall(r'^Time=\s*([0-9.eE+-]+)',tail,re.M)
    assert times and float(times[-1])<3.
    assert sha(stage/'eunha2.exe')=='6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0'
    checksum_before=(stage/'inputs.sha256').read_bytes()
    metadata_before=(stage/'stage.json').read_bytes()
    assert json.loads(metadata_before)['next_job']==str(MAIN/'to4p0/job.slurm')
    queued=command('scontrol','show','job','439115','-o')
    assert 'JobState=PENDING ' in queued and 'RunTime=00:00:00 ' in queued and f'WorkDir={BENCH} ' in queued,queued
    cases=json.loads((BENCH/'combo_staging_receipt_20261029.json').read_text())['cases']
    assert len(cases)==7 and not any((BENCH/case/'log').exists() for case in cases)
    source=DOCS/'run_kepler_checkpoint_validation_20261029.py'
    ast.parse(source.read_text())
    test=subprocess.run(['/home/kjhan/miniconda3/bin/python3',str(DOCS/'test_kepler_checkpoint_validation_20261029.py')],capture_output=True,text=True)
    assert test.returncode==0,test.stdout+test.stderr
    frozen={'audit_inputs.sha256':sha(BENCH/'audit_inputs.sha256')}
    for line in (BENCH/'audit_inputs.sha256').read_text().splitlines():
        expected,name=line.split('  ',1)
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        assert sha(BENCH/name)==expected
        frozen[name]=expected
    for case in cases:
        for name in ['run.sh','inputs.sha256']:
            frozen[case+'/'+name]=sha(BENCH/case/name)
    config={'parent_job':438957,'queued_job':439115,'stage':str(stage),'benchmark':str(BENCH),
            'expected_next_job':str(MAIN/'to4p0/job.slurm'),'cases':cases,'frozen_hashes':frozen,
            'automatic_promotion':False}
    config_path=MAIN/'checkpoint_validation_config_20261029.json'
    with config_path.open('x') as stream:json.dump(config,stream,indent=2);stream.write('\n')
    runner=MAIN/source.name
    assert not runner.exists();shutil.copy2(source,runner)
    checker=MAIN/'check_radial_and_next.py'
    before=checker.read_text()
    assert 'post_checkpoint_validation' not in before
    assert "assert len(checks)==9 and all(checks.values())" in before
    backup=MAIN/'check_radial_and_next.before_checkpoint_validation_20261029.py'
    with backup.open('x') as stream:stream.write(before)
    suffix=f'''
# The native MPI program and 15+9 gates finish before this optional dispatch.
# Only the original t=3 allocation can use its still-reserved GPU. No change
# to the next simulation, its binary, parameters, stop time or restart state.
import os as _post_os, subprocess as _post_subprocess, sys as _post_sys
if str(Path.cwd().resolve())=={str(stage)!r} and _post_os.environ.get('SLURM_JOB_ID')=='438957':
    _post_result=_post_subprocess.run([_post_sys.executable,{str(runner)!r},{str(config_path)!r}])
    if _post_result.returncode:
        raise RuntimeError('Auxiliary validation cleanup/interruption requires review; accepted t=3 checkpoint is preserved, continuation metadata unchanged')
'''
    after=before+suffix
    ast.parse(after)
    temporary=checker.with_suffix('.next.py');temporary.write_text(after);temporary.replace(checker)
    assert (stage/'inputs.sha256').read_bytes()==checksum_before
    assert (stage/'stage.json').read_bytes()==metadata_before
    result={'recorded_at':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
            'parent_job':438957,'standalone_request':439115,'installed_before_native_time':float(times[-1]),
            'checker':str(checker),'old_checker_sha256':sha(backup),'new_checker_sha256':sha(checker),
            'runner':str(runner),'runner_sha256':sha(runner),'config':str(config_path),'config_sha256':sha(config_path),
            'trigger':'After exact t=3 and native 15+radial 9 checks pass, original numeric MPI steps gone, original allocation still active and >=45 minutes remain',
            'standalone_policy':'If pending/unstarted, hold, recheck, atomic claim, pending-only cancellation and confirm cancellation before using frozen batch in existing allocation; running/completed/partial workloads skipped',
            'resource_policy':'Existing one GPU/two CPU/32 GiB allocation; no new allocation and no concurrency with original native simulation. Prior one/two-rank native MaxRSS was about 1.04/1.81 GiB per task.',
            'deadline_policy':'At most 40 minutes, absolute allocation end with cleanup reserve; terminate auxiliary group and only numeric auxiliary Slurm steps; never cancel parent/batch/extern/other jobs',
            'continuation_policy':'No automatic promotion; normal auxiliary failure with confirmed cleanup keeps original next job. Unverified residual ranks or interruption causes an operational stop at the preserved checkpoint.',
            'native_binary_parameters_IC_stop_restart_metadata_unchanged':True,
            'automatic_promotion':False,'tests':test.stdout+test.stderr,
            'astra_decision':'Inline validation approved in principle; automatic promotion deferred pending direct old-production performance evidence and full continuation contract review'}
    with receipt_path.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    (DOCS/'KEPLER_CHECKPOINT_VALIDATION_INSTALL_2026-10-29.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
