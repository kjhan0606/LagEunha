"""Replace one unstarted benchmark request with its seven-case revision."""
import ast
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import socket
import subprocess

DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/tess_padding_reuse_pair_20261028')
OLD_JOB = '439001'


def command(*args):
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if socket.gethostname().split('.')[0] != 'syntax':
        raise SystemExit('Run on syntax, the Slurm authority')
    receipt_path = ROOT / 'local_sync_extension_20261029.json'
    if receipt_path.exists():
        raise FileExistsError('Extension already recorded; do not submit a duplicate')
    before = command('scontrol', 'show', 'job', OLD_JOB, '-o')
    assert 'JobState=PENDING ' in before and 'RunTime=00:00:00 ' in before, before
    names = ['qualifier_old', 'qualifier_new', 'np1_control', 'np2_control', 'np2_reuse', 'np1_reuse']
    assert not any((ROOT / name / 'log').exists() for name in names)
    source = ROOT / 'np1_reuse'
    run = ROOT / 'np1_reuse_local_sync'
    run.mkdir(exist_ok=True)
    assert not (run / 'log').exists()
    for name in ['eunha2.exe', 'params.dat', 'khout.000000.dat', 'khout.023008.dat',
                 'audit_kepler_polar_long_gate_20261027.py', 'audit_kepler_native_radial_20261026.py']:
        link = run / name
        if link.exists() or link.is_symlink():
            assert link.is_symlink() and link.resolve() == (source / name).resolve()
        else:
            link.symlink_to((source / name).resolve())
    if (run / 'inputs.sha256').exists():
        assert (run / 'inputs.sha256').read_bytes() == (source / 'inputs.sha256').read_bytes()
    else:
        shutil.copy2(source / 'inputs.sha256', run / 'inputs.sha256')
    script = (source / 'run.sh').read_text()
    assert script.count('GFS_GPU_LOCAL_PERIODIC_SYNC=0') == 1
    assert 'GFS_REUSE_TESS_PADDING=1' in script
    script = script.replace('GFS_GPU_LOCAL_PERIODIC_SYNC=0', 'GFS_GPU_LOCAL_PERIODIC_SYNC=1')
    if (run / 'run.sh').exists():
        assert (run / 'run.sh').read_text() == script
    else:
        with (run / 'run.sh').open('x') as stream:
            stream.write(script)
    command('bash', '-n', str(run / 'run.sh'))
    subprocess.run(['sha256sum', '-c', 'inputs.sha256'], cwd=run, check=True, capture_output=True)
    job = (ROOT / 'job.slurm').read_text()
    old_cases = 'qualifier_old qualifier_new np1_control np2_control np2_reuse np1_reuse; do'
    assert job.count(old_cases) == 1
    job = job.replace(old_cases, old_cases.replace('; do', ' np1_reuse_local_sync; do'))
    job_path = ROOT / 'job_with_local_sync.slurm'
    if job_path.exists():
        assert job_path.read_text() == job
    else:
        with job_path.open('x') as stream:
            stream.write(job)
    command('bash', '-n', str(job_path))
    auditor_name = 'audit_kepler_padding_reuse_pair_20261028.py'
    auditor = ROOT / auditor_name
    ast.parse((DOCS / auditor_name).read_text())
    test = subprocess.run(['/home/kjhan/miniconda3/bin/python3',
                           str(DOCS / 'test_kepler_padding_reuse_audit_20261029.py')],
                          check=True, capture_output=True, text=True)
    # Hold first: a scheduler transition cannot start the old six-case job
    # while its shared auditor is being revised. Existing output is preserved.
    command('scontrol', 'hold', OLD_JOB)
    held = command('scontrol', 'show', 'job', OLD_JOB, '-o')
    assert 'JobState=PENDING ' in held and 'RunTime=00:00:00 ' in held, held
    prior_script = ROOT / 'prior_job_439001.slurm'
    assert not prior_script.exists()
    command('scontrol', 'write', 'batch_script', OLD_JOB, str(prior_script))
    backup = ROOT / 'audit_kepler_padding_reuse_pair_20261028.before_local_sync.py'
    with backup.open('xb') as stream:
        stream.write(auditor.read_bytes())
    shutil.copy2(DOCS / auditor_name, ROOT / (auditor_name + '.next'))
    (ROOT / (auditor_name + '.next')).replace(auditor)
    command('scancel', OLD_JOB)
    submitted = command('sbatch', '--parsable', str(job_path)).strip()
    new_job = submitted.split(';')[0]
    assert new_job.isdecimal(), submitted
    receipt = {
        'recorded_at': datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        'old_job': int(OLD_JOB), 'new_job': int(new_job),
        'old_state_before_replacement': before.strip(), 'old_state_held': held.strip(),
        'replacement_reason': 'Add one late-checkpoint verification of the existing GPU periodic alias refresh, with the same binary and geometry',
        'job_script': str(job_path), 'job_script_sha256': digest(job_path),
        'new_case': str(run), 'new_case_runner_sha256': digest(run / 'run.sh'),
        'auditor_sha256': digest(auditor), 'old_auditor_sha256': digest(backup),
        'native_binary_sha256': digest(run / 'eunha2.exe'),
        'tests': test.stdout + test.stderr,
        'existing_long_run_modified': False,
        'extra_queued_jobs': 0,
        'storage_policy': 'Seven endpoints, about 595 MB total; only one extra 85 MB snapshot; existing parent inputs linked',
        'production_promoted': False,
    }
    with receipt_path.open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    (DOCS / 'KEPLER_PADDING_REUSE_LOCAL_SYNC_STAGE_2026-10-29.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
