"""Archive and cancel an explicit list of obsolete pending Kepler jobs."""
from datetime import datetime
import json
from pathlib import Path
import re
import socket
import subprocess

DOCS = Path(__file__).resolve().parent
PROJECT = '/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/'
KEEP = {
    438957: 'Current 384^2 halo-recovered native long-time trajectory',
    438959: 'Running same-input one-rank/two-rank A100 comparison',
    437772: 'Running 420^2 inner-resolution control',
    437779: 'Conditional long-time continuation of the running resolution control',
}
GROUPS = {
    'Failed-parent q2 long chain; no executable dependency path':
        [437251, 437868, 437869, 437870, 437871, 437872, 437873, 437874,
         437897, 437909, 437912, 437915, 437916, 437932, 437933, 437934,
         437935, 437936],
    'Failed-parent t6.5 RK4/KDK comparison and its successor':
        [438324, 438326, 438327],
    'Superseded old q1 trajectory or optional execution/scaling variants':
        [438008, 438394, 438364, 438172, 437200, 437181, 436751],
    'Deferred q2 MPI/OpenMP variants while current A100 comparison runs':
        [437892, 437893, 437835, 437836],
    'Optional capacity experiment and its CPU audit; defer while pending':
        [438420, 438970],
}


def remote(args, check=True):
    # The entry point opens one SSH session. Run each archival scheduler
    # command inside that session instead of reconnecting for each job.
    if socket.gethostname().split('.')[0] != 'syntax':
        raise RuntimeError('Scheduler operations must execute inside the syntax session')
    return subprocess.run(args, text=True, capture_output=True, check=check)


def field(line, name):
    match = re.search(r'(?:^| )' + re.escape(name) + r'=([^ ]+)', line)
    if not match:
        raise ValueError(f'Missing {name} in scheduler record')
    return match[1]


def main():
    stamp = datetime.now().astimezone().strftime('%Y-%m-%d_%H%M%S')
    archive = DOCS / ('SLURM_CLEANUP_' + stamp)
    (archive / 'batch_scripts').mkdir(parents=True)
    records = remote(['scontrol', 'show', 'jobs', '-o']).stdout.splitlines()
    desired = {job: reason for reason, jobs in GROUPS.items() for job in jobs}
    metadata = {}
    for line in records:
        if not line.startswith('JobId='):
            continue
        job = int(field(line, 'JobId'))
        if job in desired or job in KEEP:
            if field(line, 'UserId') != 'kjhan(10396)':
                raise ValueError(f'Unexpected owner for job {job}')
            metadata[job] = line
    decisions = []
    cancel = []
    cap_state = field(metadata[438420], 'JobState') if 438420 in metadata else 'MISSING'
    for job, reason in desired.items():
        if job not in metadata:
            decisions.append({'job_id': job, 'action': 'skip_missing_or_terminal', 'reason': reason})
            continue
        line = metadata[job]
        state = field(line, 'JobState')
        script = remote(['scontrol', 'write', 'batch_script', str(job), '-']).stdout
        if not script.strip():
            raise ValueError(f'No batch script could be archived for {job}')
        (archive / 'batch_scripts' / f'{job}.sh').write_text(script)
        workdir = field(line, 'WorkDir')
        if not workdir.startswith(PROJECT):
            if not (job == 436751 and '/gpfs/kjhan/LagForce/' in script):
                raise ValueError(f'Unconfirmed project identity for job {job}: {workdir}')
        if state != 'PENDING':
            KEEP[job] = 'Already running/completing at cleanup; retain in-flight work'
            action = 'keep_in_flight'
        elif job == 438970 and cap_state != 'PENDING':
            KEEP[job] = 'Retain the audit of an already started capacity experiment'
            action = 'keep_audit_of_in_flight_run'
        else:
            action = 'cancel_pending'
            cancel.append(job)
        decisions.append({'job_id': job, 'name': field(line, 'JobName'),
                          'state_before': state, 'workdir': workdir,
                          'reason': reason, 'action': action})
    before = {'recorded_at': datetime.now().astimezone().isoformat(),
              'authorization': 'User explicitly instructed cleanup of unnecessary submitted jobs',
              'kept_jobs': KEEP, 'decisions': decisions,
              'scheduler_records': {str(k): v for k, v in metadata.items()},
              'data_policy': 'Preserve all source, checkpoint, audit and run files; cancel only explicit own pending scheduler jobs'}
    (archive / 'before.json').write_text(json.dumps(before, indent=2) + '\n')

    # Cancel the afterany CPU child first, so cancellation of its parent
    # cannot launch an unnecessary audit. A state filter protects races.
    results = []
    if 438970 in cancel:
        result = remote(['scancel', '--ctld', '--state=PENDING', '--user=kjhan', '438970'], check=False)
        results.append({'jobs': [438970], 'returncode': result.returncode,
                        'stdout': result.stdout, 'stderr': result.stderr})
        if result.returncode:
            raise RuntimeError('CPU audit cancellation failed; leave the parent alone')
        cancel.remove(438970)
    if cancel:
        # Descendants before parents; controller-side state and user filters
        # additionally ensure that no running or other-user job is cancelled.
        cancel_set = set(cancel)
        def depth(job, seen=None):
            seen = set() if seen is None else seen
            if job in seen:
                raise ValueError('Unexpected dependency cycle')
            parents = [int(s) for s in re.findall(r':(\d+)', field(metadata[job], 'Dependency'))
                       if int(s) in cancel_set]
            return 1 + max((depth(p, seen | {job}) for p in parents), default=0)
        cancel.sort(key=depth, reverse=True)
        result = remote(['scancel', '--ctld', '--state=PENDING', '--user=kjhan',
                         *[str(job) for job in cancel]], check=False)
        results.append({'jobs': cancel, 'returncode': result.returncode,
                        'stdout': result.stdout, 'stderr': result.stderr})
    (archive / 'cancellation_requests.json').write_text(json.dumps(results, indent=2) + '\n')
    summary = {'archive': str(archive),
               'cancellation_requests_succeeded': all(r['returncode'] == 0 for r in results),
               'requested_cancel_count': sum(len(r['jobs']) for r in results),
               'kept_jobs': KEEP,
               'requested_cancel_jobs': [job for r in results for job in r['jobs']]}
    (archive / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    if not summary['cancellation_requests_succeeded']:
        raise SystemExit('Some scheduler cancellation requests failed; inspect archived results')


if __name__ == '__main__':
    if socket.gethostname().split('.')[0] == 'syntax':
        main()
    else:
        result = subprocess.run([
            'ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', 'syntax',
            '/home/kjhan/miniconda3/bin/python3', str(Path(__file__).resolve())])
        raise SystemExit(result.returncode)
