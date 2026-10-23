#!/usr/bin/env python3
"""User-authorized one-time project stop and Git publication before power loss."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import subprocess
import traceback
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

CONTROL = Path(__file__).resolve().parent
CFG = json.loads((CONTROL / 'config.json').read_text())
ROOT = Path(CFG['root'])
PUB = Path(CFG['publish'])


def now():
    return datetime.now().astimezone().isoformat()


def save(name, data):
    path = CONTROL / name
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w') as out:
        json.dump(data, out, ensure_ascii=False, indent=2)
        out.write('\n')
        out.flush()
        os.fsync(out.fileno())
    tmp.replace(path)


def run(args, timeout=120):
    print(now(), shlex.join(args), flush=True)
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=timeout)
    if p.returncode:
        raise RuntimeError(f'{shlex.join(args)}: exit={p.returncode}\n{p.stdout}\n{p.stderr}')
    return p.stdout.strip()


def remote(args, timeout=120):
    return run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10',
                'syntax', shlex.join(args)], timeout)


def git(*args, timeout=120):
    return remote(['git', '-C', str(PUB), *args], timeout)


def preflight():
    import local_scope
    assert str(PUB) == git('rev-parse', '--show-toplevel')
    assert git('remote', 'get-url', 'origin') == CFG['remote']
    return dict(observed_at=now(), mode='read_only_no_stop_no_commit',
                deadline=CFG['deadline'], host=os.uname().nodename,
                local_project_processes=local_scope.inventory(),
                slurm=json.loads(remote(['/usr/bin/python3', str(CONTROL / 'slurm_scope.py'), 'inventory'])),
                git_status=git('status', '--porcelain'),
                local_head=git('rev-parse', 'HEAD'),
                remote_master=git('ls-remote', 'origin', 'refs/heads/master'))


def stop_work():
    import local_scope
    outcomes = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        work = {'local': pool.submit(local_scope.stop),
                'slurm': pool.submit(remote, ['/usr/bin/python3', str(CONTROL / 'slurm_scope.py'), 'stop'], 180)}
        for name, future in work.items():
            try:
                value = future.result()
                outcomes[name] = json.loads(value) if name == 'slurm' else value
            except Exception as exc:
                outcomes[name] = {'error': str(exc)}
    save('stop_result.json', outcomes)
    return outcomes


def handover(outcomes):
    evidence = ROOT / 'material_cell_reference_continuation_v1_h200_mpi4/run/continuation_result.json'
    native = json.loads(evidence.read_text()) if evidence.exists() else {}
    accepted = ROOT / 'material_cell_reference_continuation_v1_h200_mpi4/run/continuous_step7'
    checkpoints = [dict(path=str(p), bytes=p.stat().st_size) for p in sorted(accepted.glob('accepted_r*.bin'))]
    save('accepted_checkpoint_inventory.json', dict(observed_at=now(), checkpoints=checkpoints, continuation=native))
    content = f'''# Kepler power outage handover — 2026-10-23

## User instruction and stop

Stop this project's work at 10:50 Asia/Seoul, commit and push before the announced
11:00 server shutdown, then cease work. No operating-system shutdown is performed.
Recorded at: {now()}.
Only owned processes/jobs bound to `/gpfs/kjhan/LagForce/Kepler_` are selected.
No new simulation, restart, requeue or automatic continuation is authorized after
this stop until the user explicitly resumes. Preserve partial files and all raw data.

Stop observations (empty lists mean no matching work was running):

```json
{json.dumps(outcomes, ensure_ascii=False, indent=2)}
```

## Latest native continuation, superseding earlier pending/resource-wait reports

Job 435137 completed two new accepted native steps (6 and 7), with independent
actual v3 restart readbacks. Latest accepted time is
`{native.get('accepted_time', 'MISSING')}`; total accepted steps
`{native.get('total_accepted_steps', 'MISSING')}`.
The prepared genuine step-6 selected k1 was consumed exactly once; the old full
shadow and accepted intervals were not repeated. Do not reuse that consumed k1.

Latest accepted rank checkpoints remain on shared GPFS:

```json
{json.dumps(checkpoints, indent=2)}
```

Complete step-7 wall time including readback was
`{native.get('step7_complete_wall_seconds_including_readback')}` seconds.
Long release: `{native.get('long_run_released')}`.
Stop reason: `{native.get('stop_reason')}`.
The measured cost gate failed; t>6–10 Kepler stability has NOT been demonstrated.
Global map injectivity guard is also incomplete. Historical step-5/pending files
retain their observation times and do not override this actual continuation result.

Before resuming, inspect latest checkpoint/source/operator bindings and any
geometry precomputation binding. Do not silently replay accepted intervals or
apply old-state geometry to the accepted step-7 state. Preserve the consumed
`material_step6_reference_continuation.claim` and its terminal record.

## Code and data preservation

Authoritative research root: `{ROOT}`.
Original parallel solver source: `{CFG['original_source']}`.
Isolated publication worktree: `{PUB}` (Git commands run on syntax).
Remote: `{CFG['remote']}`, intended branch `{CFG['remote_branch']}`.

The commit contains source, small JSON/log evidence, exact frozen source closures,
and this handover. Native binaries, large quadrature/geometry/factor data and restart
checkpoints stay on GPFS; they are not deleted or uploaded as Git objects.
`Research/Kepler/SNAPSHOT_SOURCE_MANIFEST.json` binds copied source bytes.
`SNAPSHOT_CONTENTS.sha256` indexes the archived files.

The final commit/push hash and verification receipt are written AFTER publication
to `{CONTROL / 'publish_receipt.json'}`; this document does not claim a future push
has already succeeded. Read that receipt and `{CONTROL / 'result.json'}` before resuming.
No new physical validation was launched by the shutdown procedure.
'''
    (ROOT / 'HANDOVER_POWER_OUTAGE_20261023.md').write_text(content)
    (PUB / 'Docs/POWER_OUTAGE_HANDOVER_2026-10-23.md').write_text(content)
    old = PUB / 'Docs/HANDOVER_2026-10-23.md'
    if old.exists():
        body = old.read_text()
        if not body.startswith('> Power outage update'):
            old.write_text('> Power outage update — 2026-10-23 10:50 KST: native job 435137 accepted steps 6–7; latest t=0.00019779036422122955. Earlier pending/step-5 reports below are historical. See [current handover](POWER_OUTAGE_HANDOVER_2026-10-23.md) and the external publication receipt. Work is halted pending user resumption.\n\n' + body)


def publish():
    git('add', '--', 'Exam', 'Docs/HANDOVER_2026-10-23.md',
        'Docs/POWER_OUTAGE_HANDOVER_2026-10-23.md', 'Research/Kepler')
    paths = git('diff', '--cached', '--name-only').splitlines()
    excluded = {'.bin', '.o', '.a', '.so', '.exe', '.mod'}
    assert paths and not [p for p in paths if Path(p).suffix in excluded], paths
    message = CONTROL / 'commit_message.txt'
    message.write_text('Preserve Kepler work before scheduled power outage\n\n'
                       'Stop project work at user-requested 10:50 KST.\n'
                       'Record actual native accepted steps 6 and 7, cost gate failure,\n'
                       'latest geometry code, exact source evidence and restart locations.\n'
                       'Preserve large native data on shared GPFS.\n')
    git('commit', '--quiet', '--file', str(message), timeout=180)
    branch = CFG['remote_branch']
    fallback = None
    try:
        git('push', 'origin', 'HEAD:refs/heads/' + branch, timeout=120)
    except Exception as first_error:
        # Never force push or discard an unpublished commit. Preserve a remote
        # backup if the primary branch changed concurrently or rejects the update.
        fallback = str(first_error)
        branch = 'codex/poweroutage-20261023-' + datetime.now().strftime('%H%M%S')
        git('push', 'origin', 'HEAD:refs/heads/' + branch, timeout=120)
    head = git('rev-parse', 'HEAD')
    remote_head = git('ls-remote', 'origin', 'refs/heads/' + branch).split()[0]
    assert head == remote_head
    receipt = dict(completed_at=now(), commit=head, branch=branch,
                   remote_hash=remote_head, remote_verified=True,
                   intended_branch=CFG['remote_branch'], primary_push_error=fallback,
                   staged_files=len(paths), remaining_git_status=git('status', '--porcelain'))
    save('publish_receipt.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    if args.plan_only:
        report = preflight()
        save('preflight.json', report)
        print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
        return
    assert datetime.now().astimezone() >= datetime.fromisoformat(CFG['deadline']), 'Stop deadline has not arrived'
    with (CONTROL / 'execution.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('Another stop invocation holds the lock.', flush=True)
            return
        if (CONTROL / 'result.json').exists() and json.loads((CONTROL / 'result.json').read_text()).get('complete'):
            return
        marker = dict(CFG, halt_active=True, stop_requested_at=now(),
                      stop_reason='user_requested_pre_poweroutage')
        marker.pop('ssh_auth_sock', None)
        (ROOT / 'POWER_OUTAGE_STOP_20261023.json').write_text(json.dumps(marker, indent=2) + '\n')
        result = dict(started_at=now(), complete=False)
        try:
            outcomes = stop_work()
            result['stop'] = outcomes
            # Update only a clean, verified publication worktree, before copying.
            assert str(PUB) == git('rev-parse', '--show-toplevel')
            assert git('remote', 'get-url', 'origin') == CFG['remote']
            assert not git('status', '--porcelain'), 'Publication worktree is not clean; preserving it for review'
            git('fetch', 'origin', CFG['remote_branch'])
            git('merge', '--ff-only', 'origin/' + CFG['remote_branch'])
            handover(outcomes)
            from snapshot import make_snapshot
            result['snapshot'] = make_snapshot()
            save('result_in_progress.json', result)
            result['publication'] = publish()
            result['complete'] = bool(outcomes.get('local', {}).get('all_project_workers_stopped')
                                      and outcomes.get('slurm', {}).get('all_project_jobs_stopped'))
        except Exception:
            result['error'] = traceback.format_exc()
            print(result['error'], flush=True)
        result['finished_at'] = now()
        save('result.json', result)
        marker['shutdown_procedure_complete'] = result['complete']
        marker['result_path'] = str(CONTROL / 'result.json')
        (ROOT / 'POWER_OUTAGE_STOP_20261023.json').write_text(json.dumps(marker, indent=2) + '\n')
        print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
        if not result['complete']:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
