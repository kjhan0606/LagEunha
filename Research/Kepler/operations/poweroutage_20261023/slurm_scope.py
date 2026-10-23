#!/usr/bin/env python3
"""Inventory/stop only the user's jobs bound to the Kepler project paths."""
import json, os, shlex, subprocess, sys, time
from pathlib import Path

CONTROL = Path(__file__).resolve().parent
CFG = json.loads((CONTROL / 'config.json').read_text())


def call(args, check=True):
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=45)
    if check and p.returncode:
        raise RuntimeError(f'{args[0]} failed: {p.stderr}')
    return p


def inventory():
    assert os.getuid() == CFG['owner_uid']
    p = call(['squeue', '-h', '-u', CFG['owner_name'],
              '-o', '%i|%u|%T|%Z|%o'])
    found = []
    for line in p.stdout.splitlines():
        job, owner, state, cwd, command = line.split('|', 4)
        if owner != CFG['owner_name']:
            continue
        words = shlex.split(command)
        prefix = CFG['job_path_prefix']
        if cwd.startswith(prefix) or any(w.startswith(prefix) for w in words):
            found.append(dict(job_id=job, owner=owner, state=state,
                              cwd=cwd, command=command))
    return found


def stop():
    before = inventory()
    events = []
    for j in before:
        args = ['scancel']
        if j['state'] != 'PENDING':
            args += ['--signal=TERM', '--full']
        args += [j['job_id']]
        p = call(args, check=False)
        events.append(dict(args=args, exit=p.returncode, stderr=p.stderr))
    until = time.monotonic() + (25 if before else 0)
    while time.monotonic() < until and inventory():
        time.sleep(2)
    for j in inventory():
        p = call(['scancel', '--signal=KILL', '--full', j['job_id']], False)
        events.append(dict(job_id=j['job_id'], signal='KILL', exit=p.returncode))
        call(['scancel', j['job_id']], False)
    until = time.monotonic() + (15 if before else 0)
    after = inventory()
    while after and time.monotonic() < until:
        time.sleep(2)
        after = inventory()
    return dict(before=before, events=events, after=after,
                all_project_jobs_stopped=not after)


if __name__ == '__main__':
    action = sys.argv[1]
    if action == 'inventory':
        result = dict(jobs=inventory())
    elif action == 'stop':
        from datetime import datetime
        assert datetime.now().astimezone() >= datetime.fromisoformat(CFG['deadline'])
        result = stop()
    else:
        raise ValueError(action)
    print(json.dumps(result, ensure_ascii=False))
