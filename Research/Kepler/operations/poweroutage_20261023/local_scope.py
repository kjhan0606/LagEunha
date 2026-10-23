#!/usr/bin/env python3
"""Stop owned project workers, preserving editors, agents and unrelated jobs."""
import json, os, signal, time
from pathlib import Path

CONTROL = Path(__file__).resolve().parent
CFG = json.loads((CONTROL / 'config.json').read_text())


def process(pid):
    p = Path('/proc') / str(pid)
    try:
        if p.stat().st_uid != CFG['owner_uid']:
            return None
        fields = (p / 'stat').read_text().rsplit(')', 1)[1].split()
        argv = (p / 'cmdline').read_bytes().decode(errors='replace').split('\0')
        argv = [a for a in argv if a]
        if not argv or fields[0] == 'Z':
            return None
        try:
            cwd = os.readlink(p / 'cwd')
        except OSError:
            cwd = ''
        return dict(pid=pid, ppid=int(fields[1]), start=fields[19],
                    argv=argv, cwd=cwd)
    except (OSError, ValueError, IndexError):
        return None


def inventory():
    assert os.getuid() == CFG['owner_uid']
    processes = {}
    for p in Path('/proc').iterdir():
        if p.name.isdigit():
            row = process(int(p.name))
            if row:
                processes[row['pid']] = row
    protected = {os.getpid()}
    pid = os.getppid()
    while pid in processes and pid not in protected:
        protected.add(pid)
        pid = processes[pid]['ppid']
    # Protect all descendants of this scheduled control invocation as well.
    changed = True
    own_descendants = {os.getpid()}
    while changed:
        extra = {pid for pid, p in processes.items()
                 if p['ppid'] in own_descendants}
        changed = not extra.issubset(own_descendants)
        own_descendants |= extra
    protected |= own_descendants
    selected = set()
    for pid, p in processes.items():
        if pid in protected:
            continue
        exe = Path(p['argv'][0]).name
        worker = (exe in {'bash', 'sh', 'python', 'python3', 'mpirun',
                         'mpiexec', 'timeout', 'srun', 'g++', 'gcc', 'nvcc',
                         'continuation', 'readback', 'inspect_cell_condensation'}
                  or exe.startswith(('probe', 'hydra_', 'eunha2',
                                     'reconcile_', 'build_cell_')))
        bound = (p['cwd'].startswith(CFG['job_path_prefix'])
                 or any(a.startswith(CFG['job_path_prefix']) for a in p['argv']))
        if worker and bound:
            selected.add(pid)
    # Children of a verified project controller belong to that invocation.
    changed = True
    while changed:
        extra = {pid for pid, p in processes.items()
                 if p['ppid'] in selected and pid not in protected}
        changed = not extra.issubset(selected)
        selected |= extra
    return [processes[p] for p in sorted(selected)]


def send(row, sig):
    current = process(row['pid'])
    if current and current['start'] == row['start']:
        try:
            os.kill(row['pid'], sig)
            return True
        except ProcessLookupError:
            pass
    return False


def stop():
    before = inventory()
    selected = {p['pid'] for p in before}
    # Controller parents first so that a finished chunk cannot launch another.
    ordered = sorted(before, key=lambda p: p['ppid'] in selected)
    events = [dict(pid=p['pid'], signal='TERM', sent=send(p, signal.SIGTERM))
              for p in ordered]
    until = time.monotonic() + (20 if before else 0)
    while time.monotonic() < until and inventory():
        time.sleep(2)
    after_term = inventory()
    events += [dict(pid=p['pid'], signal='KILL', sent=send(p, signal.SIGKILL))
               for p in after_term]
    if after_term:
        time.sleep(2)
    after = inventory()
    return dict(before=before, events=events, after=after,
                all_project_workers_stopped=not after)
