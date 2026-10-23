#!/usr/bin/env python3
"""Fallback for user-manager logout; the stop script enforces single execution."""
import json
import os
from pathlib import Path
import subprocess
import time
from datetime import datetime

control = Path(__file__).resolve().parent
cfg = json.loads((control / 'config.json').read_text())
deadline = datetime.fromisoformat(cfg['deadline']).timestamp() + 5
(control / 'watchdog_pid.json').write_text(json.dumps({
    'pid': os.getpid(), 'deadline': cfg['deadline'],
    'purpose': 'Fallback only, same execution lock as systemd timer'}) + '\n')
while time.time() < deadline:
    time.sleep(min(20, deadline - time.time()))
if not (control / 'result.json').exists():
    with (control / 'watchdog.log').open('a') as out:
        subprocess.run(['/usr/bin/python3', str(control / 'stop_commit_push.py')],
                       stdout=out, stderr=subprocess.STDOUT, check=False)
