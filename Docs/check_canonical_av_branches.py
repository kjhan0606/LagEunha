"""One-shot process and bounded-log check; no polling loop."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys

root=Path(sys.argv[1])
manifest=json.loads((root/'manifest.json').read_text())
pids=subprocess.run(['pgrep','-x','eunha2.exe'],capture_output=True,text=True).stdout.split()
live={}
for pid in pids:
    try:
        proc=Path('/proc')/pid
        cwd=os.readlink(proc/'cwd')
        cmd=(proc/'cmdline').read_bytes().replace(b'\0',b' ').decode()
        if cmd.startswith(manifest['binary']+' '):live.setdefault(cwd,[]).append({'pid':int(pid),'cmdline':cmd})
    except FileNotFoundError:pass
out={'observed_utc':datetime.now(timezone.utc).isoformat(),'cases':{}}
for case in manifest['cases']:
    d=root/case['name'];entry={'processes':live.get(str(d),[])}
    if (d/'exit.code').exists():entry['exit_code']=int((d/'exit.code').read_text())
    log=d/'log'
    if log.exists():
        with log.open('rb') as f:
            f.seek(max(0,log.stat().st_size-262144));lines=f.read().decode(errors='replace').splitlines()
        accepted=[s for s in lines if s.startswith('[ACCEPTED_STATE]')]
        if accepted:
            entry['accepted']=dict(re.findall(r'(\w+)=([^\s]+)',accepted[-1]))
            entry['absolute_time']=manifest['start_time']+float(entry['accepted']['elapsed'])
        else:entry['last_lines']=lines[-4:]
    out['cases'][case['name']]=entry
print(json.dumps(out,indent=2))
