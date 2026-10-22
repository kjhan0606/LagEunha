from pathlib import Path
import json, subprocess, time
if Path('log').exists():
    raise SystemExit('Refusing to overwrite log')
t0=time.monotonic()
with open('log','x') as out:
    p=subprocess.run(['./eunha2.exe','params.dat'],stdout=out,stderr=subprocess.STDOUT)
    out.write(f'EXIT:{p.returncode}\n')
Path('timing.json').write_text(json.dumps({'elapsed_seconds':time.monotonic()-t0,'ranks':1,'exit_code':p.returncode,'shared_allocation':True},indent=2)+'\n')
raise SystemExit(p.returncode)
