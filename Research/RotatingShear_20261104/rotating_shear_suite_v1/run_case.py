from pathlib import Path
import subprocess,os,json,time,sys,struct,socket
from datetime import datetime,timezone
d=Path(sys.argv[1]).resolve();cfg=json.loads((d/'registration.json').read_text());extra=json.loads((d/'environment.json').read_text())
assert not (d/'log').exists() and not (d/'result.json').exists()
base={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','HYDRO_','SEDOV_'))};start=time.monotonic();utc=datetime.now(timezone.utc).isoformat()
with (d/'log').open('x') as f:
 p=subprocess.Popen(cfg['command'],cwd=d,env={**base,**extra},stdout=f,stderr=subprocess.STDOUT)
 (d/'process.json').write_text(json.dumps({'pid':p.pid,'node':socket.gethostname(),'start_utc':utc,'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'command':cfg['command']},indent=2)+'\n');rc=p.wait()
outputs=[]
for path in sorted(d.glob('khout.*.dat')):
 if path.is_symlink():continue
 with path.open('rb') as f:f.seek(-52,2);t=struct.unpack('<d',f.read(8))[0]
 outputs.append({'path':str(path),'time':t,'bytes':path.stat().st_size})
latest=max([r['time'] for r in outputs],default=None)
row={'case':cfg['name'],'exit':rc,'wall_seconds':time.monotonic()-start,'latest_time':latest,'target':cfg['target'],'reached_target':latest==cfg['target'],'outputs':outputs,'node':socket.gethostname(),'slurm_job_id':os.environ.get('SLURM_JOB_ID')}
(d/'result.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
sys.exit(0 if rc==0 and latest==cfg['target'] else 1)
