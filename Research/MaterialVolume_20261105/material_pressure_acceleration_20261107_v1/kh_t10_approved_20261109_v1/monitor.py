"""Finite monitor for the approved stage only. No submissions or solver changes."""
from pathlib import Path
import os,sys,json,time,re,subprocess,traceback,h5py
from gizmo_audit import audit,final
O=Path(__file__).resolve().parent;nx=int(sys.argv[1]);P=O/f'n{nx}'
def write(p,x):
 q=p.with_suffix('.tmp');q.write_text(json.dumps(x,indent=2)+'\n');q.replace(p)
def tail(p):
 if not p.exists():return ''
 with p.open('rb') as f:f.seek(max(0,p.stat().st_size-256000));return f.read().decode(errors='replace')
try:
 while True:
  status=dict(pid=os.getpid(),epoch=time.time(),target=10.,workers={})
  for kind in ['native','gizmo']:
   sp=P/kind/'status.json'
   if sp.exists():
    s=json.loads(sp.read_text());log=Path(s.get('directory',str(sp.parent)))/'log';text=tail(log)
    pat=r'\[MATERIAL_SSP_ACCEPT\] t=([^ ]+)' if kind=='native' else r'Sync-Point \d+, Time: ([\d.eE+-]+)';times=re.findall(pat,text)
    s['latest_observed_time']=float(times[-1]) if times else s.get('saved_time',s.get('start'));s['log_bytes']=log.stat().st_size if log.exists() else 0;status['workers'][kind]=s
  available={}
  for p in (P/'gizmo/output').glob('snap_*.hdf5'):
   try:
    with h5py.File(p) as f:available[float(f['Header'].attrs['Time'])]=str(p)
   except OSError:pass
   else:audit(nx,p)
  if status['workers'].get('gizmo',{}).get('state')=='completed':final(nx)
  done=[]
  for target in [2.1,4.7,9.2,10.]:
   out=P/'comparisons'/f't{target:g}';report=out/'comparison.json'
   if report.exists():done.append(target);continue
   if (P/'native'/f'endpoint_{target:g}.json').exists() and any(abs(t-target)<1e-10 for t in available):
    (P/'comparisons').mkdir(exist_ok=True)
    with (P/'comparisons'/f't{target:g}.log').open('x') as f:code=subprocess.run([sys.executable,str(O/'compare.py'),str(nx),str(target)],env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),stdout=f,stderr=f).returncode
    if code:raise RuntimeError(f'comparison failed at{target}; evolution is separate')
    done.append(target)
  status['compared_times']=done;states=[s['state'] for s in status['workers'].values()]
  status['state']='completed_awaiting_approval' if len(done)==4 and states==['completed','completed'] else ('waiting_for_slurm' if not states else 'running')
  if len(states)==2 and all(s in ['completed','failed'] for s in states) and 'failed' in states:status['state']='evolution_failed_awaiting_review'
  write(P/'progress.json',status)
  if status['state'] not in ['running','waiting_for_slurm']:break
  time.sleep(60)
except Exception as exc:
 write(P/'monitor_error.json',dict(error=repr(exc),pid=os.getpid()));traceback.print_exc();raise
