from pathlib import Path
import os,sys,json,subprocess,time,socket
O=Path(__file__).resolve().parent;nx=int(sys.argv[1]);P=O/f'n{nx}';large=nx==256
cores=sorted(os.sched_getaffinity(0))[:16] if large else list(range(12,28));assert set(cores)<=os.sched_getaffinity(0) and len(cores)==16
children=[]
for kind,c in [('native',cores[:8]),('gizmo',cores[8:])]:
 with (P/f'{kind}_supervisor.log').open('x') as f:
  p=subprocess.Popen([sys.executable,str(O/'worker.py'),str(nx),kind,','.join(map(str,c))],cwd=O,env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),stdout=f,stderr=f);children.append((kind,p))
(P/'launch.json').write_text(json.dumps(dict(host=socket.gethostname(),launcher_pid=os.getpid(),started_epoch=time.time(),workers={k:p.pid for k,p in children},cores=cores),indent=2));print('LAUNCHED',nx,[(k,p.pid) for k,p in children],flush=True)
codes={kind:p.wait() for kind,p in children};(P/'workers_exit.json').write_text(json.dumps(codes));print('EXIT',codes,flush=True)
sys.exit(0 if all(v==0 for v in codes.values()) else 2)
