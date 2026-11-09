from pathlib import Path
import os,sys,time,json,re,subprocess,shutil
from datetime import datetime
from zoneinfo import ZoneInfo
R=Path('/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1')
O=R/'Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/rt_t4_approved_20261109_v1'
shutil.copy2(__file__,O/'monitor.py')
shutil.copy2('/tmp/render_rt_t4_approved_20261109.py',O/'render_density.py')
(O/'monitor_process.json').write_text(json.dumps({'pid':os.getpid(),'scope':'Read existing approved RT progress; after successful t=4 only render density. No new evolution or retry.'},indent=2))
while True:
    s=json.loads((O/'status.json').read_text())
    if s['state'] in ['evolution_completed','failed']:break
    p=Path('/proc')/str(s.get('supervisor_pid',0))
    if not p.exists():raise RuntimeError('Evolution supervisor missing; no automatic restart')
    with (O/'log').open('rb') as f:
        f.seek(max(0,(O/'log').stat().st_size-262144));lines=f.read().decode(errors='replace').splitlines()
    accepted=[x for x in lines if x.startswith('[MATERIAL_SSP_ACCEPT]')]
    if accepted:
        v=accepted[-1];t=float(re.search(r' t=([^ ]+)',v)[1]);dt=float(re.search(r' dt=([^ ]+)',v)[1]);retry=int(re.search(r' retries=(\d+)',v)[1])
        record=dict(time_KST=datetime.now(ZoneInfo('Asia/Seoul')).isoformat(),t=t,target=4,dt=dt,last_step_retries=retry,native_pid=s.get('native_pid'),native_exists=(Path('/proc')/str(s.get('native_pid',0))).exists())
        with (O/'monitor_history.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
        print(json.dumps(record),flush=True)
    time.sleep(55)
if s['state']=='failed':
    print('EVOLUTION_FAILED',json.dumps(s),flush=True);raise SystemExit(2)
print('Final checkpoint completed; rendering actual density cells.',flush=True)
with (O/'render.log').open('x') as f:
    p=subprocess.Popen(['taskset','-c','18,19',sys.executable,str(O/'render_density.py')],cwd=R,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLCONFIGDIR=str(O/'mplconfig')),stdout=f,stderr=f)
    (O/'render_process.json').write_text(json.dumps({'pid':p.pid,'started_epoch':time.time()},indent=2));code=p.wait()
result=dict(state='completed_awaiting_approval' if code==0 else 'postprocessing_failed',updated_epoch=time.time(),saved_time=4.0,render_exit=code,result=str(O/'result.json'),figure=str(O/'figures/density_t4.png') if code==0 else None,next_stage='Not authorized; no further simulations or job release.')
(O/'status.json').write_text(json.dumps(result,indent=2)+'\n')
print('APPROVED_STAGE_TERMINAL',json.dumps(result),flush=True)
raise SystemExit(code)
