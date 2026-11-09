from pathlib import Path
import os, sys, json, hashlib, shutil, subprocess, time, re, struct, traceback
import numpy as np
R=Path('/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1')
D=R/'Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1'
H=D/'common_rt128_halo_capture_20261109_v2/recovery_continuation_v1'
O=D/'rt_t4_approved_20261109_v1'
S=D/'cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1'
sys.path.insert(0,str(S));import run_gates as gates
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def when(p):
    with p.open('rb') as f:f.seek(-52,2);return struct.unpack('d',f.read(8))[0]
def load(p):
    a=np.fromfile(p).reshape(-1,13);return a[a[:,12].argsort()]
def quantities(a):
    a=a.astype(np.longdouble);m=a[:,6]
    k=np.sum(.5*m*(a[:,2]**2+a[:,3]**2));u=np.sum(a[:,7]);g=np.sum(.5*m*a[:,1])
    return dict(mass=np.sum(m),momentum_x=np.sum(m*a[:,2]),momentum_y=np.sum(m*a[:,3]),kinetic=k,internal=u,gravitational=g,total_energy=k+u+g,total_area=np.sum(a[:,8]))
def update(state,**kw):
    f=O/'status.tmp';f.write_text(json.dumps(dict(state=state,supervisor_pid=os.getpid(),updated_epoch=time.time(),**kw),indent=2));f.replace(O/'status.json')
O.mkdir(exist_ok=False);shutil.copy2(__file__,O/'run.py')
try:
    src=H/'segment_02/rtout.004823.dat';binary=D/'common_portable_20261109_v1/native/native.exe'
    assert sha(src)=='c5f5c72a42bb2e3f22dcb328d784c5c7eb389672a4d271a89664f2756bc6aa18'
    assert sha(binary)=='11a0234acc882f6f2fc405b0053c2e5210262863be665f88cc5ae91b821ddf4e'
    assert when(src)==3.6 and json.loads((H/'segment_02/result.json').read_text())['passed']
    e=json.loads((H/'segment_02/environment.json').read_text())
    old=dict(e);e.update(HYDRO_TSTOP='4.0',GFS_MATERIAL_TRACE_FAILURE_DUMP=str(O/'failed_edges.bin'))
    changed={k:{'before':old.get(k),'after':v} for k,v in e.items() if old.get(k)!=v}
    assert set(changed)=={'HYDRO_TSTOP','GFS_MATERIAL_TRACE_FAILURE_DUMP'}
    shutil.copy2(src,O/'rtout.000001.dat');shutil.copy2(H/'segment_02/params.dat',O/'params.dat')
    cmd=['taskset','-c','14,15',str(binary),'params.dat','1']
    plan=dict(authorization='User approved only existing RT restart to t=4 and density/conservation report; wait for approval before next stage.',start=3.6,target=4.0,cells=32768,shape=[128,256],source=str(src),source_sha256=sha(src),binary=str(binary),binary_sha256=sha(binary),effective_params=str(O/'params.dat'),params_sha256=sha(O/'params.dat'),environment_changes=changed,maximum_new_checkpoints=6,checkpoint_bytes=src.stat().st_size,maximum_checkpoint_bytes_including_copy=7*src.stat().st_size,postprocessing_allowance_bytes=100*1024**2,free_bytes=shutil.disk_usage(O).free,on_failure='Stop and report; no automatic tuning or restart.',after_completion='Await next explicit stage approval. Do not release held Slurm jobs.',integrator='Existing material SSPRK2 path in native C/CUDA executable; no rebuild.',output_rule='Existing rt.c dumps at crossed multiples of 0.1, optional step21045, and final step if not already saved; .4 window.')
    assert plan['free_bytes']>plan['maximum_checkpoint_bytes_including_copy']+plan['postprocessing_allowance_bytes']
    for n,v in [('plan.json',plan),('environment.json',e),('command.json',cmd)]: (O/n).write_text(json.dumps(v,indent=2)+'\n')
    print('OUTPUT_PLAN',json.dumps(plan),flush=True)
    tick=time.monotonic()
    with (O/'log').open('x') as log:
        p=subprocess.Popen(cmd,cwd=O,env=dict(gates.ne.baseenv,**e),stdout=log,stderr=log)
        (O/'process.json').write_text(json.dumps(dict(pid=p.pid,started_epoch=time.time()),indent=2))
        update('running',native_pid=p.pid,start=3.6,target=4.0)
        print('NATIVE_STARTED',p.pid,flush=True);code=p.wait()
    wall=time.monotonic()-tick
    lines=(O/'log').read_text().splitlines();accepted=[x for x in lines if x.startswith('[MATERIAL_SSP_ACCEPT]')]
    vol=[float(re.search(r'max_relative=([^ ]+)',x)[1]) for x in lines if x.startswith('[MATERIAL_VOLUME_STEP]')]
    snaps=sorted(O.glob('rtout.*.dat'));final=snaps[-1];end=when(final)
    result=dict(exit=code,start=3.6,target=4.0,saved_time=end,steps=len(accepted),wall_seconds=wall,checkpoint=str(final),checkpoint_sha256=sha(final),maximum_volume_error=max(vol,default=None),new_checkpoint_count=len(snaps)-1,checkpoint_bytes=sum(f.stat().st_size for f in snaps),last_accepted=accepted[-1] if accepted else None)
    result['passed']=False
    if code==0 and abs(end-4)<1e-12:
        extract=R/'Research/MaterialVolume_20261105/material_kh_stop_20261106_v1/extract_ids'
        subprocess.run([str(extract),str(final),str(O/'state_ids.bin')],check=True,stdout=subprocess.DEVNULL)
        a=load(O/'state_ids.bin');a0=load(D/'common_rt128_20261109_v1/initial_ids.bin');ar=load(H/'segment_02/state_ids.bin')
        q0=quantities(a0);qr=quantities(ar);q=quantities(a);rho=a[:,6]/a[:,8]
        result.update(fixed_mass_exact=bool(np.array_equal(a[:,12],a0[:,12]) and np.array_equal(a[:,6],a0[:,6])),positive=bool(np.all(a[:,7:9]>0)),finite=bool(np.isfinite(a).all()),total_area=float(q['total_area']),gravity_energy_relative=float((q['total_energy']-q0['total_energy'])/q0['total_energy']),gravity_energy_relative_since_restart=float((q['total_energy']-qr['total_energy'])/qr['total_energy']),density_min=float(rho.min()),density_max=float(rho.max()),minimum_U=float(a[:,7].min()),minimum_V=float(a[:,8].min()),mass_relative=float((q['mass']-q0['mass'])/q0['mass']),momentum_x_change=float(q['momentum_x']-q0['momentum_x']),retries=sum(int(re.search(r'retries=(\d+)',x)[1]) for x in accepted),initial={k:float(v) for k,v in q0.items()},restart={k:float(v) for k,v in qr.items()},final={k:float(v) for k,v in q.items()},momentum_scope='Horizontal momentum is conserved by periodic x and vertical wall normals; vertical momentum is affected by gravity and wall forces and is not a conserved total.',energy_scope='U+K+sum(m*0.5*y); stationary reflecting walls, external acceleration g_y=-0.5.')
        result['passed']=bool(result['fixed_mass_exact'] and result['positive'] and result['finite'] and result['maximum_volume_error'] is not None and result['maximum_volume_error']<=1e-11 and abs(result['gravity_energy_relative'])<1e-8 and result['new_checkpoint_count']<=6)
    (O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print('EVOLUTION_RESULT',json.dumps(result),flush=True)
    if not result['passed']:
        update('failed',result=result,continuation='No automatic restart or algorithm change.');raise SystemExit(2)
    update('evolution_completed',saved_time=end,result=str(O/'result.json'),next='Density postprocessing only, within approved stage.')
except Exception as exc:
    update('failed',error=repr(exc),continuation='No automatic restart or algorithm change.');traceback.print_exc();raise
