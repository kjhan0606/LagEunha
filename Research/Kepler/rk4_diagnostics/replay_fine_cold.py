"""At most one full, two half and four quarter native attempt steps."""
from pathlib import Path
import csv,json,os,re,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
base=root/'rk4_diagnostics';dest=base/'fine_cold_interval';dest.mkdir(exist_ok=False)
source=root/'core_refinement/native_v4/fine_old_raw'
starts=list(source.glob('khdiag.step*.start.dat'));assert len(starts)==1
start=starts[0];failed=Path(str(start).replace('.start.dat','.failed.dat'))
p0,m=read_snapshot(start);_,mf=read_snapshot(failed);h=mf['dt'];target=m['t']+h
binary=root/'raw_vertex_dedup/v5/eunha2_raw_vertex_dedup.exe'
focus={68526,68592};selected=focus.copy()
initial,_=read_snapshot(root/'core_refinement/fine_integrated/khout.000001.dat')
initial_ids=(initial['u4if']&np.uint64((1<<56)-1)).astype(int)
r0=np.hypot(initial['x']-12,initial['y']-12)
selected.update(initial_ids[abs(r0-.018528970665049098)<1e-12].tolist())
manifest={'source':str(start),'failed':str(failed),'start_time':m['t'],'interval':h,'target_time':target,
          'scope':'bounded native temporal diagnosis; no timestep retry implementation or physical fix adopted',
          'cases':[]}
for subdivisions in (1,2,4):
    label={1:'full',2:'halves',4:'quarters'}[subdivisions]
    d=dest/label;d.mkdir()
    (d/'khout.000001.dat').symlink_to(start);(d/'params.dat').symlink_to(source/'params.dat')
    (d/'eunha2.exe').symlink_to(binary)
    # Allow only rounding-scale differences between the nominal equal pieces.
    planned=np.array([m['t']+k*h/subdivisions for k in range(subdivisions+1)])
    cap=float(np.diff(planned).max())
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        KEPLER_RESTART_STEP='1',EUNHA_IC='kepler',HYDRO_TSTOP=repr(target),EUNHA_DUMP_DT='1.e-6',
        GFS_EXACT_STOP='1',GFS_DIAG_DT_MAX=repr(cap),GFS_FAIL_NONPOSITIVE_U='1',GFS_RK4_FAILURE_SNAPSHOT='1',
        GFS_RK4_AUDIT_STEP_MIN='1',GFS_RK4_AUDIT_STEP_MAX='10',GFS_GPU_TESS_DUMP='0',
        GFS_GPU_FACE_DUMP_CALL='1' if subdivisions==1 else '0',GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',GFS_VELOCITY_LIMITER='0')
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT','GFS_RK4_AUDIT_STAGE','GFS_RK4_AUDIT_IDS'):env.pop(key,None)
    if subdivisions==1:env['GFS_RK4_AUDIT_STAGE']='k1_rhs'
    else:env['GFS_RK4_AUDIT_IDS']=','.join(str(i) for i in sorted(selected))
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','KEPLER_','EUNHA_'))},indent=2)+'\n')
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    row={'name':label,'subdivisions':subdivisions,'dt_ceiling':cap,'exit_code':code}
    log=(d/'log').read_text()
    row['invalid_stages']=re.findall(r'\[RK4_INVALID_STATE\].*',log)
    diagnostics=list(d.glob('khdiag.step*.failed.dat'))
    if diagnostics:
        p,meta=read_snapshot(diagnostics[0]);row['failure_snapshot']=str(diagnostics[0]);row['minimum_failed_U']=float(p['ie'].min())
        row['attempted_steps']=int(re.search(r'step=(\d+)',row['invalid_stages'][0])[1])
        if subdivisions==1:
            old,_=read_snapshot(failed);oa=np.argsort(old['u4if']);ob=np.argsort(p['u4if'])
            row['replay_max_absolute_error']={k:float(abs(old[k][oa].astype(float)-p[k][ob].astype(float)).max()) for k in ('x','y','vx','vy','volume','mass','ie')}
    else:
        outputs=[f for f in d.glob('khout.*.dat') if f.name!='khout.000001.dat']
        row['attempted_steps']=len(outputs)
        if outputs:
            last=max(outputs,key=lambda p:int(p.name.split('.')[1]));p,meta=read_snapshot(last)
            row['time']=meta['t'];row['same_endpoint']=meta['t']==target
    if row['attempted_steps']>subdivisions:raise RuntimeError('diagnostic step budget exceeded; do not continue')
    if subdivisions==1:
        f=np.loadtxt(d/'face_dump_r0_call1.tsv',skiprows=1)
        selected.update(f[np.isin(f[:,4].astype(int),list(focus)),5].astype(int).tolist())
        if len(selected)>64:raise RuntimeError('target/direct-neighbor audit selection exceeds explicit64-ID capacity')
        row['selected_ids_for_subdivision_stages']=sorted(selected)
    manifest['cases'].append(row)
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(label,'exit',code,'attempts',row['attempted_steps'],'negative',row['invalid_stages'][:2],flush=True)
    if code==0:break
    if code!=92:raise RuntimeError(f'unexpected diagnostic exit {code}; do not continue')
