"""One captured interval only: unchanged full step versus capped half steps."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
args=argparse.ArgumentParser();args.add_argument('--run-dir',default='failure_replay');args.add_argument('--binary',default='rk4_diagnostics/v3/eunha2_rk4_diagnostics.exe');args=args.parse_args()
base=root/'rk4_diagnostics';dest=base/args.run_dir;dest.mkdir(exist_ok=False)
capture=json.loads((base/'validation/manifest.json').read_text())['sedov_failure']
start=Path(capture['valid_step_start_snapshot']);failed=Path(str(start).replace('.start.dat','.failed.dat'))
_,m=read_snapshot(start);_,mf=read_snapshot(failed)
h=mf['dt'];target=m['t']+h
binary=root/args.binary
cases=[]
for label,strict,cap in (('full_strict',1,h),('full_observe',0,h),('halves_strict',1,h/2)):
    d=dest/label;d.mkdir();(d/'khout.000001.dat').symlink_to(start)
    (d/'params.dat').symlink_to(start.parent/'params.dat');(d/'eunha2.exe').symlink_to(binary)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        KEPLER_RESTART_STEP='1',EUNHA_IC='shock_probe',HYDRO_TSTOP=repr(target),EUNHA_DUMP_DT='1.e-6',
        GFS_EXACT_STOP='1',GFS_DIAG_DT_MAX=repr(cap),GFS_FAIL_NONPOSITIVE_U=str(strict),GFS_RK4_FAILURE_SNAPSHOT=str(strict),
        GFS_RK4_AUDIT_STEP_MIN='1',GFS_RK4_AUDIT_STEP_MAX='100',GFS_GPU_FACE_DUMP_CALL='0',
        GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',GFS_VELOCITY_LIMITER='0')
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT','GFS_RK4_AUDIT_STAGE'):
        env.pop(key,None)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','EUNHA_','HYDRO_','KEPLER_'))},indent=2)+'\n')
    cases.append((d,env))
def execute(item):
    d,env=item
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    return d,code
result={'source_step_start':str(start),'original_failed_state':str(failed),'start_time':m['t'],
        'original_step':h,'target_time':target,
        'scope':'One common interval. Actual RK4 evolves internal Udot; conservative-total-energy FE bound is diagnostic only.',
        'cases':{}}
with ThreadPoolExecutor(max_workers=2) as pool:
    for future in as_completed([pool.submit(execute,c) for c in cases]):
        d,code=future.result();row={'exit_code':code,'stages':[],'focus_cells':{}}
        rows=np.genfromtxt(d/'rk4_stage_r0.tsv',names=True,dtype=None,encoding=None,delimiter='\t')
        keys=list(dict.fromkeys(zip(rows['step'],rows['stage'])))
        for step,stage in keys:
            z=rows[(rows['step']==step)&(rows['stage']==stage)]
            valid=z['rhs_valid']>0
            row['stages'].append({'step':int(step),'stage':str(stage),'dt':float(z['h'][0]),'min_U':float(z['U_raw'].min()),
                                  'negative_cells':int(np.sum(z['U_raw']<=0)),
                                  'linear_positive_dt_min':float(z['linear_U_positive_dt'][valid].min()) if valid.any() else None,
                                  'conserved_FE_positive_dt_min':float(z['conserved_FE_positive_dt'][valid].min()) if valid.any() else None})
        for i in (1699,1829):
            z=rows[rows['id']==i]
            row['focus_cells'][str(i)]=[{field:(str(v[field]) if field=='stage' else float(v[field])) for field in ('step','stage','h','U_raw','V','rho_stored','P_stored','Udot','ax_hydro','ay_hydro','die_stored','linear_U_positive_dt','conserved_FE_positive_dt')} for v in z]
        outputs=[p for p in d.glob('khout.*.dat') if p.name!='khout.000001.dat']
        if outputs:
            last=max(outputs,key=lambda p:int(p.name.split('.')[1]));p,meta=read_snapshot(last)
            row['final_time']=meta['t'];row['final_min_U']=float(p['ie'].min());row['native_steps_completed']=len(outputs)
            row['same_endpoint']=meta['t']==target
        diag=sorted(d.glob('khdiag.step*.failed.dat'))
        if diag:
            p,meta=read_snapshot(diag[0]);row['failure_snapshot']=str(diag[0]);row['failure_min_U']=float(p['ie'].min())
            if d.name=='full_strict':
                a,_=read_snapshot(failed)
                oa=np.argsort(a['u4if']);op=np.argsort(p['u4if'])
                errors={field:float(np.max(abs(a[field][oa].astype(float)-p[field][op].astype(float)))) for field in ('x','y','vx','vy','ie','volume','mass')}
                row['replay_vs_captured_failed_state_max_absolute']=errors
        result['cases'][d.name]=row
        (dest/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(d.name,'exit',code,'negative_stages',[(z['step'],z['stage'],z['negative_cells']) for z in row['stages'] if z['negative_cells']],flush=True)
