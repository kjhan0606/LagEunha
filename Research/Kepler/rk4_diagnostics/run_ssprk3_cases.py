"""Continue only outstanding native interval cases; preserve prior proofs."""
from pathlib import Path
import argparse,csv,hashlib,json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['intervals','fine_active_half','fresh'],default='intervals')
args=parser.parse_args()
dest=root/'rk4_diagnostics'/dict(intervals='ssprk3_intervals_v2',fine_active_half='ssprk3_fine_active_half_v2',fresh='ssprk3_fresh_v2')[args.phase]
dest.mkdir(exist_ok=False)
binary=root/'rk4_diagnostics/ssprk3_v2_build/eunha2_ssprk3.exe'
assert (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
manifest={'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'phase':args.phase,'cases':[],
    'existing_proofs':str(root/'rk4_diagnostics/ssprk3_v1/manifest.json'),
    'scope':'SSPRK3 plus current-stage EOS/gradient preparation; not an isolated same-RHS temporal comparison or long stability proof',
    'output_policy':'one final native checkpoint per successful case; explicit exact-endpoint output. About0.1GiB for four interval cases, about39MiB fresh.'}
def save(): (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def energy(p,kepler):
    e=p['ie'].astype(float)+.5*p['mass']*(p['vx']**2+p['vy']**2)
    if kepler:e-=p['mass']/np.sqrt((p['x']-12)**2+(p['y']-12)**2+.047**2)
    return float(e.sum())

def execute(name,initial,params,stop,cap,focus,max_rhs=256,max_steps=64):
    d=dest/name;d.mkdir();(d/'OWNER.txt').write_text('LagForce Codex native SSPRK3 bounded verification. Preserve outputs.\n')
    for key,path in [('khout.000001.dat',initial),('params.dat',params),('eunha2.exe',binary)]: (d/key).symlink_to(path)
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('GFS_'):env.pop(key)
    kepler=name.startswith('fine')
    env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        EUNHA_IC='kepler' if kepler else 'shock_probe',KEPLER_RESTART_STEP='1',
        HYDRO_TSTOP=repr(stop),EUNHA_DUMP_DT=repr(stop),GFS_EXACT_STOP='1',GFS_DIAG_DT_MAX=repr(cap),
        GFS_SSPRK3='1',GFS_SSP_MAX_RHS=str(max_rhs),GFS_SSP_MAX_STEPS=str(max_steps),
        GFS_FAIL_NONPOSITIVE_U='1',GFS_GPU_FINAL_VOLUME='1',GFS_GPU_RAW_FACES='1',
        GFS_VELOCITY_LIMITER='0',GFS_KEPLER_HALO_CELLS='8',GFS_E_AUDIT='1',GFS_FLOOR_LOG='1',
        GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_TESS_DUMP='0',GFS_RK4_AUDIT_STEP_MIN='1',
        GFS_RK4_AUDIT_STEP_MAX=str(max_steps+1),GFS_RK4_AUDIT_IDS=','.join(map(str,focus)))
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','EUNHA_','KEPLER_','NATIVE_'))},indent=2)+'\n')
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    log=(d/'log').read_text();files=sorted(f for f in d.glob('khout.*.dat') if f.name!='khout.000001.dat')
    row=dict(name=name,exit_code=code,stop=stop,cap=cap,initial=str(initial),effective_params=str(d/'params.dat'),
             rollbacks=log.count('[SSP_ROLLBACK]'),accepted=log.count('[SSP_ACCEPT]'),
             elapsed_seconds=json.loads((d/'timing.json').read_text())['elapsed_seconds'])
    steps=np.genfromtxt(d/'ssprk3_steps.tsv',names=True,delimiter='\t',ndmin=1) if (d/'ssprk3_steps.tsv').exists() else None
    if files:
        p,m=read_snapshot(files[-1]);q,_=read_snapshot(initial);ia=np.argsort(p['u4if']);ib=np.argsort(q['u4if'])
        row.update(final=str(files[-1]),time=m['t'],same_endpoint=m['t']==stop,min_U=float(p['ie'].min()),
            individual_mass_exact=np.array_equal(p['mass'][ia],q['mass'][ib]),W0_exact=bool(np.all(p['w2']==0)),
            relative_total_energy_change=abs(energy(p,kepler)-energy(q,kepler))/max(abs(energy(q,kepler)),1e-300))
        if kepler:
            ic,_=read_snapshot(root/'core_refinement/fine_integrated/khout.000001.dat')
            icorder=np.argsort(ic['u4if']);gamma=float(np.float32(1.666667))
            Uad=ic['ie'][icorder].astype(float)*(ic['volume'][icorder].astype(float)/p['volume'][ia])**(gamma-1)
            ids=(p['u4if'][ia]&np.uint64((1<<56)-1)).astype(int)
            row['focus_entropy_ratio']={str(ident):float((p['ie'][ia]/Uad)[ids==ident][0]) for ident in focus}
    else:row['missing_checkpoint']=True
    if steps is not None and len(steps):
        row.update(rhs_total=int(steps['rhs_total'][-1]),floor_injection=float(steps['floor_injection'].sum()),
            net_hydro_work=float(steps['hydro_integral'].sum()),sum_hydro_abs_integral=float(steps['hydro_abs_integral'].sum()),
            finite_step_hydro_residual=float(steps['finite_step_hydro_residual'].sum()),
            finite_step_potential_residual=float(steps['finite_step_potential_residual'].sum()))
    manifest['cases'].append(row);save();print(json.dumps(row),flush=True)
    if code or not row.get('same_endpoint') or not row.get('individual_mass_exact') or row.get('min_U',-1)<=0 or row.get('floor_injection',-1)!=0:
        manifest['classification']='STOP: failed native interval, admissibility or output gate';save();raise SystemExit(2)
    return row

if args.phase=='intervals':
    prior=json.loads((root/'rk4_diagnostics/ssprk3_v1/manifest.json').read_text())
    for key in ('default_comparison','rollback_mpi1','rollback_mpi4','MPI1_MPI4'):
        assert prior[key]['all_recorded_fields_bitwise_equal'] and prior[key]['metadata_equal']
    sources=[('fine',Path(json.loads((root/'rk4_diagnostics/fine_cold_interval/manifest.json').read_text())['source']),[67935,68001,68526,68592]),
             ('sedov',root/'rk4_diagnostics/validation/sedov_first_failure/khdiag.step000256.start.dat',[1699,1829])]
    for name,initial,focus in sources:
        _,m=read_snapshot(initial);_,mf=read_snapshot(str(initial).replace('.start.dat','.failed.dat'))
        h=mf['dt'];stop=m['t']+h;rows=[]
        for factor in (1.,.5):
            rows.append(execute(name+('_base' if factor==1 else '_half'),initial,initial.parent/'params.dat',stop,h*factor,focus))
        p,mp=read_snapshot(rows[0]['final']);q,mq=read_snapshot(rows[1]['final'])
        a=np.argsort(p['u4if']);b=np.argsort(q['u4if'])
        diff={k:float(abs(p[k][a].astype(float)-q[k][b].astype(float)).max()) for k in ('x','y','vx','vy','volume','ie')}
        manifest[name+'_ceiling_comparison']=dict(max_absolute_errors=diff,endpoint_time=mp['t'],same_endpoint=mp['t']==mq['t'])
        save()
elif args.phase=='fine_active_half':
    prior=root/'rk4_diagnostics/ssprk3_intervals_v2/fine_base'
    steps=np.genfromtxt(prior/'ssprk3_steps.tsv',names=True,delimiter='\t',ndmin=1)
    cap=.5*float(steps['h'].max())
    old=json.loads((root/'rk4_diagnostics/fine_cold_interval/manifest.json').read_text())
    initial=Path(old['source'])
    execute('fine_active_half',initial,initial.parent/'params.dat',old['target_time'],cap,[67935,68001,68526,68592])
    manifest['reason']='Nominal half-ceiling did not bind in prior comparison; this is the one effective temporal refinement comparison, not a scan.'
else:
    prior=json.loads((root/'rk4_diagnostics/ssprk3_intervals_v2/manifest.json').read_text())
    assert len(prior['cases'])==4 and all(r['exit_code']==0 and r['same_endpoint'] for r in prior['cases'])
    fine=root/'core_refinement/fine_integrated'
    execute('fine_fresh',fine/'khout.000001.dat',fine/'params.dat',.2,1.,[67935,68001,68526,68592],3072,768)
manifest['classification']='completed; inspect local accuracy, entropy and energy before any long continuation'
save();print(manifest['classification'],flush=True)
