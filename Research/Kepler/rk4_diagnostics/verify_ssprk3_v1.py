"""Native RK default, periodic MPI rollback, and bounded failure intervals."""
from pathlib import Path
import hashlib,json,os,struct,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
dest=root/'rk4_diagnostics/ssprk3_v1';dest.mkdir(exist_ok=False)
binary=root/'rk4_diagnostics/ssprk3_v1_build/eunha2_ssprk3.exe'
assert (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
manifest={'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
    'scope':'default-regression, periodic MPI rollback, two captured intervals; not a long Kepler pass',
    'output_policy':'one final native snapshot per successful branch, initial IC and failed normal-start retained; about0.3GiB expected',
    'cases':[]}

def save(): (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def execute(name,initial,params,stop,cap,ssp=True,ranks=1,force=False,audit_ids=None):
    d=dest/name;d.mkdir();(d/'OWNER.txt').write_text('LagForce Codex bounded native SSPRK3 verification. Preserve results.\n')
    (d/'khout.000001.dat').symlink_to(initial);(d/'params.dat').symlink_to(params)
    (d/'eunha2.exe').symlink_to(binary)
    env=os.environ.copy()
    for key in list(env):
        if key.startswith('GFS_'):env.pop(key)
    env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS=str(ranks),I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
        KEPLER_RESTART_STEP='1',HYDRO_TSTOP=repr(stop),EUNHA_DUMP_DT=repr(stop),
        GFS_EXACT_STOP='1',GFS_DIAG_DT_MAX=repr(cap),GFS_SSPRK3=str(int(ssp)),
        GFS_SSP_MAX_RHS='256',GFS_SSP_MAX_STEPS='64',GFS_SSP_FORCE_REJECT_ONCE=str(int(force)),
        GFS_FAIL_NONPOSITIVE_U='1',GFS_RK4_FAILURE_SNAPSHOT='0',GFS_GPU_RAW_FACES='1',
        GFS_GPU_FINAL_VOLUME='1',GFS_VELOCITY_LIMITER='0',GFS_KEPLER_HALO_CELLS='8',
        GFS_E_AUDIT='1',GFS_FLOOR_LOG='1',GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_TESS_DUMP='0')
    for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
    if name.startswith('uniform'):env['EUNHA_IC']='affine_probe'
    elif name.startswith('sedov'):env['EUNHA_IC']='shock_probe'
    else:env['EUNHA_IC']='kepler'
    if audit_ids:
        env.update(GFS_RK4_AUDIT_STEP_MIN='1',GFS_RK4_AUDIT_STEP_MAX='100000',
                   GFS_RK4_AUDIT_IDS=','.join(map(str,audit_ids)))
    (d/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','EUNHA_','HYDRO_','KEPLER_','NATIVE_'))},indent=2)+'\n')
    code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
    row={'name':name,'exit_code':code,'ranks':ranks,'ssp':ssp,'cap':cap,'stop':stop,'initial':str(initial),'effective_params':str(d/'params.dat')}
    outputs=sorted(f for f in d.glob('khout.*.dat') if f.name!='khout.000001.dat')
    log=(d/'log').read_text()
    row['rollbacks']=log.count('[SSP_ROLLBACK]');row['accepted_steps']=log.count('[SSP_ACCEPT]')
    if outputs:
        p,m=read_snapshot(outputs[-1]);row.update(final=str(outputs[-1]),time=m['t'],min_U=float(p['ie'].min()),same_endpoint=m['t']==stop)
        q,_=read_snapshot(initial);a=np.argsort(p['u4if']);b=np.argsort(q['u4if'])
        row['individual_mass_exact']=np.array_equal(p['mass'][a],q['mass'][b])
        row['W0_exact']=bool(np.all(p['w2']==0))
    if ssp and code==0 and outputs:
        steps=np.genfromtxt(d/'ssprk3_steps.tsv',names=True,delimiter='\t',ndmin=1)
        row['rhs_total']=int(steps['rhs_total'][-1]);row['floor_injection']=float(steps['floor_injection'].sum())
        row['max_relative_mass_ledger_error']=float(abs(steps['mass']-steps['mass'][0]).max()/steps['mass'][0])
        row['net_hydro_work']=float(steps['hydro_integral'].sum())
        row['sum_hydro_abs_integral']=float(steps['hydro_abs_integral'].sum())
        row['finite_step_hydro_residual']=float(steps['finite_step_hydro_residual'].sum())
        row['finite_step_potential_residual']=float(steps['finite_step_potential_residual'].sum())
        def energy(a):
            vel2=a['vx']**2+a['vy']**2; e=a['ie'].astype(float)+.5*a['mass']*vel2
            if env['EUNHA_IC']=='kepler':e-=a['mass']/np.sqrt((a['x']-12)**2+(a['y']-12)**2+.047**2)
            return float(e.sum())
        row['relative_total_energy_change']=abs(energy(p)-energy(q))/max(abs(energy(q)),1e-300)
    if code==0 and not outputs:row['checkpoint_missing']=True
    manifest['cases'].append(row);save();print(json.dumps(row),flush=True)
    return d,row

def compare(a,b):
    p,mp=read_snapshot(a);q,mq=read_snapshot(b);ia=np.argsort(p['u4if']);ib=np.argsort(q['u4if'])
    assert np.array_equal(p['u4if'][ia],q['u4if'][ib])
    errors={k:float(abs(p[k][ia].astype(float)-q[k][ib].astype(float)).max()) for k in p.dtype.names if k!='u4if'}
    return dict(max_absolute_errors=errors,all_recorded_fields_bitwise_equal=all(v==0 for v in errors.values()),metadata_equal=mp==mq)

# Default RK4 physics must retain its frozen fine preflight result.
fine=root/'core_refinement/fine_integrated'
d,r=execute('default_rk4',fine/'khout.000001.dat',fine/'params.dat',1e-7,1e-7,ssp=False)
assert r['exit_code']==0
manifest['default_comparison']=compare(Path(r['final']),root/'core_refinement/native_v4/preflight_fine/khout.000002.dat')
save();assert manifest['default_comparison']['all_recorded_fields_bitwise_equal'] and manifest['default_comparison']['metadata_equal']

# Uniform flow crosses periodic boundaries and MPI ownership cuts in the
# forced, discarded FE endpoint. Its accepted trajectory must be unaffected.
icdir=dest/'uniform_ic';icdir.mkdir()
template,_=read_snapshot(fine/'khout.000001.dat');n=64;N=n*n
raw=np.zeros(N*template.dtype.itemsize,np.uint8);p=np.ndarray(N,dtype=template.dtype,buffer=raw)
ids=np.arange(N);h=1/n;gamma=float(np.float32(1.666667))
p['u4if']=ids;p['x']=(ids%n+.5)*h;p['y']=(ids//n+.5)*h
p['vx']=.2;p['vy']=.1;p['mass']=h*h;p['volume']=h*h;p['den']=1
p['pressure']=1e-6;p['ie']=float(np.float32(1e-6))*h*h/(gamma-1)
initial=icdir/'khout.000001.dat'
with initial.open('wb') as fh:
    fh.write(struct.pack('<4i',-1,N,template.dtype.itemsize,1));fh.write(raw.tobytes())
    fh.write(struct.pack('<ddfiiddff',0.,1e-7,1e-7,n,n,1.,1.,0.,0.))
params=root/'rk4_diagnostics/validation/sedov_first_failure/params.dat'
uniform={}
for ranks in (1,4):
    for force in (False,True):
        name=f'uniform_mpi{ranks}_'+('reject' if force else 'reference')
        d,r=execute(name,initial,params,.05,.05 if force else .025,ranks=ranks,force=force)
        assert r['exit_code']==0 and r['same_endpoint'] and r['individual_mass_exact'] and r['W0_exact']
        assert r['rollbacks']==int(force) and r['accepted_steps']==2 and r['floor_injection']==0
        pp,_=read_snapshot(Path(r['final']));oo=np.argsort(pp['u4if'])
        exact=np.column_stack(((p['x']+.05*.2)%1,(p['y']+.05*.1)%1))
        delta=np.column_stack((pp['x'][oo],pp['y'][oo]))-exact;delta-=np.rint(delta)
        assert np.max(abs(delta))<1e-12
        uniform[name]=Path(r['final'])
manifest['rollback_mpi1']=compare(uniform['uniform_mpi1_reference'],uniform['uniform_mpi1_reject'])
manifest['rollback_mpi4']=compare(uniform['uniform_mpi4_reference'],uniform['uniform_mpi4_reject'])
manifest['MPI1_MPI4']=compare(uniform['uniform_mpi1_reference'],uniform['uniform_mpi4_reference'])
save()
for key in ('rollback_mpi1','rollback_mpi4','MPI1_MPI4'):
    assert manifest[key]['all_recorded_fields_bitwise_equal'] and manifest[key]['metadata_equal'],(key,manifest[key])

# The two physics intervals use only one coarse/half ceiling pair each.
sources=[('fine',json.loads((root/'rk4_diagnostics/fine_cold_interval/manifest.json').read_text())['source'],[67935,68001,68526,68592]),
         ('sedov',str(root/'rk4_diagnostics/validation/sedov_first_failure/khdiag.step000256.start.dat'),[1699,1829])]
for label,src,focus in sources:
    initial=Path(src);failed=Path(str(initial).replace('.start.dat','.failed.dat'))
    _,m=read_snapshot(initial);_,mf=read_snapshot(failed);interval=mf['dt'];target=m['t']+interval
    paths=[]
    for factor in (1.,.5):
        d,r=execute(label+('_base' if factor==1 else '_half'),initial,initial.parent/'params.dat',target,interval*factor,audit_ids=focus)
        paths.append(r)
        if r['exit_code'] or not r.get('same_endpoint'):
            manifest['classification']='STOP: bounded interval failed; no fresh/long continuation';save();raise SystemExit(2)
        assert r['individual_mass_exact'] and r['W0_exact'] and r['min_U']>0 and r['floor_injection']==0
    manifest[label+'_ceiling_comparison']=compare(Path(paths[0]['final']),Path(paths[1]['final']))
    save()
manifest['classification']='bounded group completed; evaluate local time sensitivity/thermal accuracy before fresh or long continuation'
save();print(manifest['classification'],flush=True)
