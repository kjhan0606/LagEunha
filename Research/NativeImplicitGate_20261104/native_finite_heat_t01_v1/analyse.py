from pathlib import Path
import json,re,struct,numpy as np
root=Path(__file__).resolve().parent;pkg=root.parent
initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_w0_initial_20261101_v1/initial696.bin')
gamma=float(np.float32(1.666667))
def total(a):return float(np.sum(a,dtype=np.longdouble))
def read(p):
    with p.open('rb') as f:
        h=f.read(64);f.seek(-52,2);footer=f.read(52)
    magic,version,n,stride,prefix,features=struct.unpack('<iIQIIQ',h[:32])
    if (magic,version,n,stride,prefix)!=(-2147483607,2,90352,696,576):raise ValueError('Unexpected state ABI')
    if p.stat().st_size!=64+n*stride+52:raise ValueError('Incomplete state')
    raw=np.memmap(p,dtype='u1',mode='r',offset=64,shape=(n,stride))
    def col(off,typ='f8'):return np.ndarray((n,),dtype=typ,buffer=raw,offset=off,strides=(stride,)).copy()
    ids=col(8,'u8')&np.uint64(2**56-1);order=np.argsort(ids)
    if len(np.unique(ids))!=n:raise ValueError('Nonunique IDs')
    d={k:col(off)[order] for k,off in [('x',24),('y',32),('mass',48),('vx',56),('vy',64),('W',576),('pi',584),('U',672),('V',688)]}
    d['id']=ids[order];d['time']=struct.unpack('<d',footer[:8])[0]
    d['r']=np.hypot(d['x']-3,d['y']-3);d['rho']=d['mass']/d['V'];d['P']=(gamma-1)*d['U']/d['V']
    d['cs']=np.sqrt(gamma*d['P']/d['rho']);d['ell']=(d['x']-3)*d['vy']-(d['y']-3)*d['vx']
    d['vr']=np.divide((d['x']-3)*d['vx']+(d['y']-3)*d['vy'],d['r'],out=np.zeros_like(d['r']),where=d['r']>0)
    d['E']=total(d['U']+.5*d['mass']*(d['vx']**2+d['vy']**2)-d['mass']/np.sqrt(d['r']**2+.01175**2))
    return d
ic=read(initial);R0=ic['r']
cells=np.fromfile(pkg/'native_rhs_comparison/t0_unified/weight_native_cells.bin',dtype=[('id','i8'),('v','f8',(10,))]);order=np.argsort(cells['id']);assert np.array_equal(cells['id'][order],ic['id'])
V0=cells['v'][order,9]
summary=[];states={}
for d in sorted(root.glob('*_to*')):
    result=d/'result.json'
    if not result.exists():continue
    result=json.loads(result.read_text())
    if result['exit']:summary.append(result);continue
    outputs=[p for p in d.glob('khout.*.dat') if not p.is_symlink()]
    if len(outputs)!=1:raise ValueError('Expected one output')
    a=read(outputs[0]);assert np.array_equal(a['id'],ic['id'])
    log=(d/'log').read_text();q=np.array([[float(v) for v in row] for row in re.findall(r'UNIFIED_COHORT_HEAT\] interval=\S+ core=(\S+) original_disk=(\S+) original_ambient=(\S+)',log)])
    heat=dict(zip(['original_core','original_disk','original_ambient'],q.sum(axis=0)))
    implicit_steps=re.findall(r'UNIFIED_KDK\] step=\d+ t=(\S+) dt=(\S+).*? total_s=(\S+)',log)
    adiabat=[float(v) for v in re.findall(r'max_adiabat_identity_error=(\S+)',log)]
    row={**result,'actual_time':a['time'],'final':str(outputs[0]),'ids_unchanged':True,'max_mass_change':float(np.max(abs(a['mass']-ic['mass']))),'max_W':float(abs(a['W']).max()),'max_pi':float(abs(a['pi']).max()),'min_U':float(a['U'].min()),'min_V':float(a['V'].min()),'energy_relative_to_t0':(a['E']-ic['E'])/abs(ic['E']), 'angular_momentum_relative_to_t0':(total(a['mass']*a['ell'])-total(ic['mass']*ic['ell']))/abs(total(ic['mass']*ic['ell'])), 'actual_viscous_heat_in_segment':heat,'accepted_steps':len(implicit_steps),'max_adiabat_identity_error':max(adiabat) if adiabat else None,'integrated_step_wall_seconds':sum(float(r[2]) for r in implicit_steps),'cohorts':{}}
    entropy=a['U']/ic['U']*(a['V']/V0)**(gamma-1)
    for label,mask in [('original_core',R0<.5),('original_disk',(R0>=.5)&(R0<=2)),('original_ambient',R0>2),('paper_current_annulus',(a['r']>=.8)&(a['r']<=1.2))]:
        mass=total(a['mass'][mask]);vk=a['r'][mask]/(a['r'][mask]**2+.01175**2)**.75
        row['cohorts'][label]={'N':int(mask.sum()),'U':total(a['U'][mask]),'entropy_equivalent_heat_at_initial_V':total(ic['U'][mask]*(entropy[mask]-1)),'min_entropy_ratio':float(entropy[mask].min()),'max_entropy_ratio':float(entropy[mask].max()), 'vr_mass_rms':float(np.sqrt(total(a['mass'][mask]*a['vr'][mask]**2)/mass)),'max_abs_vr':float(abs(a['vr'][mask]).max()),'L1_relative_vphi':float(np.mean(abs(a['ell'][mask]/a['r'][mask]-vk)/vk)), 'max_relative_density_change_by_id':float(abs(a['V'][mask]/V0[mask]-1).max()),'mean_abs_relative_volume_change_by_id':total(a['mass'][mask]*abs(a['V'][mask]/V0[mask]-1))/mass}
    summary.append(row);states[(result['branch'],result['target'])]=a
    bins=[]
    for left,right in zip(np.linspace(0,2.75,45)[:-1],np.linspace(0,2.75,45)[1:]):
        mask=(a['r']>=left)&(a['r']<right)
        if not mask.any():continue
        mass=total(a['mass'][mask]);vol=total(a['V'][mask])
        bins.append({'r_lo':float(left),'r_hi':float(right),'N':int(mask.sum()),'rho_volume_mean':mass/vol,'P_volume_mean':total(a['P'][mask]*a['V'][mask])/vol,'cs_mass_mean':total(a['cs'][mask]*a['mass'][mask])/mass,'vr_mass_rms':float(np.sqrt(total(a['mass'][mask]*a['vr'][mask]**2)/mass))})
    (d/'radial_profile.json').write_text(json.dumps(bins,indent=2)+'\n')
comparisons=[]
for target in [.01,.03,.1]:
    if ('base',target) not in states or ('half',target) not in states:continue
    a,b=states['base',target],states['half',target];mask=(R0>=.5)&(R0<=2)
    comparisons.append({'time':target,'scope':'same-ID original-disk difference, not an Eulerian exact-solution norm', 'mass_weighted_abs_velocity_difference':total(a['mass'][mask]*np.hypot(a['vx'][mask]-b['vx'][mask],a['vy'][mask]-b['vy'][mask]))/total(a['mass'][mask]), 'mass_weighted_abs_relative_density_difference':total(a['mass'][mask]*abs(a['rho'][mask]/b['rho'][mask]-1))/total(a['mass'][mask])})
(root/'assessment.json').write_text(json.dumps({'scope':'bounded finite-heat gate; no full 2D/3D or matched speed qualification','gamma':gamma,'initial_stored_volume_refresh_max_relative':float(abs(V0/ic['V']-1).max()),'cases':summary,'time_refinement':comparisons},indent=2)+'\n')
for row in summary:
    print(row['branch'],row['target'],'exit',row['exit'],'steps',row.get('accepted_steps'),'dE',row.get('energy_relative_to_t0'))
    if 'cohorts' in row:
        for key in ['original_core','original_disk','paper_current_annulus']:
            v=row['cohorts'][key];print(key,'U',v['U'],'entropy_heat',v['entropy_equivalent_heat_at_initial_V'],'vphiL1',v['L1_relative_vphi'],'vr_rms',v['vr_mass_rms'])
print(comparisons)
