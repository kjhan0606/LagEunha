"""Read-only observations for the native self-consistent midpoint transaction.

No RHS, fit, correction, acceptance, or replacement of solver state. Region
masks use the original exact-support radial envelope and both-state union.
"""
import hashlib
import json
import shutil
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from full_moment_state_io import load,dump,digest,HEADER

REGIONS=('all','inner','cavity','disk_core','inner_contact','outer_contact','dynamic_buffer','wall')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):
    with Path(path).open('x') as f:json.dump(data,f,indent=2);f.write('\n')

def radial_bounds(cell,radius):
    minimum=[];maximum=[];inside=True
    for line in cell.lines:
        A=line.A;d=tuple(b-a for a,b in zip(A,line.B));norm2=sum(x*x for x in d);normal=(d[1],-d[0])
        if sum(n*(-a) for n,a in zip(normal,A))>0:inside=False
        t=max(Fraction(0),min(Fraction(1),-sum(a*x for a,x in zip(A,d))/norm2))
        point=tuple(a+t*x for a,x in zip(A,d))
        minimum.append(sum(x*x for x in point));maximum.extend(sum(x*x for x in q) for q in (A,line.B))
    if cell.arcs:minimum.append(Fraction(float(radius))**2);maximum.append(Fraction(float(radius))**2)
    if inside:minimum.append(Fraction(0))
    return float(np.sqrt(np.longdouble(min(minimum)))),float(np.sqrt(np.longdouble(max(maximum))))

def summary(state_path,geometry_path,initial_path,gamma):
    state_path,geometry_path=Path(state_path),Path(geometry_path)
    state=load(state_path);initial=load(initial_path);ids=state['ids'];n=len(ids)
    with np.load(geometry_path/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    a=state['state'][:-4].reshape(n,20);area=mesh['moments'][:,0]
    assert np.array_equal(ids,initial['ids']) and np.array_equal(state['mass'],initial['mass'])
    assert np.array_equal(mesh['ids'],ids) and np.array_equal(mesh['sites'],a[:,:2]) and np.array_equal(mesh['weights'],a[:,2])
    thermal=a[:,19].astype(np.longdouble)-.5*(a[:,15].astype(np.longdouble)+a[:,18].astype(np.longdouble))
    assert np.isfinite(gamma) and gamma>1 and (area>0).all() and (thermal>0).all()
    sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_circular_curved_modules_v3_20261019')
    from curved_disk_cell import DiskCell
    cells=json.loads((geometry_path/'full_exact_supports.json').read_text())['cells']
    bounds=np.array([radial_bounds(DiskCell.from_native(cell,float(mesh['radius'])),float(mesh['radius'])) for cell in cells])
    lo,hi=bounds.T;wall=mesh['arc_counts']>0;pad=32*np.finfo(float).eps*np.maximum(1,hi)
    ci=(lo-pad<=.5)&(hi+pad>=.5);co=(lo-pad<=2)&(hi+pad>=2)
    r0=np.linalg.norm(initial['state'][:-4].reshape(n,20)[:,:2],axis=1)
    masks=dict(all=np.ones(n,bool),inner=r0<.4,cavity=(hi<.5)&~ci,
        disk_core=(lo>.5)&(hi<2)&~ci&~co,inner_contact=ci,outer_contact=co,
        dynamic_buffer=(lo>2)&~wall&~co,wall=wall)
    return dict(ids=ids,mass=state['mass'],area=area,rho=state['mass'].astype(np.longdouble)/area,
        pressure=(np.longdouble(gamma)-1)*thermal/area,time=np.array(state['time']),
        state_digest=np.array(digest(state_path)),geometry_sha256=np.array(sha(geometry_path/'full_exact_supports.json')),
        face_keys=mesh['face_keys'],arc_counts=mesh['arc_counts'],line_counts=mesh['line_counts'],**masks)

def agreement(provider,proof,output):
    output=Path(output);output.mkdir(exist_ok=False)
    for key in ('ids','mass','time'):assert np.array_equal(provider[key],proof[key]),key
    assert all((provider[key]>0).all() and np.isfinite(provider[key]).all() and (proof[key]>0).all() and np.isfinite(proof[key]).all() for key in ('rho','pressure'))
    rows={};passed=True
    for name in REGIONS:
        mask=provider[name]|proof[name];w=provider['mass'][mask].astype(np.longdouble)
        row=dict(cells=int(mask.sum()),membership_changes=int(np.count_nonzero(provider[name]!=proof[name])),fields={})
        for field in ('rho','pressure'):
            error=(proof[field][mask]-provider[field][mask])/provider[field][mask]
            rms=float(np.sqrt(np.sum(w*error*error)/w.sum())) if len(w) else 0.
            maximum=float(abs(error).max(initial=0.));ok=rms<=1e-8 and maximum<=1e-7
            passed=passed and ok;row['fields'][field]=dict(RMS=rms,maximum=maximum,passed=ok)
        rows[name]=row
    report=dict(status='passed_actual_midpoint_agreement' if passed else 'failed_actual_midpoint_agreement',
        regions=rows,RMS_limit=1e-8,maximum_limit=1e-7,union_of_both_actual_region_masks=True,
        provider_state_digest=str(provider['state_digest']),proof_state_digest=str(proof['state_digest']),
        provider_geometry_sha256=str(provider['geometry_sha256']),proof_geometry_sha256=str(proof['geometry_sha256']),
        physical_providers_added=0,accepted_time_steps=0)
    with (output/'proof_summary.npz').open('xb') as f:np.savez(f,**proof)
    write(output/'report.json',report)
    assert passed,'Actual provider midpoint and accepted arithmetic midpoint disagree'
    return output

def capture_provider(service,request):
    ctx=service.ctx;root=service.root;output=Path(bytes(request.output_dir).decode())
    initial=Path(bytes(request.aux_path).decode());base_path=Path(bytes(request.fits_path).decode())
    gamma=json.loads(Path(ctx['physics_file']).read_text())['physics']['gamma']
    values=summary(ctx['state'],ctx['geometry'],initial,gamma)
    output.mkdir(exist_ok=False)
    with (output/'summary.npz').open('xb') as f:np.savez(f,**values)
    # All arrays are actual values, captured before the RELEASE lifetime ends.
    for name,source in [('base_state',base_path),('provider_state',Path(ctx['state']))]:
        state=load(source);header=HEADER.unpack((source/'header.bin').read_bytes())
        dump(output/name,state,header[3],header[4])
    copied={}
    sources={
        'base_areas.bin':Path(bytes(request.bounds_path).decode())/'areas.bin',
        'actual_stage_report.json':root/'actual_stage_report.json',
        'cell_transport_journal.npz':root/'actual_stage_audit/cell_transport_journal.npz',
        'cell_transport_journal.json':root/'actual_stage_audit/cell_transport_journal.json',
        'generator.bin':root/'rates/generator.bin','weight_rates.bin':root/'rates/weight_rates.bin',
        'Psi_rates.bin':root/'rates/Psi_rates.bin','ledgers.bin':root/'rates/ledgers.bin',
        'stage.bin':root/'native_stage/outputs/stage.bin','density_fits.json':root/'global_fits.json'}
    for name,source in sources.items():
        assert not (output/name).exists();shutil.copy2(source,output/name);copied[name]=sha(output/name)
    journal=json.loads((output/'cell_transport_journal.json').read_text())
    assert journal['stage_key']==ctx['stage_key'] and journal['time']==float(values['time'])
    write(output/'capture.json',dict(status='captured_final_actual_midpoint_before_release',
        state_digest=ctx['state_sha'],geometry_sha256=ctx['geometry_sha'],stage_key=ctx['stage_key'],
        copied_sha256=copied,physical_providers_added=0,accepted_time_steps=0))
    return output

def moment_kinematics(raw):
    """Cellwise affine/spin energies and angular momenta from raw moments."""
    raw=np.asarray(raw,dtype=np.longdouble);m=raw[:,0];q=raw[:,1:3];p=raw[:,7:9]
    Q=raw[:,3:7].reshape(-1,2,2);D=raw[:,9:13].reshape(-1,2,2);Pi=raw[:,13:17].reshape(-1,2,2)
    outer=lambda a,b:np.einsum('ni,nj->nij',a,b)
    trace=lambda a:np.trace(a,axis1=1,axis2=2)
    S=Q-outer(q,q)/m[:,None,None];S=.5*(S+S.transpose(0,2,1))
    C=D-outer(q,p)/m[:,None,None]
    det=S[:,0,0]*S[:,1,1]-S[:,0,1]*S[:,1,0]
    if np.any(m<=0) or np.any(det<=0):raise ValueError('non-positive mass/covariance in spin journal')
    invS=np.empty_like(S);invS[:,0,0]=S[:,1,1]/det;invS[:,1,1]=S[:,0,0]/det
    invS[:,0,1]=-S[:,0,1]/det;invS[:,1,0]=-S[:,1,0]/det
    B=np.einsum('nji,njk->nik',C,invS)
    ell=C[:,0,1]-C[:,1,0];inertia=trace(S);centroid=q/m[:,None]
    orbital=centroid[:,0]*p[:,1]-centroid[:,1]*p[:,0]
    total=D[:,0,1]-D[:,1,0]
    Kaffine=.5*np.einsum('nij,njk,nik->n',B,S,B)
    Krot=ell*ell/(2*inertia);T=Pi-outer(p,p)/m[:,None,None]
    Kcm=.5*np.sum(p*p,axis=1)/m;Kresidual=.5*trace(T)-Kaffine
    thermal=raw[:,17]-.5*trace(Pi)
    return dict(S=S,C=C,B=B,ell=ell,inertia=inertia,orbital_L=orbital,total_L=total,
        Kaffine=Kaffine,Krot=Krot,Kaffine_nonspin=Kaffine-Krot,
        Kcm=Kcm,Kresidual=Kresidual,thermal=thermal)

def moment_tangent(raw,rate):
    """Midpoint derivatives of affine, spin, orbital, and thermal ledgers."""
    raw=np.asarray(raw,dtype=np.longdouble);rate=np.asarray(rate,dtype=np.longdouble)
    kin=moment_kinematics(raw);m=raw[:,0];md=rate[:,0]
    q=raw[:,1:3];qd=rate[:,1:3];p=raw[:,7:9];pd=rate[:,7:9]
    Qd=rate[:,3:7].reshape(-1,2,2);Dd=rate[:,9:13].reshape(-1,2,2);Pid=rate[:,13:17].reshape(-1,2,2)
    outer=lambda a,b:np.einsum('ni,nj->nij',a,b)
    trace=lambda a:np.trace(a,axis1=1,axis2=2)
    S=kin['S'];C=kin['C'];B=kin['B']
    Sdot=Qd-(outer(qd,q)+outer(q,qd))/m[:,None,None]
    Sdot+=outer(q,q)*md[:,None,None]/m[:,None,None]**2
    Cdot=Dd-(outer(qd,p)+outer(q,pd))/m[:,None,None]
    Cdot+=outer(q,p)*md[:,None,None]/m[:,None,None]**2
    inertia_dot=trace(Sdot)
    centroid=q/m[:,None];centroid_dot=(qd-centroid*md[:,None])/m[:,None]
    orbital_dot=(centroid_dot[:,0]*p[:,1]-centroid_dot[:,1]*p[:,0]
        +centroid[:,0]*pd[:,1]-centroid[:,1]*pd[:,0])
    total_dot=Dd[:,0,1]-Dd[:,1,0]
    ell_dot=total_dot-orbital_dot
    omega=kin['ell']/kin['inertia']
    Krot_dot=omega*ell_dot-.5*omega*omega*inertia_dot
    velocity=p/m[:,None]
    Kcm_dot=np.sum(velocity*pd,axis=1)-.5*np.sum(velocity*velocity,axis=1)*md
    Kaffine_dot=np.einsum('nji,nij->n',B,Cdot)-.5*np.einsum('nij,njk,nik->n',B,Sdot,B)
    Tdot=(Pid-(outer(pd,p)+outer(p,pd))/m[:,None,None]
        +outer(p,p)*md[:,None,None]/m[:,None,None]**2)
    Kresidual_dot=.5*trace(Tdot)-Kaffine_dot
    thermal_dot=rate[:,17]-.5*trace(Pid)
    return dict(ell=ell_dot,orbital_L=orbital_dot,total_L=total_dot,
        inertia=inertia_dot,Krot=Krot_dot,Kaffine=Kaffine_dot,
        Kaffine_nonspin=Kaffine_dot-Krot_dot,Kcm=Kcm_dot,
        Kresidual=Kresidual_dot,thermal=thermal_dot,
        energy_partition=Kcm_dot+Kaffine_dot+Kresidual_dot+thermal_dot)

def weighted_rms(value,mass,mask):
    v=np.asarray(value,dtype=np.longdouble)[mask];w=np.asarray(mass,dtype=np.longdouble)[mask]
    return float(np.sqrt(np.sum(w*v*v)/np.sum(w)))

def final_journal(service,request):
    provider_path=Path(bytes(request.bounds_path).decode());output=Path(bytes(request.output_dir).decode())
    base=load(provider_path/'base_state');mid=load(provider_path/'provider_state');end=load(service.ctx['state'])
    n=len(base['ids']);h=float(end['time'])-float(base['time'])
    # The exact integrator h is passed as binary64, not inferred from clocks.
    import struct
    exact_h=struct.unpack('<d',bytes(request.reserved)[:8])[0]
    assert h>0 and abs(h-exact_h)<=8*np.finfo(float).eps*max(1.,abs(h))
    h=exact_h
    with np.load(provider_path/'summary.npz') as f:p={k:f[k] for k in f.files}
    with np.load(provider_path/'cell_transport_journal.npz') as f:journal={k:f[k] for k in f.files}
    with np.load(Path(service.ctx['geometry'])/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    assert np.array_equal(base['ids'],end['ids']) and np.array_equal(journal['ids'],end['ids'])
    base_area=np.fromfile(provider_path/'base_areas.bin',dtype=np.float64);area=mesh['moments'][:,0]
    delta=area.astype(np.longdouble)-base_area.astype(np.longdouble)-np.longdouble(h)*journal['volume']
    stage=np.fromfile(provider_path/'stage.bin',dtype=np.float64).reshape(n,94)
    ar=lambda s:s['state'][:-4].reshape(n,20)
    base_weight=ar(base)[:,2].copy();accepted_weight=ar(end)[:,2].copy()
    weight_rate=np.fromfile(provider_path/'weight_rates.bin',dtype=np.float64)
    assert weight_rate.shape==(n,) and np.isfinite(weight_rate).all()
    weight_delta=accepted_weight.astype(np.longdouble)-base_weight.astype(np.longdouble)
    integrated_weight_rate=np.longdouble(h)*weight_rate.astype(np.longdouble)
    weight_update_defect=weight_delta-integrated_weight_rate
    weight_update_bound=4*np.finfo(float).eps*(
        abs(base_weight.astype(np.longdouble))+abs(accepted_weight.astype(np.longdouble))
        +abs(integrated_weight_rate))
    weight_update_bound=np.maximum(weight_update_bound,np.finfo(float).tiny)
    raw_change=ar(end)[:,3:].astype(np.longdouble)-ar(base)[:,3:].astype(np.longdouble)
    raw_rate=stage[:,1:18].astype(np.longdouble)
    raw_base=np.column_stack((base['mass'],ar(base)[:,3:])).astype(np.longdouble)
    raw_mid=np.column_stack((mid['mass'],ar(mid)[:,3:])).astype(np.longdouble)
    raw_end=np.column_stack((end['mass'],ar(end)[:,3:])).astype(np.longdouble)
    midpoint_rate=np.column_stack((stage[:,0],stage[:,1:18])).astype(np.longdouble)
    kb=moment_kinematics(raw_base);km=moment_kinematics(raw_mid);ke=moment_kinematics(raw_end)
    kd=moment_tangent(raw_mid,midpoint_rate)
    delta={key:ke[key]-kb[key] for key in ('ell','orbital_L','total_L','inertia',
        'Krot','Kaffine','Kaffine_nonspin','Kcm','Kresidual','thermal')}
    rate_residual={key:delta[key]/np.longdouble(h)-kd[key]
        for key in ('ell','orbital_L','total_L','inertia','Krot','Kaffine',
                    'Kaffine_nonspin','Kcm','Kresidual','thermal')}
    spin_regions={};masks={key:np.asarray(p[key],dtype=bool) for key in
        ('all','cavity','inner_contact','disk_core')}
    for name,mask in masks.items():
        row={'cells':int(mask.sum()),'fields':{}}
        for key in delta:
            row['fields'][key]=dict(
                finite_change_sum=float(np.sum(delta[key][mask],dtype=np.longdouble)),
                midpoint_rate_sum=float(np.sum(kd[key][mask],dtype=np.longdouble)),
                finite_change_over_h_minus_midpoint_rate_RMS=weighted_rms(rate_residual[key],base['mass'],mask))
        spin_regions[name]=row
    angular_exchange_closure=delta['ell']+delta['orbital_L']-delta['total_L']
    tangent_energy_closure=kd['energy_partition']-midpoint_rate[:,17]
    spin_rate_residual_summary={key:dict(
        maximum_absolute=str(abs(value).max()),
        mass_weighted_RMS=weighted_rms(value,base['mass'],masks['all']))
        for key,value in rate_residual.items()}
    ledger_L_change=np.longdouble(end['state'][-1])-np.longdouble(base['state'][-1])
    R=lambda s:np.asarray(s['factors'],np.longdouble)@np.asarray(s['factors'],np.longdouble).transpose(0,2,1)
    output.mkdir(exist_ok=False)
    arrays=dict(ids=end['ids'],base_area=base_area,provider_area=p['area'],accepted_area=area,
        Vdot=journal['volume'],delta_V_map=delta,relative_delta_V_map=delta/base_area,
        base_weight=base_weight,accepted_weight=accepted_weight,solved_weight_rate=weight_rate,
        weight_delta=weight_delta,integrated_weight_rate=integrated_weight_rate,
        weight_update_defect=weight_update_defect,weight_update_roundoff_bound=weight_update_bound,
        internal_spin_base=kb['ell'],internal_spin_midpoint=km['ell'],internal_spin_accepted=ke['ell'],
        internal_spin_change=delta['ell'],midpoint_internal_spin_rate=kd['ell'],
        centroid_orbital_L_base=kb['orbital_L'],centroid_orbital_L_accepted=ke['orbital_L'],
        centroid_orbital_L_change=delta['orbital_L'],midpoint_centroid_orbital_L_rate=kd['orbital_L'],
        total_cell_L_base=kb['total_L'],total_cell_L_accepted=ke['total_L'],
        total_cell_L_change=delta['total_L'],midpoint_total_cell_L_rate=kd['total_L'],
        inertia_base=kb['inertia'],inertia_accepted=ke['inertia'],inertia_change=delta['inertia'],
        rigid_spin_energy_base=kb['Krot'],rigid_spin_energy_accepted=ke['Krot'],
        rigid_spin_energy_change=delta['Krot'],midpoint_rigid_spin_energy_rate=kd['Krot'],
        affine_energy_base=kb['Kaffine'],affine_energy_accepted=ke['Kaffine'],
        affine_energy_change=delta['Kaffine'],midpoint_affine_energy_rate=kd['Kaffine'],
        affine_nonspin_energy_base=kb['Kaffine_nonspin'],affine_nonspin_energy_accepted=ke['Kaffine_nonspin'],
        affine_nonspin_energy_change=delta['Kaffine_nonspin'],midpoint_affine_nonspin_energy_rate=kd['Kaffine_nonspin'],
        centroid_kinetic_change=delta['Kcm'],midpoint_centroid_kinetic_rate=kd['Kcm'],
        residual_kinetic_change=delta['Kresidual'],midpoint_residual_kinetic_rate=kd['Kresidual'],
        thermal_energy_change=delta['thermal'],midpoint_thermal_energy_rate=kd['thermal'],
        **{f'spin_rate_residual_{key}':value for key,value in rate_residual.items()},
        raw_change=raw_change,midpoint_raw_rate=raw_rate,raw_change_minus_h_rate=raw_change-h*raw_rate,
        factor_covariance_change=R(end)-R(base),
        Psi_change=end['Psi']-base['Psi'],Psi_rate=np.fromfile(provider_path/'Psi_rates.bin',dtype=np.float64),
        ledger_change=end['state'][-4:]-base['state'][-4:])
    with (output/'cells.npz').open('xb') as f:np.savez(f,**arrays)
    oldkeys={tuple(row) for row in p['face_keys']};newkeys={tuple(row) for row in mesh['face_keys']}
    report=dict(status='captured_verified_map_cell_journal',h=h,base_time=float(base['time']),time=float(end['time']),
        base_state_digest=digest(provider_path/'base_state'),provider_state_digest=digest(provider_path/'provider_state'),
        mapped_state_digest=digest(service.ctx['state']),provider_geometry_sha256=str(p['geometry_sha256']),
        accepted_geometry_sha256=service.ctx['geometry_sha'],faces_removed=len(oldkeys-newkeys),faces_added=len(newkeys-oldkeys),
        wall_arc_count_changes=int(np.count_nonzero(p['arc_counts']!=mesh['arc_counts'])),
        signed_delta_V_map=str(delta.sum()),maximum_absolute_delta_V_map=str(abs(delta).max()),
        weight_rate_gauge_mean=float(np.mean(weight_rate,dtype=np.longdouble)),
        maximum_absolute_weight_update_defect=str(abs(weight_update_defect).max()),
        maximum_weight_update_roundoff_bound=str(weight_update_bound.max()),
        maximum_weight_update_bound_fraction=float(np.max(abs(weight_update_defect)/weight_update_bound)),
        weight_update_roundoff_passed=bool(np.all(abs(weight_update_defect)<=weight_update_bound)),
        angular_exchange_closure=dict(
            maximum_absolute_spin_plus_orbit_minus_total_change=str(abs(angular_exchange_closure).max()),
            summed_spin_change=str(np.sum(delta['ell'],dtype=np.longdouble)),
            summed_orbital_L_change=str(np.sum(delta['orbital_L'],dtype=np.longdouble)),
            summed_total_cell_L_change=str(np.sum(delta['total_L'],dtype=np.longdouble)),
            global_ledger_L_change=str(ledger_L_change),
            total_cell_L_minus_global_ledger_L=str(np.sum(delta['total_L'],dtype=np.longdouble)-ledger_L_change)),
        midpoint_energy_partition=dict(
            maximum_absolute_rate_closure=str(abs(tangent_energy_closure).max()),
            summed_partition_rate=str(np.sum(kd['energy_partition'],dtype=np.longdouble)),
            summed_hydro_energy_rate=str(np.sum(midpoint_rate[:,17],dtype=np.longdouble))),
        spin_midpoint_rate_residuals=spin_rate_residual_summary,spin_regions=spin_regions,
        provider_density_fit_sha256=sha(provider_path/'density_fits.json'),accepted_density_fit_sha256=sha(service.root/'global_fits.json'),
        physical_providers_added=0,accepted_time_steps=0,
        limits='Actual signed data, not a certified geometry error enclosure. Pi map corrections retained, not removed.')
    write(output/'report.json',report);return output
