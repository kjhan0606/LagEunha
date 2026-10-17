"""Read-only higher-quadrature, pointwise and kinetic/entropy accuracy audit.

L1 acceptance retains the established physical scales. Sampled maxima are
reported separately and are not certified continuum supremum bounds.
"""
import hashlib,json,os
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from simple_wave_evolution import install_domain
from compare_physical_budget import scales_and_goals
from analyze_simple_wave_time import sampled_fields
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def environment(cfg):
    from density_quadrature import validate_order
    density_order=validate_order(cfg.get('density_quadrature_order',6))
    name=cfg.get('velocity_reconstruction','affine')
    modes=dict(quadratic_budget_continuous='dominant',quadratic_budget_projected='projected',quadratic_budget_projected_regularized='projected_regularized')
    if name!='affine' and name not in modes:raise ValueError('unsupported audit velocity closure')
    pressure=cfg.get('pressure_reconstruction','linear');degrees=dict(linear=2,quadratic_cell_mean=2,polynomial_cell_mean_degree3=3,polynomial_cell_mean_degree4=4)
    if pressure not in degrees:raise ValueError('unsupported audit pressure closure')
    geometry=cfg.get('power_cell_evaluation_backend','active')
    if geometry not in ('active','seeded_active','exhaustive'):raise ValueError('unsupported audit power cell backend')
    reflection=cfg.get('pressure_reflection_y','')
    velocity_fit=cfg.get('velocity_fit','radial_fixed_mean')
    if velocity_fit not in ('radial_fixed_mean','tensor_free_mean','tensor_free_affine'):raise ValueError('unsupported audit velocity fit')
    return dict(LAGFORCE_DENSITY_QUADRATURE_ORDER=str(density_order),LAGFORCE_PRESSURE_FIT=cfg.get('pressure_fit','radial_fixed_mean'),LAGFORCE_PRESSURE_REFLECT_Y=reflection,LAGFORCE_QUADRATIC_VELOCITY='0' if name=='affine' else '1',LAGFORCE_QUADRATIC_BUDGET=modes.get(name,'dominant'),
        LAGFORCE_VELOCITY_FIT=velocity_fit,LAGFORCE_QUADRATIC_PRESSURE='0' if pressure=='linear' else '1',LAGFORCE_PRESSURE_DEGREE=str(degrees[pressure]),
        LAGFORCE_FAST_PRESSURE_POLYNOMIAL='1' if cfg.get('pressure_evaluation_backend','monomial')=='horner' else '0',
        LAGFORCE_ACTIVE_POWER_CELLS='0' if geometry=='exhaustive' else '1',LAGFORCE_SEEDED_POWER_CELLS='1' if geometry=='seeded_active' else '0',
        LAGFORCE_VECTOR_SMOOTH_PATH='1' if cfg.get('mass_path_evaluation_backend','scalar')=='cell_batched' else '0',LAGFORCE_VECTOR_MOMENT_FLUX='1' if cfg.get('moment_flux_evaluation_backend','scalar')=='batched' else '0',
        LAGFORCE_VECTOR_QUADRATIC_BUDGET='1' if cfg.get('quadratic_budget_backend','scalar')=='batched' else '0')


def errors(numerical,reference,weights,scales):
    weights=np.asarray(weights);out={}
    if np.any(weights<0):raise ValueError('positive audit quadrature required')
    for key,scale in scales.items():
        delta=np.asarray(numerical[key])-np.asarray(reference[key]);norm=np.linalg.norm(delta,axis=1) if delta.ndim==2 else abs(delta)
        if scale<=0 or not np.all(np.isfinite(norm)):raise ValueError('invalid audit field or scale')
        out[key]=dict(L1=float(weights@norm/scale),L2_squared=float(weights@(norm*norm)/scale**2),maximum_sampled=float(np.max(norm)/scale))
    return out


def merge(left,right):
    for key,row in right.items():
        if key not in left:left[key]=dict(L1=0.,L2_squared=0.,maximum_sampled=0.)
        left[key]['L1']+=row['L1'];left[key]['L2_squared']+=row['L2_squared'];left[key]['maximum_sampled']=max(left[key]['maximum_sampled'],row['maximum_sampled'])


def finish_norms(rows):
    return {k:dict(L1=v['L1'],L2=float(np.sqrt(v['L2_squared'])),maximum_sampled=v['maximum_sampled']) for k,v in rows.items()}


def reconstructed_raw(points,weights,rho,velocity,pressure,residual_covariance,gamma):
    dm=np.asarray(weights)*rho;m=float(dm.sum());p=dm@velocity;Pi=velocity.T@(dm[:,None]*velocity)+m*residual_covariance
    return np.r_[m,dm@points,(points.T@(dm[:,None]*points)).ravel(),p,(points.T@(dm[:,None]*velocity)).ravel(),Pi.ravel(),weights@pressure/(gamma-1)+.5*np.trace(Pi)]


def inspect(entry):
    install_domain();root=Path(entry['root']);cfg=json.loads((root/'simple_wave_config.json').read_text());cp=json.loads((root/'checkpoint.json').read_text());initial=json.loads((root/'initial_state.json').read_text());report=json.loads((root/'simple_wave_report.json').read_text())
    identity=hashlib.sha256((root/'simple_wave_config.json').read_bytes()).hexdigest()
    if cp['config_sha256']!=identity or initial['config_sha256']!=identity or report['config']!=cfg:raise ValueError('audit source identity mismatch')
    if not report['completed'] or not report['numerical_passed'] or abs(cp['t']-cfg['end'])>1e-14:raise ValueError('audit endpoint incomplete or numerical gate failed')
    np.testing.assert_array_equal(cp['mass'],initial['mass']);os.environ.update(environment(cfg))
    mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20);L=[np.asarray(v) for v in cp['factors']]
    cells,edges=model.geometry(a[:,:2],a[:,2],cp['t']);states=[model.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    quadratic=os.environ['LAGFORCE_QUADRATIC_VELOCITY']=='1'
    from quadratic_velocity_reconstruction import prepare,correction
    if quadratic:prepare(states,cells,edges,mass,L)
    centers,grad=kinetic.pressure_gradients(states,edges);pressure_polys=None
    if os.environ['LAGFORCE_QUADRATIC_PRESSURE']=='1':
        from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value
        pressure_polys,_=prepare_pressure(states,cells)
    problem=cfg.get('problem','nonlinear_simple_wave')
    if problem=='steady_vortex':
        from steady_vortex_exact import physical_raw,fields
    elif problem=='nonlinear_simple_wave':
        from simple_wave_exact import physical_raw,fields
    else:raise ValueError('unknown physical audit benchmark')
    scales,goals=scales_and_goals(cfg);all_scales=dict(scales,energy_with_residual=scales['pressure'],energy_mean=scales['pressure'],entropy_relative=1.)
    residual_cov=[];stored_kinetic=.5*np.sum(a[:,15]+a[:,18]);resolved_kinetic=0.;residual_kinetic=0.
    for s,l,m in zip(states,L,mass):
        T=s.get('velocity_quadratic',{}).get('residual_map',np.eye(2));factor=l@T;cov=factor@factor.T/m;residual_cov.append(cov);residual_kinetic+=.5*m*np.trace(cov)
        points=s['density'].points;v=s['v']+(points-s['c'])@s['B'].T
        if quadratic:v+=correction(s,points)
        resolved_kinetic+=.5*(s['density'].weighted_mass@np.sum(v*v,axis=1))
    quadratures=[]
    for order in entry.get('orders',[16,24]):
        totals={};represented_mass=exact_kinetic=tensor_error=represented_kinetic=represented_energy=0.;min_density=min_pressure=float('inf');vertex_min_density=float('inf');identities=dict(density_moments=0.,momentum=0.,mixed_moment=0.,kinetic_tensor=0.)
        for i,(poly,s) in enumerate(zip(cells,states)):
            raw,points,w,rho,v,P=physical_raw(poly,cp['t'],cfg['amplitude'],order)
            rn=s['density'].value(points);vn=s['v']+(points-s['c'])@s['B'].T
            if quadratic:vn+=correction(s,points)
            pn=s['P']+(points-centers[i])@grad[i] if pressure_polys is None else pressure_value(pressure_polys[i],points)
            rebuilt=reconstructed_raw(points,w,rn,vn,pn,residual_cov[i],model.GAMMA);stored=np.r_[mass[i],a[i,3:]];c=s['c'];h=np.sqrt(s['density'].area)
            def centered_density(r):return r[3:7].reshape(2,2)-np.outer(c,r[1:3])-np.outer(r[1:3],c)+r[0]*np.outer(c,c)
            identities['density_moments']=max(identities['density_moments'],abs(rebuilt[0]-stored[0])/mass[i],np.linalg.norm((rebuilt[1:3]-rebuilt[0]*c)-(stored[1:3]-stored[0]*c))/(mass[i]*h),np.linalg.norm(centered_density(rebuilt)-centered_density(stored))/(mass[i]*h*h))
            identities['momentum']=max(identities['momentum'],float(np.linalg.norm(rebuilt[7:9]-stored[7:9])/(mass[i]*scales['velocity'])))
            mixed=(rebuilt[9:13]-stored[9:13]).reshape(2,2)-np.outer(c,rebuilt[7:9]-stored[7:9]);identities['mixed_moment']=max(identities['mixed_moment'],float(np.linalg.norm(mixed)/(mass[i]*h*scales['velocity'])))
            identities['kinetic_tensor']=max(identities['kinetic_tensor'],float(np.linalg.norm(rebuilt[13:17]-stored[13:17])/(mass[i]*scales['velocity']**2)))
            represented_kinetic+=.5*np.trace(rebuilt[13:17].reshape(2,2));represented_energy+=rebuilt[17]
            en=pn/(model.GAMMA-1)+.5*rn*np.sum(vn*vn,axis=1);er=P/(model.GAMMA-1)+.5*rho*np.sum(v*v,axis=1)
            entropy=pn/rn**model.GAMMA;entropy_reference=P/rho**model.GAMMA
            numerical=dict(density=rn,velocity=vn,pressure=pn,energy_mean=en,energy_with_residual=en+.5*rn*np.trace(residual_cov[i]),entropy_relative=entropy/entropy_reference)
            reference=dict(density=rho,velocity=v,pressure=P,energy_mean=er,energy_with_residual=er,entropy_relative=np.ones_like(rho))
            merge(totals,errors(numerical,reference,w,all_scales));represented_mass+=w@rn;exact_kinetic+=.5*np.trace(raw[13:17].reshape(2,2));tensor_error+=np.linalg.norm(a[i,15:19].reshape(2,2)-raw[13:17].reshape(2,2))
            min_density=min(min_density,float(rn.min()));min_pressure=min(min_pressure,float(pn.min()));vertex_min_density=min(vertex_min_density,float(s['density'].value(poly).min()))
        quadratures.append(dict(order=order,norms=finish_norms(totals),density_min_sampled=min_density,density_min_on_vertices=vertex_min_density,pressure_min_sampled=min_pressure,
            represented_density_mass_relative_error=float(abs(represented_mass-mass.sum())/mass.sum()),kinetic_tensor_L1_over_velocity_scale_squared=float(tensor_error/scales['velocity']**2),reference_kinetic_energy=float(exact_kinetic),
            reconstruction_maximum_cell_identity_errors=identities,represented_kinetic_energy=float(represented_kinetic),high_quadrature_kinetic_relative_identity_error=float(abs(represented_kinetic-stored_kinetic)/max(abs(stored_kinetic),1e-100)),represented_total_energy=float(represented_energy),high_quadrature_total_energy_relative_identity_error=float(abs(represented_energy-a[:,19].sum())/abs(a[:,19].sum()))))
    probes=[]
    for n in entry.get('probe_grids',[256,512]):
        x,y=np.meshgrid((np.arange(n)+.3819660112501051)/n,(np.arange(n)+.6180339887498949)/n);points=np.column_stack([x.ravel(),y.ravel()]);numerical=sampled_fields(cp,points);rho,v,P=fields(points,cp['t'],cfg['amplitude'])
        norms=finish_norms(errors(numerical,dict(density=rho,velocity=v,pressure=P),np.full(len(points),1/len(points)),scales));probes.append(dict(grid=n,norms=norms))
    result=dict(completed=True,label=entry['label'],root=str(root),n=cfg['n'],t=cp['t'],problem=problem,physical_scales=scales,physical_goals=goals,
        original_driver_physical_accuracy_passed=report['physical_accuracy_passed'],higher_quadrature=quadratures,eulerian_probes=probes,
        primary_L1_targets_passed_in_all_audits=bool(all(row['norms'][k]['L1']<goals[k] for row in quadratures+probes for k in goals)),
        fixed_original_mass_array_exact=True,stored_kinetic_energy=float(stored_kinetic),numerical_mean_kinetic_energy=float(resolved_kinetic),unresolved_sample_kinetic_energy=float(residual_kinetic),
        unresolved_fraction_of_reference_kinetic_energy=float(residual_kinetic/quadratures[-1]['reference_kinetic_energy']),
        kinetic_closure_relative_identity_error=float(abs(stored_kinetic-resolved_kinetic-residual_kinetic)/max(abs(stored_kinetic),1e-100)),
        checkpoint_sha256=hashlib.sha256((root/'checkpoint.json').read_bytes()).hexdigest(),
        note='Completed immutable numerical states only; no evolution repeated. Established acceptance is L1 of density/velocity/pressure. Higher Gaussian orders and separate Eulerian probe grids are numerical cross-checks, not analytic integration enclosures. Sampled maxima and vertex minima expose local errors; they do not assert a continuum supremum bound. Entropy and residual kinetic energy are additional diagnostics, not silently changed acceptance targets.')
    atomic(entry['label']+'_field_accuracy_audit.json',result);print(entry['label'],result['primary_L1_targets_passed_in_all_audits'],flush=True);return result


if __name__=='__main__':
    cfg=json.loads(Path('physical_field_audit_config.json').read_text());rows=[inspect(v) for v in cfg['cases']];atomic('physical_field_accuracy_audit_report.json',dict(completed=True,rows=rows))
