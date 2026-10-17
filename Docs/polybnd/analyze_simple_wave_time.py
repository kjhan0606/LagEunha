"""Read-only uniform-time comparison, including reconstructed physical fields."""
import hashlib,json
import os
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.strtree import STRtree
from adaptive_coupled_mass import atomic
from simple_wave_evolution import install_domain
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def check_environment(cfg):
    from quadratic_velocity_reconstruction import closure_name,budget_backend_name
    from quadratic_pressure_reconstruction import closure_name as pressure_name
    from polynomial_pressure_reconstruction import evaluation_backend_name
    velocity=closure_name() if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1' else 'affine'
    active=os.environ.get('LAGFORCE_ACTIVE_POWER_CELLS','0')=='1';default_geometry='active' if active else 'exhaustive';geometry='seeded_active' if os.environ.get('LAGFORCE_SEEDED_POWER_CELLS','0')=='1' else default_geometry
    for key,actual,default in [('velocity_reconstruction',velocity,'affine'),('pressure_reconstruction',pressure_name(),'linear'),('quadratic_budget_backend',budget_backend_name(),'scalar'),('pressure_evaluation_backend',evaluation_backend_name(),'monomial'),('power_cell_evaluation_backend',geometry,default_geometry)]:
        if cfg.get(key,default)!=actual:raise ValueError('time field evaluation environment mismatch: '+key)


def observed(differences,screen):
    if min(differences)<=screen:return None
    return float(np.log2(differences[0]/differences[1]))


def sampled_fields(cp,points):
    mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],cp['t']);polys=np.asarray([Polygon(v) for v in cells],dtype=object)
    states=[model.decode(v,np.r_[m,r[3:]]) for v,m,r in zip(cells,mass,a)]
    quadratic=os.environ.get('LAGFORCE_QUADRATIC_VELOCITY','0')=='1'
    if quadratic:
        from quadratic_velocity_reconstruction import prepare,correction
        prepare(states,cells,edges,mass,[np.asarray(v) for v in cp['factors']])
    centers,grad=kinetic.pressure_gradients(states,edges)
    from quadratic_pressure_reconstruction import enabled as quadratic_pressure_enabled
    pressure_polynomials=None
    if quadratic_pressure_enabled():
        from quadratic_pressure_reconstruction import prepare as prepare_pressure,value as pressure_value
        pressure_polynomials,_=prepare_pressure(states,cells)
    pairs=STRtree(polys).query(shapely.points(points),predicate='covered_by')
    coverage=np.bincount(pairs[0],minlength=len(points))
    if not np.all(coverage==1):raise RuntimeError('time comparison probe partition failed')
    rho=np.zeros(len(points));velocity=np.zeros_like(points);pressure=np.zeros(len(points))
    for i,s in enumerate(states):
        ids=pairs[0,pairs[1]==i];p=points[ids]
        rho[ids]=s['density'].value(p);velocity[ids]=s['v']+(p-s['c'])@s['B'].T
        if quadratic:velocity[ids]+=correction(s,p)
        pressure[ids]=s['P']+(p-centers[i])@grad[i]
        if pressure_polynomials is not None:pressure[ids]=pressure_value(pressure_polynomials[i],p)
    return dict(density=rho,velocity=velocity,pressure=pressure)


def run():
    install_domain();cfg=json.loads(Path('time_analysis_config.json').read_text());states=[];reports=[];configs=[];cps=[]
    for root in map(Path,cfg['roots']):
        c=json.loads((root/'simple_wave_config.json').read_text());r=json.loads((root/'simple_wave_report.json').read_text());cp=json.loads((root/'checkpoint.json').read_text())
        check_environment(c)
        if not r['numerical_passed'] or not c['fixed_dt'] or r['rejected']!=0:raise ValueError('nonuniform or failed time branch')
        start=r['source_time'];end=c['end'];tol=64*np.finfo(float).eps*max(abs(start),abs(end))
        if abs(cp['t']-end)>tol:raise ValueError('time comparison endpoint mismatch')
        for k in range(1,r['steps']+1):
            step=json.loads((root/('step_'+str(k)+'_report.json')).read_text())
            if abs(step['dt']-c['max_dt'])>tol:raise ValueError('actual accepted dt was adaptive')
        states.append(np.asarray(cp['state'])[:-4].reshape(-1,20));reports.append(r);configs.append(c);cps.append(cp)
    if len(states)!=3:raise ValueError('three time levels required')
    for key in ('n','amplitude','end','preserved_initial_state_sha256','preserved_source_rhs_sha256','mass_root_backend','velocity_reconstruction','pressure_reconstruction','geometry_path'):
        if len({c.get(key,'quadratic_midpoint' if key=='geometry_path' else None) for c in configs})!=1:raise ValueError('time branch source/model mismatch: '+key)
    for a,b in zip(configs[:-1],configs[1:]):
        if a['max_dt']!=2*b['max_dt']:raise ValueError('dt halving required')
    mass=np.asarray(cps[0]['mass']);amp=configs[0]['amplitude']
    for cp in cps:np.testing.assert_array_equal(cp['mass'],mass)
    screen=256*np.finfo(float).eps/amp;rows={}
    for name,sl in dict(sites=slice(0,2),weights=slice(2,3),q=slice(3,5),Q=slice(5,9),p=slice(9,11),D=slice(11,15),Pi=slice(15,19),E=slice(19,20)).items():
        diffs=[]
        for a,b in zip(states[:-1],states[1:]):
            delta=a[:,sl]-b[:,sl]
            if name=='weights':delta-=delta[-1]
            value=float(np.max(abs(delta))/amp) if name in ('sites','weights') else float(np.sum(abs(delta))/(mass.sum()*amp))
            diffs.append(value)
        rows[name]=dict(differences=diffs,roundoff_screen=screen,observed_order=observed(diffs,screen))
    x,y=np.meshgrid((np.arange(128)+.5)/128,(np.arange(128)+.5)/128);points=np.column_stack([x.ravel(),y.ravel()])
    fields=[sampled_fields(cp,points) for cp in cps];field_rows={}
    for key in fields[0]:
        diffs=[]
        for a,b in zip(fields[:-1],fields[1:]):
            delta=a[key]-b[key]
            value=np.mean(np.linalg.norm(delta,axis=1)) if delta.ndim==2 else np.mean(abs(delta))
            diffs.append(float(value/amp))
        field_rows[key]=dict(differences_over_amplitude=diffs,roundoff_screen=screen,observed_order=observed(diffs,screen))
    defects=[r['integrated_absolute_mass_defect'] for r in reports]
    result=dict(completed=True,source_time=reports[0]['source_time'],end=configs[0]['end'],actual_steps=[r['steps'] for r in reports],
        actual_dt=[c['max_dt'] for c in configs],all_numerical_gates_passed=True,state=rows,fields=field_rows,
        integrated_absolute_mass_defects=defects,mass_defect_adjacent_ratios=[a/b for a,b in zip(defects[:-1],defects[1:])],
        fine_temporal_field_difference_below_physical_budget=bool(max(v['differences_over_amplitude'][-1] for v in field_rows.values())<.001),
        note='All actual accepted dt checked; common preserved state/RHS; no reintegration. Physical fields sampled on the same16384 Eulerian points. Orders below conservative roundoff screen omitted. This is temporal self-convergence, not spatial Euler accuracy.')
    atomic('simple_wave_time_analysis.json',result);print(result,flush=True)


if __name__=='__main__':run()
