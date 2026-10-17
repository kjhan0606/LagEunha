"""Transport a numerical reconstruction, then solve its fixed-mass capacities.

This stores the full subcell field and is a distinct affine-contact comparison,
not an implementation of the six-moment closure or a general pressure solver.
"""
import hashlib,json,time
from pathlib import Path
import numpy as np
import shapely
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from predictive_boundary_density import BoundaryDensity
from saved_curved_contact import checked_json,initialize_saved,reference_raw,advected_boundary
from transported_density_capacity import PolygonDensity,solve_capacity
import moment_ale_evolve as model


def numerical_affine_flow(state,mass,areas):
    a=np.asarray(state)[:-4].reshape(-1,20);mass=np.asarray(mass);total=mass.sum()
    q=a[:,3:5].sum(0);Q=a[:,5:9].reshape(-1,2,2).sum(0)
    p=a[:,9:11].sum(0);D=a[:,11:15].reshape(-1,2,2).sum(0)
    S=Q-np.outer(q,q)/total;C=D-np.outer(q,p)/total
    B=np.linalg.solve(S,C).T;b=(p-B@q)/total
    qi=a[:,3:5];Qi=a[:,5:9].reshape(-1,2,2)
    expected_p=mass[:,None]*b+qi@B.T
    expected_D=np.einsum('ni,j->nij',qi,b)+Qi@B.T
    pressure=(model.GAMMA-1)*(a[:,19]-.5*(a[:,15]+a[:,18]))/areas
    diag=dict(momentum_affine_defect=float(np.max(abs(a[:,9:11]-expected_p)/mass[:,None])),
        mixed_affine_defect=float(np.max(abs(a[:,11:15].reshape(-1,2,2)-expected_D)/mass[:,None,None])),
        pressure_variation=float(np.ptp(pressure)),nilpotent_defect=float(np.max(abs(B@B))),
        stationary_defect=float(np.max(abs(B@b))),divergence=float(abs(np.trace(B))))
    if max(diag.values())>1e-10 or np.min(pressure)<=0:raise ValueError('not a uniform-pressure stationary affine shear: '+str(diag))
    return b,B,float(np.mean(pressure)),diag


def density_error(field,high,total):
    inside=shapely.area(shapely.intersection(field.polygons,high));area=shapely.area(field.polygons)
    error=np.abs(field.densities-100.)@inside+np.abs(field.densities-1.)@np.maximum(0.,area-inside)
    return float(error/total)


def run():
    started=time.monotonic();entry=json.loads(Path('curved_transport_config.json').read_text())
    cfg={k:entry[k] for k in ('initial_state','initial_state_sha256','reference','reference_sha256')}
    cfg.update(kind='saved_curved_contact',n=entry['n'])
    target,_=initialize_saved(cfg);source=checked_json(entry['source'],entry['source_sha256'])
    np.testing.assert_array_equal(source['mass'],target)
    state=np.asarray(source['state']);a=state[:-4].reshape(-1,20);t0=float(source['t'])
    cells,edges=model.geometry(a[:,:2],a[:,2],t0)
    b,B,P,affine=numerical_affine_flow(state,target,np.asarray([Polygon(v).area for v in cells]))
    reference=checked_json(entry['reference'],entry['reference_sha256'])
    field_path=Path('retained_field.npz');metadata_path=Path('retained_field.json')
    signature=hashlib.sha256(json.dumps(dict(source=entry['source_sha256'],options=entry['density_provider']),sort_keys=True).encode()).hexdigest()
    if field_path.exists():
        metadata=json.loads(metadata_path.read_text())
        if metadata['signature']!=signature or metadata['sha256']!=hashlib.sha256(field_path.read_bytes()).hexdigest():
            raise ValueError('retained field identity mismatch')
        data=np.load(field_path);field=PolygonDensity(data['vertices'],data['densities'])
        closure=metadata['closure']
    else:
        provider=BoundaryDensity(**entry['density_provider'])
        _,closure=provider(cells,edges,np.column_stack([target,a[:,3:]]))
        vertices=np.concatenate([v for v,rho in provider.last_cells]);densities=np.concatenate([rho for v,rho in provider.last_cells])
        field=PolygonDensity(vertices,densities)
        np.savez_compressed(field_path,vertices=vertices,densities=densities)
        atomic(metadata_path,dict(signature=signature,sha256=hashlib.sha256(field_path.read_bytes()).hexdigest(),
            closure=closure,source_sha256=entry['source_sha256'],triangles=len(vertices)))
    validation_path=Path('retained_field_validation.json')
    validation_identity=dict(signature=signature,field_sha256=hashlib.sha256(field_path.read_bytes()).hexdigest())
    if validation_path.exists():
        validation=json.loads(validation_path.read_text())
        if validation['identity']!=validation_identity:raise ValueError('retained field validation identity mismatch')
        initial_mass_error=validation['mass_relative'];initial_raw_error=validation['raw_absolute']
        initial_density=validation['density_L1_relative']
    else:
        initial_mass,initial_coverage,initial_raw=field.integrate([Polygon(v) for v in cells],moments=True)
        initial_mass_error=float(np.max(abs(initial_mass-target)/target))
        initial_raw_error=float(np.max(abs(initial_raw-np.column_stack([target,a[:,3:9]]))))
        initial_density=density_error(field,Polygon(advected_boundary(reference['reference_boundary'],t0)),target.sum())
        atomic(validation_path,dict(identity=validation_identity,mass_relative=initial_mass_error,
            raw_absolute=initial_raw_error,density_L1_relative=initial_density))
    if initial_mass_error>1e-9 or initial_raw_error>1e-10:raise RuntimeError('retained field does not reproduce the source moments')
    weights=a[:,2].copy();rows=[]
    for t in entry['times']:
        if t<=t0:raise ValueError('transport targets must follow source')
        output=Path('t'+format(t,'.8g')+'_report.json');checkpoint=Path('t'+format(t,'.8g')+'_checkpoint.json')
        if output.exists():
            row=json.loads(output.read_text())
            if not row['numerical_passed']:raise RuntimeError('failed transport target requires diagnosis')
            old=json.loads(checkpoint.read_text());weights=np.asarray(old['weights']);rows.append(row);continue
        dt=t-t0;F=np.eye(2)+dt*B;shift=dt*b
        moved=field.affine_image(F,shift);sites=a[:,:2]@F.T+shift;domain=model.domain(t0)@F.T+shift
        pending=Path('t'+format(t,'.8g')+'_pending.json')
        if pending.exists():
            saved=json.loads(pending.read_text())
            if saved['t']!=t or saved['signature']!=signature:raise ValueError('transport pending identity mismatch')
            weights=np.asarray(saved['weights'])
        def save(w,c,h):
            atomic(pending,dict(t=t,weights=w,history=h,signature=signature))
            print('capacity',entry['n'],t,h[-1]['iteration'],h[-1]['mass_relative'],flush=True)
        weights,newcells,history=solve_capacity(sites,weights,domain,target,moved,callback=save)
        mass,coverage,raw=moved.integrate(newcells,moments=True)
        exact,geometry,high=reference_raw([np.asarray(p.exterior.coords[:-1]) for p in newcells],reference['reference_boundary'],t)
        volume=density_error(moved,high,target.sum());physical_mass=float(abs(exact[:,0]-target).sum()/target.sum())
        density_drift=abs(volume-initial_density)
        global_error=float(abs(mass.sum()-field.total_mass)/target.sum())
        reference_mass_error=float(abs(exact[:,0].sum()-target.sum())/target.sum())
        generator_defect=float(np.max(abs((sites-a[:,:2])/dt-(a[:,:2]@B.T+b))))
        numerical=history[-1]['mass_relative']<5e-11 and history[-1]['coverage_relative']<1e-9 and global_error<1e-11 and reference_mass_error<1e-11 and density_drift<1e-10 and generator_defect<1e-10
        accuracy=volume<entry['density_goal'] and physical_mass<entry['cell_mass_goal']
        row=dict(n=entry['n'],t=t,source_time=t0,numerical_passed=bool(numerical),physical_accuracy_passed=bool(accuracy),
            density_L1_relative=volume,cell_mass_L1_relative=physical_mass,
            q_L1_over_total_mass=float(abs(raw[:,1:3]-exact[:,1:3]).sum()/target.sum()),
            Q_L1_over_total_mass=float(abs(raw[:,3:7]-exact[:,3:7]).sum()/target.sum()),
            density_L1_change_from_source=float(density_drift),capacity_mass_relative=history[-1]['mass_relative'],
            coverage_relative=history[-1]['coverage_relative'],generator_velocity_defect=generator_defect,
            numerical_global_mass_relative=global_error,reference_global_mass_relative=reference_mass_error,
            retained_triangles=len(field.polygons),newton_iterations=len(history)-1,pressure=P,
            jacobian_roundoff_enclosures=len(moved.jacobian_enclosures),
            maximum_jacobian_spatial_enclosure=max((v['spatial_enclosure'] for v in moved.jacobian_enclosures),default=0.),
            note='Numerical subcell field retained and transported under its own stationary affine velocity. Original fixed capacities; no generator steering or physical interface in transport/capacity. Pressure-force evolution outside this equilibrium is not tested.')
        atomic(checkpoint,dict(t=t,sites=sites,weights=weights,domain=domain,target_mass=target,
            integrated_mass=mass,position_raw_moments=raw,source_field=str(field_path),source_signature=signature))
        atomic(Path('t'+format(t,'.8g')+'_jacobian_enclosures.json'),moved.jacobian_enclosures)
        atomic(output,row);rows.append(row)
        if not numerical:raise RuntimeError('numerical transport gate failed: '+str(row))
        pending.replace(Path('t'+format(t,'.8g')+'_solved.json'))
    report=dict(completed=True,n=entry['n'],numerical_passed=all(r['numerical_passed'] for r in rows),
        physical_accuracy_passed=all(r['physical_accuracy_passed'] for r in rows),
        source_time=t0,source_density_L1_relative=initial_density,source_mass_relative=initial_mass_error,
        source_raw_absolute_error=initial_raw_error,numerical_affine_flow=dict(b=b,B=B,pressure=P,diagnostic=affine),
        closure=closure,rows=rows,seconds=time.monotonic()-started,
        density_goal=entry['density_goal'],cell_mass_goal=entry['cell_mass_goal'],
        note='Independent retained-density affine-contact benchmark. Query times are capacity solves, not RK steps. This is a different state representation from the evolving six-moment closure.')
    atomic('curved_transport_report.json',report)
    print('transport complete',entry['n'],report['numerical_passed'],report['physical_accuracy_passed'],flush=True)
    raise SystemExit(0 if report['numerical_passed'] else 2)


if __name__=='__main__':run()
