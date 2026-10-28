"""Reconstruct full weighted-minus-Voronoi face-speed difference."""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351

def relative_weight_speed(delta_w,delta_wdot,d,ddot):
    """Normal power-face speed relative to the Voronoi face at same sites."""
    return delta_wdot/(2*d)-delta_w*ddot/(2*d*d)

def verify_derivative():
    rng=np.random.default_rng(20261027)
    d=rng.uniform(.02,.2,5000);ddot=rng.uniform(-.8,.8,5000)
    dw=rng.uniform(-.0002,.0002,5000);dwdot=rng.uniform(-.01,.01,5000)
    h=1e-7
    finite=((dw+h*dwdot)/(2*(d+h*ddot))-dw/(2*d))/h
    analytic=relative_weight_speed(dw,dwdot,d,ddot)
    scale=np.maximum(1.,np.abs(analytic))
    error=float(np.max(np.abs(finite-analytic)/scale))
    if error>2e-5:raise ValueError(f'Power-plane derivative check failed: {error}')
    return error

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('receipt',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    audit=json.loads(args.receipt.read_text())
    if len(audit['checks'])!=21 or not all(audit['checks'].values()):
        raise ValueError('Native weighted step did not pass')
    for key in ('initial','parent','endpoint'):
        if digest(Path(audit[key]))!=audit[key+'_sha256']:raise ValueError(f'Changed {key}')
    initial,_=ordered(Path(audit['initial']),N)
    p0,m0=ordered(Path(audit['parent']),N);p1,m1=ordered(Path(audit['endpoint']),N)
    dt=float(m1['t'])-float(m0['t'])
    if abs(dt-.001)>1e-10 or not np.array_equal(p0['mass'],p1['mass']):
        raise ValueError('Expected accepted fixed-mass native step')
    xy0=np.column_stack((p0['x'].astype(float),p0['y'].astype(float)))
    xy1=np.column_stack((p1['x'].astype(float),p1['y'].astype(float)))
    xy=.5*(xy0+xy1);vi=.5*np.column_stack((p0['vx'].astype(float)+p1['vx'].astype(float),
                                           p0['vy'].astype(float)+p1['vy'].astype(float)))
    w0=p0['w2'].astype(float);w1=p1['w2'].astype(float)
    w=.5*(w0+w1);wdot=(w1-w0)/dt
    r0=np.hypot(initial['x'].astype(float)-3.,initial['y'].astype(float)-3.)
    active=np.flatnonzero((r0>=.5)&(r0<2.));flag=np.zeros(N,bool);flag[active]=True
    lookup=np.full(N,-1,dtype=int);lookup[active]=np.arange(len(active))
    area,meta,tri,circ=power_areas(xy,w,active,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid midpoint power cells')
    i,j,face,length=finite_faces_touching_active(tri,circ,flag)
    delta=xy[j]-xy[i];d=np.linalg.norm(delta,axis=1);n=delta/d[:,None]
    dv=vi[j]-vi[i];ddot=np.einsum('ij,ij->i',n,dv)
    delta_w=w[i]-w[j];delta_wdot=wdot[i]-wdot[j]
    speed_rate=delta_wdot/(2*d)
    speed_separation=-delta_w*ddot/(2*d*d)
    speed_full=speed_rate+speed_separation
    secant_d=(np.linalg.norm(xy1[j]-xy1[i],axis=1)-np.linalg.norm(xy0[j]-xy0[i],axis=1))/dt
    speed_secant=delta_wdot/(2*d)-delta_w*secant_d/(2*d*d)
    vdot_error=float(np.sqrt(np.mean((ddot-secant_d)**2)))
    rho=.5*(p0['den'].astype(float)+p1['den'].astype(float));mass=p0['mass'].astype(float)
    vol0=p0['volume'].astype(float);vol1=p1['volume'].astype(float)
    disk_mass=mass[active]
    def sweep(speed,donor):
        rf=.5*(rho[i]+rho[j]) if donor=='arithmetic' else np.where(speed>=0,rho[j],rho[i])
        fm=dt*length*speed*rf
        dm=np.zeros(N);np.add.at(dm,i,fm);np.add.at(dm,j,-fm)
        gross=float(np.sum(np.abs(fm),dtype=np.longdouble));imbalance=float(np.sum(dm,dtype=np.longdouble))
        if np.any(mass+dm<=0) or abs(imbalance)>1e-12*max(gross,1e-30):
            raise ValueError(f'{donor} transfer fails positivity or conservation')
        change=(mass+dm)/vol1/(mass/vol0)-1.
        return {'disk_mean_abs_density_change':float(np.average(np.abs(change[active]),weights=disk_mass)),
                'gross_face_mass_transfer':gross,
                'disk_mass_fraction_change_p99':float(np.quantile(np.abs(dm[active])/mass[active],.99)),
                'sum_cell_mass_change':imbalance}
    fixed=(mass/vol1)/(mass/vol0)-1.
    result={'initial_sha256':audit['initial_sha256'],'start_sha256':audit['parent_sha256'],
            'end_sha256':audit['endpoint_sha256'],'start_time':float(m0['t']),
            'end_time':float(m1['t']),'dt':dt,'disk_cells':len(active),
            'faces_touching_disk':len(i),'analytic_plane_derivative_test_max_relative_error':verify_derivative(),
            'generator_separation_speed_rms':float(np.sqrt(np.mean(ddot*ddot))),
            'separation_secant_speed_rms':float(np.sqrt(np.mean(secant_d*secant_d))),
            'separation_speed_velocity_vs_secant_rms_difference':vdot_error,
            'weight_rate_term_speed_rms':float(np.sqrt(np.mean(speed_rate*speed_rate))),
            'finite_weight_times_separation_term_speed_rms':float(np.sqrt(np.mean(speed_separation*speed_separation))),
            'full_relative_weight_face_speed_rms':float(np.sqrt(np.mean(speed_full*speed_full))),
            'rate_only_to_full_speed_rms_ratio':float(np.sqrt(np.mean(speed_rate**2))/np.sqrt(np.mean(speed_full**2))),
            'velocity_secant_full_speed_rms_difference':float(np.sqrt(np.mean((speed_full-speed_secant)**2))),
            'fixed_mass_disk_mean_abs_density_change':float(np.average(np.abs(fixed[active]),weights=disk_mass)),
            'rate_only_arithmetic_remap':sweep(speed_rate,'arithmetic'),
            'full_relative_speed_arithmetic_remap':sweep(speed_full,'arithmetic'),
            'full_relative_speed_upwind_remap':sweep(speed_full,'upwind'),
            'scope':'First-order midpoint mass-sweep diagnostic. Includes both weight-rate face motion and finite-weight correction from changing generator separation, relative to the Voronoi face at identical sites. It is not a native RK4 ALE flux or evolved variable-mass solution.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
