"""Analytic potential derivative of the exponential quadratic density closure.

For fixed local polynomial basis b, rho=exp(theta.b), A=int rho b b^T,
lambda=A^-1 int rho Phi b. Reynolds transport gives
Udot=lambda.Mdot + integral_boundary rho (Phi-lambda.b) w_normal.
This is diagnostic only; it does not alter the force, flux, mass or energy.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from adaptive_coupled_mass import restore
from gravity_coupled_2d import hydro_dense_state
from audit_kepler_energy_step import configure,phi,atomic
from kepler_disk_reference import polar_rule


def basis(x,density):
    return np.column_stack([np.ones(len(x)),base.features(x,density.origin,density.scale)])


def local_moment_rate(da,origin,scale):
    ox,oy=origin;qx,qy=da[3:5];Q=da[5:9].reshape(2,2)
    return np.array([0.,qx/scale,qy/scale,(Q[0,0]-2*ox*qx)/scale**2,
        (.5*(Q[0,1]+Q[1,0])-ox*qy-oy*qx)/scale**2,(Q[1,1]-2*oy*qy)/scale**2])


def run(seed,enrichment=False):
    root=Path(seed);names=('source.json','trial.json','configuration.json')
    data={n:(root/n).read_bytes() for n in names};configure(json.loads(data['configuration.json']))
    mass=np.asarray(json.loads(data['source.json'])['mass'])
    record=restore(json.loads(data['trial.json'])['record']);y,L,dy=hydro_dense_state(record,mass,0.)
    a=y[:-4].reshape(-1,20);da=dy[:-4].reshape(-1,20);mask=np.asarray(record['frozen_first'],dtype=bool)
    cells,edges=base.geometry(a[:,:2],a[:,2],record['t'])
    states=[base.decode(p,np.r_[m,row[3:]]) for p,m,row in zip(cells,mass,a)]
    density=[s['density'] for s in states]
    coefficients=[];rate=[];variance=[];potential_mean=[]
    for i,d in enumerate(density):
        x,w=polar_rule(cells[i],16,radii=(.25,.5));dm=w*d.value(x);b=basis(x,d)
        A=b.T@((dm/mass[i])[:,None]*b);rhs=b.T@(dm*phi(x)/mass[i])
        coefficient=np.linalg.solve(A,rhs);coefficients.append(coefficient)
        remainder=phi(x)-b@coefficient
        variance.append(float(dm@(remainder**2)/mass[i]))
        potential_mean.append(float(dm@phi(x)/mass[i]))
        rate.append(float(coefficient@local_moment_rate(da[i],d.origin,d.scale)))
    rate=np.asarray(rate);rows=[]
    for order in (8,16,24):
        z,w=np.polynomial.legendre.leggauss(order);z=(z+1)/2;w=w/2;boundary=np.zeros(len(mass))
        for i,j,start,end,n,length in edges:
            if j<0:continue # fixed reflecting outer walls
            x=start+z[:,None]*(end-start);delta=a[j,:2]-a[i,:2]
            wn=(delta@da[j,:2]-(x-a[i,:2])@(da[j,:2]-da[i,:2])+.5*(da[i,2]-da[j,2]))/np.linalg.norm(delta)
            for k,sign in ((i,1),(j,-1)):
                d=density[k];residual=phi(x)-basis(x,d)@coefficients[k]
                boundary[k]+=sign*length*float(w@(d.value(x)*residual*wn))
        total=rate+boundary
        rows.append(dict(face_order=order,moment_part=float(rate[~mask].sum()),
            moving_boundary_part=float(boundary[~mask].sum()),potential_rate=float(total[~mask].sum())))
    report=dict(identity={n:hashlib.sha256(v).hexdigest() for n,v in data.items()},
        rows=rows,full_validation=False,scope='analytic derivative of reconstructed U on one saved hydro stage')
    if enrichment:
        import os
        from smooth_mass_path import evaluate_smooth_path
        from quadratic_velocity_reconstruction import prepare,velocity_samples
        from kinetic_factor_rhs import ALPHA
        from kepler_spline_gravity import acceleration
        from kepler_disk_reference import CENTER
        if os.environ.get('LAGFORCE_QUADRATIC_VELOCITY')=='1':prepare(states,cells,edges,mass,L)
        work=np.zeros(len(mass))
        for i,s in enumerate(states):
            d=density[i];x=d.points;v,_=velocity_samples(s,L[i],mass[i],x,ALPHA)
            work[i]=d.weighted_mass@np.sum(v.mean(1)*acceleration(x,center=CENTER,rc=.5),axis=1)
        _,diag=evaluate_smooth_path(record['t'],y,L,mass,(da[:,:2],da[:,2]),frozen_mask=mask,potential=phi)
        missing=diag['potential_mass_rates']-work-total
        sigma=np.sqrt(np.asarray(variance));active=np.flatnonzero(~mask)
        if np.any(sigma[active]<=0):raise ValueError('potential enrichment has zero independent variance')
        # rho_dot=rho*(Phi-projection(Phi))*missing/(m*variance)
        # Its dimensionless RMS fractional-density rate is |missing|/(m*sigma).
        stiffness=abs(missing)/(mass*sigma)
        worst=int(active[np.argmax(stiffness[active])])
        report['enrichment']=dict(sum_missing_rate=float(missing[active].sum()),
            maximum_rms_fractional_density_rate=float(stiffness[worst]),worst_cell=worst,
            predicted_rms_change_at_dt00025=float(.00025*stiffness[worst]))
        atomic('enrichment_cells.json',dict(missing=missing,variance=variance,mass=mass,frozen=mask,rms_rate=stiffness))
        atomic('worst_cell.json',dict(cell=worst,polygon=cells[worst],raw=np.r_[mass[worst],a[worst,3:]],
            coefficient=coefficients[worst],sigma=sigma[worst],potential_mean=potential_mean[worst],
            missing_rate=missing[worst],density_quadrature_order=density[worst].quadrature_order))
    atomic('report.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('seed');p.add_argument('--enrichment',action='store_true');a=p.parse_args();run(a.seed,a.enrichment)
