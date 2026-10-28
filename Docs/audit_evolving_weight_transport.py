"""Sample a geometry-driven weight target and its associated remap cost.

Prescribed circular generator trajectories only. This is a candidate selection
audit, not a fluid solver; donor mixing is a first-order diagnostic, not a
claim about every possible higher-order transport discretization.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import lsmr


def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args()
    args.output.mkdir(exist_ok=False)
    code=Path('/gpfs/kjhan/LagForce/Codex_weighted_compression_geometry_20261017/audit.py')
    spec=importlib.util.spec_from_file_location('weighted_geometry',code)
    geometry_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(geometry_module)
    geom=geometry_module.geometry
    source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator/rhs.dat')
    a=np.loadtxt(source);a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;v=a[:,3:5];n=len(x);r=np.linalg.norm(x,axis=1)
    ids=np.flatnonzero(r<.7);unknown=np.flatnonzero(r<.5);inner=r[ids]<.5
    check=set(ids[::max(1,len(ids)//16)].tolist())
    omega=(x[:,0]*v[:,1]-x[:,1]*v[:,0])/r**2;omega[r>=2.5]=0
    assert np.max(abs(np.sum(x*v,axis=1)[r<.7]))<1e-12
    rho=.01;pressure=2.7777796321446218e-5;gamma=1.6666669845581055
    records=[];old_weight=None;old_time=None;started=time.monotonic()
    for stamp in (0.,.001,.002,.005,.01,.02,.05,.1,.25,.5,1.):
        c,s=np.cos(omega*stamp),np.sin(omega*stamp)
        pos=np.column_stack((c*x[:,0]-s*x[:,1],s*x[:,0]+c*x[:,1]))
        vel=np.column_stack((-omega*pos[:,1],omega*pos[:,0]))
        zero=np.zeros(n);V,D,L,_=geom(pos,zero,vel,ids)
        epsilon=1e-5
        _,_,Lp,_=geom(pos+epsilon*vel,zero,vel,ids)
        _,_,Lm,_=geom(pos-epsilon*vel,zero,vel,ids)
        J=(Lp-Lm)/(2*epsilon)
        scaled=diags(1/V)@J[:,unknown]@diags(a[unknown,6])
        sol=lsmr(scaled,-D/V,damp=1.,atol=1e-10,btol=1e-10,maxiter=2000)
        assert sol[1] in (1,2),('unconverged target',stamp,sol[1])
        weight=np.zeros(n);weight[unknown]=a[unknown,6]*sol[0]
        weight*=min(1.,.25/max(np.max(abs(sol[0])),1e-30))
        volume,rate,lap,certificate=geom(pos,weight,vel,ids,check)
        rms=lambda z:float(np.sqrt(np.average(z[inner]**2,weights=a[ids[inner],6])))
        record={'elapsed':stamp,'lsmr_iterations':int(sol[2]),'max_weight_over_source_volume':float(np.max(abs(weight[unknown]/a[unknown,6]))),
                'voronoi_generator_compression_rms':rms(D/V),'weighted_generator_compression_rms':rms(rate/volume),
                'weighted_to_voronoi_compression_ratio':rms(rate/volume)/rms(D/V),
                'min_volume_over_voronoi':float(np.min(volume/V)),'max_volume_over_voronoi':float(np.max(volume/V)),
                'sample_all_site_plane_violation':certificate}
        if old_weight is not None:
            dt=stamp-old_time
            previous_volume,_,_,_=geom(pos,old_weight,vel,ids)
            change=rho*(volume-previous_volume)
            # At identical positions, the patch perimeter is unaffected by
            # the supported central weights; signed mass changes must close.
            assert abs(change.sum())<1e-12
            weight_rate=(weight-old_weight)/dt
            coo=lap.tocoo();i=ids[coo.row];j=coo.col
            take=i<j;i=i[take];j=j[take];conductance=-coo.data[take]
            assert np.all(conductance>0)
            flux=rho*conductance*(weight_rate[i]-weight_rate[j]) # positive: j -> i
            support=abs(flux)>1e-30
            assert np.all((r[i[support]]<.7)&(r[j[support]]<.7))
            mdot=np.bincount(i,weights=flux,minlength=n)-np.bincount(j,weights=flux,minlength=n)
            assert np.max(abs(mdot[ids]-rho*(lap@weight_rate)))<1e-12
            donor=np.where(flux>=0,j,i)
            pdot=np.zeros((n,2));np.add.at(pdot,i,flux[:,None]*vel[donor]);np.add.at(pdot,j,-flux[:,None]*vel[donor])
            specific_energy=pressure/((gamma-1)*rho)+.5*np.sum(vel*vel,axis=1)
            edot=np.bincount(i,weights=flux*specific_energy[donor],minlength=n)-np.bincount(j,weights=flux*specific_energy[donor],minlength=n)
            kinetic_rate=float(np.sum(vel*pdot)-.5*np.dot(np.sum(vel*vel,axis=1),mdot))
            mixing=.5*abs(flux)*np.sum((vel[i]-vel[j])**2,axis=1)
            assert abs(kinetic_rate+mixing.sum())<1e-12
            inner_fraction=.5*((r[i]<.5).astype(float)+(r[j]<.5).astype(float))
            Uinner=pressure/(gamma-1)*volume[inner].sum()
            record['weight_exchange']={'interval':dt,'uniform_density':rho,
                'exact_frozen_position_mass_change_sum':float(change.sum()),
                'exact_mass_change_fraction_rms':rms(change/(rho*previous_volume)),
                'exact_mass_change_fraction_max':float(np.max(abs(change/(rho*previous_volume)))),
                'endpoint_secant_fractional_mass_rate_rms':rms(mdot[ids]/(rho*volume)),
                'endpoint_secant_absolute_face_mass_flux_sum':float(abs(flux).sum()),
                'global_mass_rate':float(mdot.sum()),'global_momentum_rate':pdot.sum(axis=0).tolist(),
                'global_hydrodynamic_energy_rate':float(edot.sum()),
                'first_order_donor_mixing_power':float(mixing.sum()),
                'first_order_inner_mixing_over_thermal_energy_per_time':float(np.dot(inner_fraction,mixing)/Uinner),
                'warning':'Secant weight rate sampled at endpoint geometry; not a time-integrated remap. Uniform density and point-sampled donor velocities, not evolved fluid. Global energy excludes external gravitational potential.'}
        np.savez(args.output/f'weights_{stamp:g}.npz',weights=weight,ids=ids,volume=volume,generator_volume_rate=rate)
        records.append(record);old_weight=weight;old_time=stamp
        print(json.dumps(record),flush=True)
    out={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'geometry_source_sha256':hashlib.sha256(code.read_bytes()).hexdigest(),
         'scope':'Prescribed circular orbit samples of an independently recomputed regularized weight target. This is not continuous fluid evolution, a smooth weight ODE, or long-term stability evidence.',
         'candidate':'At each sample solve damped least squares of the W=0 generator-compression Jacobian, damping=1; bound max |w|/source_volume by .25. Locations and velocities are not optimized.',
         'mass_coupling':'Only weight-driven face displacement produces remap mass transfer. Equal fixed weights give zero such exchange. Baseline hydrodynamic motion is not simulated.',
         'records':records,'seconds':time.monotonic()-started,
         'limitations':['Target can be nonsmooth across topology events; small-interval rates must be examined before choosing a weight law.',
                        'First-order donor mixing is not a lower bound on high-order conservative reconstruction.',
                        'Pressure, gravity, density contrast, angular momentum moments and coupled energy evolution remain unverified.']}
    (args.output/'result.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
