"""Physical planar IC deliberately inside the finite phase blend band."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
from common_interface_contact import initialize
from discontinuous_density_audit import clip
from predictive_boundary_density import BoundaryDensity
from predictive_evolution_batch import density_options
from predictive_coupling_batch import metrics
from continue_common_contact import score,valid_metric
from neighbor_density_audit import phase_floor_blend
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def choose_angle(n,fraction):
    _,y=setup('shear_contact',n,100.);a=y[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],0.);center=np.array([.5,.5])
    for i,vertices in enumerate(cells):
        area=Polygon(vertices).area
        for vertex in vertices:
            if np.min(vertex)<=1e-4 or np.max(vertex)>=1-1e-4:continue
            delta=vertex-center
            if np.linalg.norm(delta)<1e-10:continue
            critical=float(np.degrees(np.arctan2(-delta[0],delta[1])))
            normal=np.array([np.cos(np.radians(critical)),np.sin(np.radians(critical))])
            values=(vertices-center)@normal
            if not (values.min()>=-1e-12 or values.max()<=1e-12):continue
            def minority(angle):
                theta=np.radians(angle);normal=np.array([np.cos(theta),np.sin(theta)])
                high=clip(vertices,normal,normal@center).area/area
                return min(high,1-high)
            for sign in (-1,1):
                if minority(critical+sign)<fraction:continue
                lo,hi=0.,1.
                for _ in range(45):
                    mid=(lo+hi)/2
                    if minority(critical+sign*mid)<fraction:lo=mid
                    else:hi=mid
                angle=critical+sign*(lo+hi)/2
                if abs(minority(angle)-fraction)>1e-12:continue
                return angle,i,critical
    raise RuntimeError('no physical near-corner phase configuration found')


def run():
    n=3;fraction=2.5e-7;angle,cell,critical=choose_angle(n,fraction)
    mass,y,_,_=initialize(n,100.,angle);a=y[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],0.)
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=means.min();high=means.max()
    actual=float(min(means[cell]-low,high-means[cell])/(high-low))
    if not 1e-7<actual<4e-7:raise RuntimeError('target cell is outside the actual inferred blend band')
    selected=dict(n=n,angle=angle,cell=cell,critical_angle=critical,
        minority_fraction=actual,blend=phase_floor_blend(means[cell],low,high,1e-7))
    atomic('selected_physical_initial_condition.json',selected);results=[]
    for epsilon in (0.,1e-8,1e-7,1e-6):
        cfg=dict(kind='common_planar_contact',n=n,amplitude=100.,angle=angle,smooth_phase_floor=True,
            path_control_horizon=1e-4,density_provider=dict(mode='moment',segments=64,radial=1e-7,
            core_level=2,align_core=True,construct_prior_witness=True,lp_retry=True,
            phase_fraction_floor=epsilon,refine_moment_direction=False))
        cp=dict(config=cfg,t=0.,mass=mass,state=y,factors=[np.zeros((2,2)) for _ in mass])
        tag='epsilon_'+format(epsilon,'.0e');atomic(tag+'_initial.json',cp)
        try:
            provider=BoundaryDensity(**density_options(cfg,True))
            dy,diag=kinetic.rhs(0.,y,mass,covariance_factors=cp['factors'],face_density_provider=provider,
                face_segments=64,face_order=4,exact_bin_transport=True)
            m=metrics(dy,diag,mass);row=score(cp,cfg,mass,tag)
            result=dict(epsilon=epsilon,passed=bool(valid_metric(m) and row['reference_valid']),
                metrics=m,score=row,target_prior=provider.last_priors[cell],
                weight_rate=np.asarray(dy)[:-4].reshape(-1,20)[:,2])
        except Exception as exc:result=dict(epsilon=epsilon,passed=False,error=repr(exc))
        results.append(result);atomic('phase_band_progress.json',dict(selected=selected,results=results))
        print('phase band',epsilon,result['passed'],result.get('error',''),flush=True)
    report=dict(completed=True,passed=next(r for r in results if r['epsilon']==1e-7)['passed'],
        all_variants_passed=all(r['passed'] for r in results),selected=selected,results=results,
        note='Exact moments of a new physical planar IC place one cell inside the actual finite blend band. This is a static RHS/reconstruction probe, not a temporal band-crossing or epsilon->0 convergence test. References only initialize and score; each variant solves its own instantaneous weight rate.')
    atomic('phase_band_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
