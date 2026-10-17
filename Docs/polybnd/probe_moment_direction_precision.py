"""Optional numerical angle polish: initial common IC and rejected path nodes."""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from common_interface_contact import initialize
from predictive_evolution_batch import density_options,initialize_case,path_sample
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
import kinetic_factor_rhs as kinetic


def run():
    config=json.loads(Path('precision_config.json').read_text());rows=[]
    for angle in (0.,23.):
        mass,y,_,_=initialize(8,100.,angle);L=[np.zeros((2,2)) for _ in mass]
        options=dict(mode='moment',segments=64,radial=1e-7,core_level=2,align_core=True,
            vector_geometry=True,construct_prior_witness=True,lp_retry=True,smooth_phase_floor=True,phase_fraction_floor=1e-7)
        for refine in (False,True):
            started=time.monotonic();provider=BoundaryDensity(**options,refine_moment_direction=refine)
            dy,d=kinetic.rhs(0.,y,mass,covariance_factors=L,face_density_provider=provider,face_segments=64,face_order=4,exact_bin_transport=True)
            target=np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle))])
            angle_error=max((float(np.degrees(np.arctan2(np.cross(target,p['normal']),target@p['normal']))) for p in provider.last_priors if not p['flat']),key=abs)
            row=dict(kind='initial',angle=angle,refine=refine,metrics=metrics(dy,d,mass),max_signed_normal_error_degrees=angle_error,priors=provider.last_priors,seconds=time.monotonic()-started)
            rows.append(row);atomic('precision_progress.json',dict(rows=rows));print('precision initial',angle,refine,angle_error,flush=True)
    root=Path(config['source']);cp=json.loads((root/'inputs/n8_angle0/checkpoint.json').read_text());mass,_=initialize_case(cp['config'])
    data=json.loads((root/'evolution/n8_angle0/short/rejection_5.json').read_text())
    for node in (0,32,64):
        for refine in (False,True):
            options=density_options(cp['config'],True);options['refine_moment_direction']=refine
            value=path_sample(data['record'],mass,node,BoundaryDensity(**options))
            row=dict(kind='rejected_path',index=node,refine=refine,normalized_mass_rate=np.asarray(value['mass_rate'])/mass)
            rows.append(row);atomic('precision_progress.json',dict(rows=rows))
            print('precision path',node,refine,float(np.max(abs(row['normalized_mass_rate']))),flush=True)
    passed=all(r['metrics']['closure_moment_error']<1e-8 and r['metrics']['mass_flux']<1e-10 and max(r['metrics']['affine_rhs_errors'].values())<1e-8 for r in rows if r['kind']=='initial')
    atomic('precision_report.json',dict(completed=True,passed=passed,rows=rows,
        note='Optional least-squares residual polish in the same selected angular basin. Exact normal only scores the initial fit. Stored rejected path rates are not re-solved. This diagnostic has not changed active BK/BJ evolution.'))
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()
