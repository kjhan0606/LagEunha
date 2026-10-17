"""Static common-field checks of moment normals and aligned bounded density."""
import time,json,argparse
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from common_interface_contact import initialize
from predictive_boundary_density import BoundaryDensity
from predictive_coupling_batch import metrics
from predictive_resolution_audit import trace_error
from face_density_selection_audit import one_sided_interval
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def reference_bin_error(provider,edges,refs,phases):
    error=norm=maximum=jump=0.;S=provider.segments;left=np.arange(S)/S;right=(np.arange(S)+1)/S
    for edge_id,(i,j,a,b,n,length) in enumerate(edges):
        for cell,normal in ((i,n),(j,-n)):
            if cell<0:continue
            exact=np.zeros(S);coverage=0.
            for ref,rho in zip(refs,phases):
                interval=one_sided_interval(ref,a,b,-normal)
                if interval is None:continue
                lo,hi=interval;coverage+=hi-lo
                exact+=rho*S*np.maximum(0.,np.minimum(right,hi)-np.maximum(left,lo))
            if abs(coverage-1)>1e-8:raise RuntimeError('common-reference face coverage failed')
            difference=abs(provider.last_profiles[(cell,edge_id)]-exact)
            maximum=max(maximum,float(difference.max()/max(phases)))
            error+=length*float(difference.mean());norm+=length*float(exact.mean())
        if j>=0:jump=max(jump,float(np.max(abs(provider.last_profiles[(i,edge_id)]-provider.last_profiles[(j,edge_id)])))/max(phases))
    return dict(relative_L1=error/norm,max_error_over_high_density=maximum,max_shared_face_jump_over_high_density=jump)


def run():
    parser=argparse.ArgumentParser();parser.add_argument('--reuse');parser.add_argument('--retry-lp',action='store_true');parser.add_argument('--prior-witness',action='store_true');args=parser.parse_args()
    previous=json.loads(Path(args.reuse).read_text())['results'] if args.reuse else []
    reusable={(r['n'],r['angle'],r['segments']):r for r in previous if r['passed']}
    rows=[]
    for n in (3,4,8):
        for angle in (0.,23.):
            mass,y,refs,phases=initialize(n,100.,angle);a=y[:-4].reshape(-1,20);_,edges=model.geometry(a[:,:2],a[:,2],0.)
            L=[np.zeros((2,2)) for _ in mass];derivatives={}
            for S in (16,64,256):
                if (n,angle,S) in reusable:
                    row=dict(reusable[(n,angle,S)],reused=True);rows.append(row)
                    atomic('common_resolution_progress.json',dict(results=rows));print(n,angle,S,'reused',flush=True);continue
                started=time.monotonic();provider=BoundaryDensity('moment',segments=S,vector_geometry=True,align_core=True,lp_retry=args.retry_lp,construct_prior_witness=args.prior_witness)
                try:
                    dy,d=kinetic.rhs(0.,y,mass,covariance_factors=L,face_density_provider=provider,face_segments=S,face_order=2,exact_bin_transport=True)
                    m=metrics(dy,d,mass);face=trace_error(provider,edges,refs,phases);means=reference_bin_error(provider,edges,refs,phases)
                    expected=np.array([np.cos(np.radians(angle)),np.sin(np.radians(angle))])
                    angles=[abs(float(np.degrees(np.arctan2(p['normal'][0]*expected[1]-p['normal'][1]*expected[0],np.dot(p['normal'],expected))))) for p in provider.last_priors if not p['flat']]
                    angle_error=max(angles,default=0.)
                    quadrature_difference=0.
                    if S==64:
                        higher,_=kinetic.rhs(0.,y,mass,covariance_factors=L,face_density_provider=provider,face_segments=S,face_order=4,exact_bin_transport=True)
                        quadrature_difference=float(np.max(abs(dy-higher)/np.maximum(1.,abs(higher))))
                    passed=m['mass_flux']<1e-10 and m['closure_moment_error']<1e-8 and m['density_min']>0 and m['closure_volume_cap_excess']<1e-8 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 and m['ledger_rate_error']<1e-10 and means['max_error_over_high_density']<1e-6 and means['max_shared_face_jump_over_high_density']<1e-6 and angle_error<1e-4
                    passed=passed and quadrature_difference<1e-8
                    derivative_difference=float(np.max(abs(dy-derivatives[S//4]))) if S//4 in derivatives else None
                    derivatives[S]=dy
                    row=dict(n=n,angle=angle,segments=S,passed=bool(passed),status='validated' if passed else 'gate_failure',metrics=m,
                        face_pointwise_L1_relative=face['relative_L1'],bin_means=means,max_normal_angle_error_degrees=angle_error,
                        rhs_difference_from_previous_resolution=derivative_difference,seconds=time.monotonic()-started)
                    row['split_order2_vs4_difference']=quadrature_difference if S==64 else None
                    row['solver_attempts']=[f['solver_attempts'] for f in provider.last_fits if f['solver_attempts']]
                    row['prior_witnesses']=d.get('closure_prior_witnesses',0);row['lp_calls']=d['closure_lp_calls']
                    if args.prior_witness:
                        atomic(f'witness_n{n}_angle{angle:g}_S{S}.json',dict(n=n,angle=angle,segments=S,core_level=2,align_core=True,radial=1e-7,
                            priors=provider.last_priors,densities=[rho for _,rho in provider.last_cells],fits=provider.last_fits,
                            note='Numerical prior plus moment correction. Objective <=1e-5 and bulk optimum >=0 give original near-optimal budget without an LP optimum.'))
                except Exception as exc:row=dict(n=n,angle=angle,segments=S,passed=False,status='failed',error=repr(exc),seconds=time.monotonic()-started)
                rows.append(row);atomic('common_resolution_progress.json',dict(results=rows));print(n,angle,S,row['status'],row.get('face_pointwise_L1_relative'),flush=True)
    report=dict(completed=True,passed=all(r['passed'] for r in rows),results=rows,
        note='Common physical planar initial fields, exact m,q,Q only at initialization. Moment-based normal, inferred-interface aligned core, bounded numerical density. Reference only scores completed reconstruction. Two-point Gauss is exact for cubic affine-velocity/constant-pressure kinetic flux within sign-split bins. Static normal/shared-trace/flux tests, not spatial evolution convergence.')
    atomic('common_resolution_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
