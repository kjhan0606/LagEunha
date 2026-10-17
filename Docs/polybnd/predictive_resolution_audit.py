"""Predictive boundary-trace resolution and layer-thickness diagnostics.
Reference density is used only after reconstruction/RHS have been evaluated.
"""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
from predictive_coupling_batch import saved_case,metrics
from predictive_boundary_density import BoundaryDensity
from face_density_selection_audit import one_sided_interval
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic


def references(name,t):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    _,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);a=y0[:-4].reshape(-1,20)
    initial,_=model.geometry(a[:,:2],a[:,2],0.);refs=[]
    for v in initial:
        r=v.copy();r[:,0]+=t*(.2+1.1*r[:,1]);refs.append(r)
    return refs,np.where(a[:,0]<.5,1.,cfg['amplitude'])


def trace_error(provider,edges,refs,phases):
    error=norm=coverage_error=0.
    for edge_id,(i,j,a,b,n,length) in enumerate(edges):
        for cell,normal in ((i,n),(j,-n)):
            if cell<0:continue
            intervals=[];coverage=0.
            for ref,rho in zip(refs,phases):
                interval=one_sided_interval(ref,a,b,-normal)
                if interval is not None:
                    lo,hi=interval;coverage+=hi-lo;intervals.append((lo,hi,rho))
            coverage_error=max(coverage_error,abs(coverage-1))
            profile=provider.last_profiles[(cell,edge_id)];S=len(profile)
            for k,value in enumerate(profile):
                for lo,hi,rho in intervals:
                    overlap=max(0.,min((k+1)/S,hi)-max(k/S,lo))
                    error+=length*overlap*abs(value-rho);norm+=length*overlap*abs(rho)
    if coverage_error>1e-8:raise RuntimeError('reference face coverage failed')
    return dict(relative_L1=error/norm,absolute_L1=error,normalization=norm,coverage_error=coverage_error)


def run():
    if not json.loads(Path('static_report.json').read_text())['passed']:raise RuntimeError('static gates failed')
    all_cases={};passed=True
    for name in ('n3_r10','n3_r100','n4_r100'):
        cp,mass,y,L=saved_case(name);a=y[:-4].reshape(-1,20);_,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
        refs,phases=references(name,cp['t']);rows=[];derivatives={}
        for S,radial in [(s,r) for s in (16,32,64,128) for r in (1e-6,1e-7)]+[(256,1e-7)]:
            tag=f'{name}_n{S}_r{radial:.0e}';path=Path(tag+'.json');started=time.monotonic()
            if path.exists():
                row=json.loads(path.read_text())
            else:
                provider=BoundaryDensity(segments=S,radial=radial)
                try:
                    dy,d=kinetic.rhs(cp['t'],y,mass,covariance_factors=L,face_density_provider=provider,
                        face_segments=S,face_order=4,exact_bin_transport=True)
                    m=metrics(dy,d,mass)
                    gates=m['mass_flux']<1e-10 and m['closure_moment_error']<1e-8 and m['density_min']>0 and m['closure_volume_cap_excess']<1e-8 and m['pressure_rate']<1e-8 and max(m['affine_rhs_errors'].values())<1e-8 and m['ledger_rate_error']<1e-10
                    row=dict(segments=S,radial=radial,status='validated' if gates else 'gate_failure',passed=gates,
                        metrics=m,derivative=dy,face_density_error=trace_error(provider,edges,refs,phases),seconds=time.monotonic()-started)
                except Exception as exc:
                    row=dict(segments=S,radial=radial,status='closure_unavailable',passed=False,error=repr(exc),seconds=time.monotonic()-started)
                atomic(path,row)
            if row['status']=='validated':derivatives[(S,radial)]=np.asarray(row['derivative'])
            rows.append({k:v for k,v in row.items() if k!='derivative'})
            print(tag,row['status'],row.get('face_density_error',{}).get('relative_L1'),flush=True)
        differences=[]
        for S in (16,32,64,128):
            key=(S,1e-7);fine=(2*S,1e-7)
            if key in derivatives and fine in derivatives:
                d=derivatives[key]-derivatives[fine];r=d[:-4].reshape(-1,20)
                differences.append(dict(segments=[S,2*S],max_abs_rhs=float(np.max(abs(d))),max_weight_rate=float(np.max(abs(r[:,2]))),
                    max_raw_moment_rate_over_mass=float(np.max(abs(r[:,3:])/mass[:,None]))))
        thickness=[]
        for S in (16,32,64,128):
            if (S,1e-6) in derivatives and (S,1e-7) in derivatives:
                thickness.append(dict(segments=S,max_abs_rhs=float(np.max(abs(derivatives[(S,1e-6)]-derivatives[(S,1e-7)])))))
        # Coarser unavailable reconstructions are outcomes, not accepted evolution states.
        case_pass=(64,1e-7) in derivatives and (128,1e-7) in derivatives and (256,1e-7) in derivatives
        passed=passed and case_pass
        all_cases[name]=dict(passed=case_pass,rows=rows,successive_rhs_differences=differences,thickness_rhs_differences=thickness)
        atomic('resolution_progress.json',dict(cases=all_cases))
    report=dict(completed=True,passed=bool(passed),cases=all_cases,
        note='Numerical-only inferred trace before reference scoring. Fixed saved t=.01 state. Tangential-bin refinement is not a spatial-grid convergence test. Layer-thickness invariance is expected for the RHS once the same prescribed face profiles are feasible; it does not certify the interior density. No evolution replay.')
    atomic('resolution_report.json',report)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()
