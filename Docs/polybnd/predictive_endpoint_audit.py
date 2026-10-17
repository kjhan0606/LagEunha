"""Read-only density accuracy after accepted predictive evolution.
Exact material polygons are constructed only for scoring reconstructions.
"""
import json,time
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup,diagnostic
from shear_contact_check import affine_errors
from shear_interface_audit import exact_on_cells
from predictive_boundary_density import BoundaryDensity
from predictive_resolution_audit import trace_error
import moment_ale_evolve as model


def run():
    previous=json.loads(Path('evolution_report.json').read_text())
    if not previous['passed']:raise RuntimeError('evolution not complete/passed')
    results={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        source=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=source['config']
        mass=np.asarray(source['mass']);_,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);a0=y0[:-4].reshape(-1,20)
        initial,_=model.geometry(a0[:,:2],a0[:,2],0.);phases=np.where(a0[:,0]<.5,1.,cfg['amplitude']);rows=[]
        short=json.loads((Path('endpoints')/name/'short.json').read_text())
        final=json.loads((Path('endpoints')/name/'extended.json').read_text())
        duration=final['t']-source['t']
        absolute=np.asarray(short['absolute'])+final['absolute'];signed=np.asarray(short['signed'])+final['signed'];quad=np.asarray(short['quad_delta'])+final['quad_delta']
        peak_rate=max(cp['peak']['rate']*duration/cp['config']['duration'] for cp in (short,final))
        corrected=np.asarray(final['state']);corrected[-4:]-=np.asarray(source['state'])[-4:]
        final_diagnostic=diagnostic(final['t'],corrected,mass,np.asarray(source['state']))
        affine=affine_errors(np.asarray(final['state']),mass)
        combined=dict(duration=duration,normalized_peak_mass_rate=peak_rate,
            integrated_absolute_mass_defect=float(np.max(absolute/mass)),integrated_signed_mass_defect=float(np.max(abs(signed)/mass)),
            quadrature_mass_difference=float(np.max(abs(quad)/mass)),ledger_error=final_diagnostic['ledger_error'],affine_errors=affine,
            pressure_error=final_diagnostic['pressure_departure'],generator_defect=max(cp['peak']['generator'] for cp in (short,final)))
        combined['passed']=bool(peak_rate<1e-8 and combined['integrated_absolute_mass_defect']<1e-10 and combined['quadrature_mass_difference']<1e-11 and combined['ledger_error']<1e-10 and max(affine.values())<1e-8 and combined['pressure_error']<1e-8 and combined['generator_defect']<1e-8)
        for interval in ('original','short','extended'):
            cp=source if interval=='original' else json.loads((Path('endpoints')/name/(interval+'.json')).read_text())
            np.testing.assert_allclose(cp['mass'],mass,rtol=0,atol=0)
            t=cp['t'];a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],t)
            raw=np.column_stack([mass,a[:,3:]]);exact,coverage,shape=exact_on_cells(cells,initial,phases,t)
            bound=float(np.sum(abs(mass-exact[:,0]))/mass.sum());global_mass=float(abs(exact[:,0].sum()-mass.sum())/mass.sum())
            refs=[]
            for v in initial:
                r=v.copy();r[:,0]+=t*(.2+1.1*r[:,1]);refs.append(r)
            high_phase=unary_union([Polygon(r) for r,rho in zip(refs,phases) if rho==cfg['amplitude']])
            for mode in ('neighbor','moment'):
                started=time.monotonic();provider=BoundaryDensity(mode,vector_geometry=True,smooth_phase_floor=True)
                _,diag=provider(cells,edges,raw)
                error=0.
                for parts,densities in provider.last_cells:
                    for triangle,rho in zip(parts,densities):
                        p=Polygon(triangle);A=p.area;high_area=p.intersection(high_phase).area
                        error+=abs(rho-1.)*max(0.,A-high_area)+abs(rho-cfg['amplitude'])*high_area
                volume=error/float(mass.sum());face=trace_error(provider,edges,refs,phases)
                valid=coverage<1e-11 and global_mass<1e-11 and volume+1e-8>=bound and diag['closure_moment_error']<1e-8 and diag['closure_min_density']>0 and diag['closure_volume_cap_excess']<1e-8
                row=dict(interval=interval,t=t,mode=mode,smooth_phase_floor=True,reference_valid=bool(valid),density_L1_relative=volume,
                    cell_average_L1_lower_bound=bound,face_pointwise_L1_relative=face['relative_L1'],
                    q_L1_over_total_mass=float(np.sum(abs(raw[:,1:3]-exact[:,1:3]))/mass.sum()),
                    Q_L1_over_total_mass=float(np.sum(abs(raw[:,3:7]-exact[:,3:7]))/mass.sum()),
                    max_partition_area_error=coverage,exact_global_mass_relative_error=global_mass,
                    max_cell_material_symmetric_difference_fraction=shape,closure=diag,inferred_priors=provider.last_priors,seconds=time.monotonic()-started)
                rows.append(row);atomic(name+'_progress.json',dict(rows=rows));print(name,interval,mode,'L1',volume,'face',face['relative_L1'],flush=True)
        results[name]=dict(rows=rows,combined_evolution=combined)
    report=dict(completed=True,passed=all(c['combined_evolution']['passed'] and all(r['reference_valid'] for r in c['rows']) for c in results.values()),cases=results,
        note='Reference-validity gate is not a density-accuracy gate. Exact sheared initial phases used only after numerical reconstruction. Fixed mass; original and continued checkpoints. Moment prior at continued endpoints is diagnostic only; evolution uses neighbor prior. Initial interface differs across n, so no spatial convergence order.')
    atomic('endpoint_density_report.json',report)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
