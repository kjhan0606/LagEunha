"""Independent prescribed-path mass audit on actual C 3D power cells.

The audit uses the derivative of the recorded quadratic midpoint path, never
another mass-root solve. Nested Simpson estimates are numerical checks, not a
continuum enclosure or full acoustic-run certificate.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from check_3d_coupled_rhs import initial
from kinetic_moment_rhs_3d import rhs
from moment_evolution_3d import step,dense_state


def run(library):
    geo=Geometry3D(library);rng=np.random.default_rng(271)
    sites=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])+rng.uniform(-.035,.035,(8,3));weights=np.zeros(8)
    def velocity(p):return .02*np.column_stack([np.sin(np.pi*p[:,0])*np.cos(np.pi*p[:,1]),np.sin(np.pi*p[:,1])*np.cos(np.pi*p[:,2]),np.sin(np.pi*p[:,2])*np.cos(np.pi*p[:,0])])
    raw,L=initial(geo,sites,weights,velocity,dispersion=.003);mass=raw[:,0]
    dr,dx,dw,d=rhs(geo,sites,weights,raw,L,face_order=3)
    ar,ax,aw,a=rhs(geo,sites,weights,raw,L,face_order=3,audit_geometry_rates=(dx,dw))
    np.testing.assert_array_equal(dr,ar);np.testing.assert_array_equal(dw,aw);assert not a['mass_root_solved']
    perturbed=dw+np.linspace(-.01,.01,8);br,bx,bw,b=rhs(geo,sites,weights,raw,L,face_order=3,audit_geometry_rates=(dx,perturbed))
    np.testing.assert_array_equal(bw,perturbed);assert b['mass_flux']>1e-5;np.testing.assert_array_equal(b['mass_rates'],br[:,0])
    rows=[];bounds=lambda t:np.array([[0.,0.,0.],[1.,1.,1.]])
    for dt in (.01,.005,.0025):
        end=step(geo,0.,sites,weights,raw,L,dt,bounds,lambda p:np.zeros(3),return_record=True);record=end[-1]['dense_record'];dense=dense_state(record,1.)
        for actual,expected in zip(dense[:4],end[:4]):np.testing.assert_allclose(actual,expected,atol=3e-14,rtol=0)
        rates=[];maximum_generator=0.;start=time.monotonic()
        for s in np.linspace(0.,1.,17):
            x,w,r,F,v,z=dense_state(record,float(s));np.testing.assert_array_equal(r[:,0],mass)
            _,_,_,audit=rhs(geo,x,w,r,F,bounds(s*dt),lambda p:np.zeros(3),face_order=3,audit_geometry_rates=(v,z))
            assert not audit['mass_root_solved'];rates.append(audit['mass_rates']);maximum_generator=max(maximum_generator,audit['generator_velocity_defect'])
        rates=np.asarray(rates)
        def integrate(values):return dt/(3*(len(values)-1))*(values[0]+values[-1]+4*values[1:-1:2].sum(axis=0)+2*values[2:-1:2].sum(axis=0))
        fine=integrate(rates);coarse=integrate(rates[::2]);absolute=integrate(abs(rates))
        row=dict(dt=dt,stage_mass_flux=end[-1]['stage_mass_flux'],integrated_relative_defect=float(np.max(abs(fine)/mass)),integrated_absolute_relative_defect=float(np.max(absolute/mass)),nested_quadrature_difference=float(np.max(abs(fine-coarse)/mass)),maximum_relative_mass_rate=float(np.max(abs(rates)/mass)),generator_path_defect=maximum_generator,seconds=time.monotonic()-start)
        rows.append(row);Path('progress.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
    report=dict(completed=True,diagnostic_rates_preserved=True,perturbed_rate_residual=b['mass_flux'],rows=rows,full_acoustic_path_certified=False)
    Path('dense_mass_path_report.json').write_text(json.dumps(report,indent=2));return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');a=p.parse_args();print(json.dumps(run(a.library)),flush=True)
