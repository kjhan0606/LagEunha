"""Reconstruction-independent density L1 lower bounds from six stored moments.

A quadratic polynomial is bounded on a triangle by its Bernstein coefficients.
The LP constructs a continuous bounded test function, not a density fit to the
exact solution. Reference fields enter this diagnostic only.
"""
import json,os
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from adaptive_coupled_mass import atomic
from simple_wave_exact import physical_raw
from simple_wave_evolution import install_domain
import moment_ale_evolve as model


def reference_spec(report):
    problem=report.get('problem','nonlinear_simple_wave')
    if problem=='nonlinear_simple_wave':
        return physical_raw,report['amplitude'],report['config']['physical_goal'],16,20
    if problem=='steady_vortex':
        from steady_vortex_exact import physical_raw as vortex_raw
        return vortex_raw,1.,report['config']['physical_goals']['density_field'],16,24
    raise ValueError('unknown density reference problem: '+problem)


def basis(points,origin,scale):
    z=(np.asarray(points,dtype=np.longdouble)-np.asarray(origin,dtype=np.longdouble))/np.asarray(scale,dtype=np.longdouble)
    return np.column_stack([np.ones(len(z)),z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])


def bernstein_matrix(poly,origin,scale):
    poly=np.asarray(poly,dtype=np.longdouble);center=poly.mean(0);rows=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        tri=np.array([center,a,b]);values=basis(tri,origin,scale);rows.extend(values)
        for i,j in [(0,1),(0,2),(1,2)]:rows.append(2*basis([(tri[i]+tri[j])/2],origin,scale)[0]-.5*(values[i]+values[j]))
    return np.asarray(rows,dtype=np.longdouble)


def moment_vector(raw,origin,scale):
    raw=np.asarray(raw,dtype=np.longdouble);o=np.asarray(origin,dtype=np.longdouble);h=np.asarray(scale,dtype=np.longdouble)
    m=raw[0];q=raw[1:3];Q=raw[3:7].reshape(2,2)
    first=q-o*m;second=Q-np.outer(o,q)-np.outer(q,o)+m*np.outer(o,o)
    return np.r_[m,first/h,second[0,0]/h[0]**2,(second[0,1]+second[1,0])/(2*h[0]*h[1]),second[1,1]/h[1]**2]


def dual_bound(poly,discrepancy,origin,scale):
    B=bernstein_matrix(poly,origin,scale);delta=np.asarray(discrepancy,dtype=np.longdouble);size=float(np.max(abs(delta)))
    if size==0:return 0.,np.zeros(6)
    result=linprog(-np.asarray(delta/size,dtype=float),A_ub=np.asarray(np.vstack([B,-B]),dtype=float),b_ub=np.ones(2*len(B)),bounds=[(None,None)]*6,method='highs',options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    if not result.success:raise RuntimeError('bounded polynomial dual failed: '+result.message)
    coefficient=np.asarray(result.x,dtype=np.longdouble)
    # Verify the certificate with extended arithmetic. This rescales the test
    # function only; the numerical density and its moments are never changed.
    envelope=max(np.longdouble(1),np.max(abs(B@coefficient)))+64*np.finfo(np.longdouble).eps*np.max(np.sum(abs(B*coefficient),axis=1))
    coefficient/=envelope
    return float(abs(coefficient@delta)),coefficient


def cellwise_supremum_lower_bound(l1_bound,poly,density_scale):
    p=np.asarray(poly,dtype=np.longdouble);p=p-p.mean(0);q=np.roll(p,-1,axis=0)
    area=.5*abs(np.sum(p[:,0]*q[:,1]-p[:,1]*q[:,0]))
    if area<=0 or density_scale<=0:raise ValueError('invalid local density bound scale')
    return float(l1_bound/(area*density_scale))


def inspect(entry):
    root=Path(entry['root']);cp=json.loads((root/'checkpoint.json').read_text());report=json.loads((root/'simple_wave_report.json').read_text());mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    if cp['t']!=report['t'] or cp['t']!=report['config']['end']:raise ValueError('density bound needs the reported endpoint')
    from physical_field_accuracy_audit import environment
    os.environ.update(environment(report['config']))
    raw_function,scale_density,goal,low_order,high_order=reference_spec(report)
    cells,_=model.geometry(a[:,:2],a[:,2],cp['t']);total=quad=0.;rows=[]
    for i,(poly,m,row) in enumerate(zip(cells,mass,a)):
        lo=np.min(poly,axis=0);hi=np.max(poly,axis=0);origin=(lo+hi)/2;scale=(hi-lo)/2
        ref,*_=raw_function(poly,cp['t'],report['amplitude'],high_order);coarse,*_=raw_function(poly,cp['t'],report['amplitude'],low_order)
        numerical=np.r_[m,row[3:9]];reference=moment_vector(ref,origin,scale);delta=moment_vector(numerical,origin,scale)-reference
        value,coefficient=dual_bound(poly,delta,origin,scale);gap=float(abs(coefficient@(reference-moment_vector(coarse,origin,scale))))
        total+=value;quad+=gap;rows.append(dict(cell=i,bound=value,reference_difference=gap,cellwise_supremum_lower_bound_over_scale=cellwise_supremum_lower_bound(max(0.,value-gap),poly,scale_density),coefficient=np.asarray(coefficient,float),moment_discrepancy=np.asarray(delta,float)))
    result=dict(label=entry['label'],problem=report.get('problem','nonlinear_simple_wave'),n=report['n'],t=cp['t'],density_scale=scale_density,physical_density_goal=goal,density_moment_lower_bound_over_scale=total/scale_density,
        reference_witness_difference_over_scale=quad/scale_density,density_field_error_over_scale=report['final_score']['errors_over_amplitude']['density_field'],reference_orders=[low_order,high_order],
        lower_bound_exceeds_physical_goal=bool((total-quad)/scale_density>goal),worst_cells=sorted(rows,key=lambda v:v['bound'],reverse=True)[:12],
        note='Any integrable density matching the stored m,q,Q obeys this L1 lower bound, conditional on reference moments. Polynomial bounds verified through Bernstein coefficients on all polygon fan triangles. Reference quadrature difference is reported; no rigorous analytic quadrature enclosure is claimed. No numerical moment or evolution is modified.')
    if result['problem']=='nonlinear_simple_wave':
        result.update(density_moment_lower_bound_over_amplitude=total/scale_density,reference_witness_difference_over_amplitude=quad/scale_density,density_field_error_over_amplitude=result['density_field_error_over_scale'])
    result['cellwise_supremum_lower_bound_over_scale']=max(v['cellwise_supremum_lower_bound_over_scale'] for v in rows)
    result['worst_local_cells']=sorted(rows,key=lambda v:v['cellwise_supremum_lower_bound_over_scale'],reverse=True)[:12]
    result['local_bound_note']='L_infinity on a cell is at least its certified L1 moment bound divided by cell area. Reference quadrature difference is subtracted as an observed screen, not a rigorous analytic quadrature enclosure. No new pointwise acceptance target is imposed.'
    atomic(entry['label']+'_obstruction.json',result);print(result,flush=True);return result


def run():
    install_domain();cfg=json.loads(Path('moment_obstruction_config.json').read_text());rows=[inspect(v) for v in cfg['cases']]
    atomic('moment_obstruction_report.json',dict(completed=True,rows=rows))


if __name__=='__main__':run()
