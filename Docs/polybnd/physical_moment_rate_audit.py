"""Read-only consistency of centered density and velocity moments.

All numerical derivatives come from preserved initial RHS files. Exact Euler
Reynolds derivatives are diagnostic references in the same moving cells; they
never enter evolution. Moment centering is differentiated, including dm/dt.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from simple_wave_evolution import install_domain
from simple_wave_rhs_consistency import exact_rates
import moment_ale_evolve as model


def centered_rates(raw,rate):
    r=np.asarray(raw,dtype=np.longdouble);d=np.asarray(rate,dtype=np.longdouble)
    m=r[0];dm=d[0];q=r[1:3];dq=d[1:3];Q=r[3:7].reshape(2,2);dQ=d[3:7].reshape(2,2)
    p=r[7:9];dp=d[7:9];D=r[9:13].reshape(2,2);dD=d[9:13].reshape(2,2)
    Pi=r[13:17].reshape(2,2);dPi=d[13:17].reshape(2,2)
    S=Q-np.outer(q,q)/m;C=D-np.outer(q,p)/m
    dS=dQ-(np.outer(dq,q)+np.outer(q,dq))/m+np.outer(q,q)*dm/m**2
    dC=dD-(np.outer(dq,p)+np.outer(q,dp))/m+np.outer(q,p)*dm/m**2
    X=np.asarray(np.linalg.solve(np.asarray(S,float),np.asarray(C,float)),dtype=np.longdouble)
    B=X.T
    dB=np.linalg.solve(np.asarray(S,float),np.asarray(dC-dS@X,float)).T
    R=Pi-np.outer(p,p)/m-C.T@X
    dR=dPi-(np.outer(dp,p)+np.outer(p,dp))/m+np.outer(p,p)*dm/m**2-dC.T@X-X.T@dC+X.T@dS@X
    return dict(centroid=np.asarray((dq-q*dm/m)/m,float),velocity=np.asarray((dp-p*dm/m)/m,float),
        density_covariance=np.asarray(dS/m-S*dm/m**2,float),velocity_gradient=dB,
        residual_velocity_covariance=np.asarray(dR/m-R*dm/m**2,float))


def inspect(entry):
    initial=json.loads(Path(entry['initial']).read_text());mass=np.asarray(initial['mass']);y=np.asarray(initial['state']);a=y[:-4].reshape(-1,20);t=initial['t']
    source=json.loads(Path(entry['rhs']).read_text())
    if entry.get('trial_record'):
        record=restore(source['record']);np.testing.assert_array_equal(record['y'],y);dy=np.asarray(record['d0'])
    else:
        if source['t']!=t or source['state_sha256']!=hashlib.sha256(y.tobytes()).hexdigest():raise ValueError('initial RHS state identity mismatch')
        dy=np.asarray(source['dy'])
    derivative=dy[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],t)
    if entry['problem']=='steady_vortex':
        from steady_vortex_exact import physical_raw
        from vortex_rhs_consistency import time_fields as fields
    elif entry['problem']=='nonlinear_simple_wave':
        from simple_wave_exact import physical_raw,fields
    else:raise ValueError('unsupported Euler reference')
    refs=[]
    for order in [12,16]:
        rates,areas,_=exact_rates(cells,edges,a[:,:2],derivative[:,:2],derivative[:,2],t,entry['amplitude'],order,field_function=fields,raw_function=physical_raw)
        refs.append((rates,areas))
    numerical=np.column_stack([np.zeros(len(mass)),derivative[:,3:]])
    metrics={};quadrature={};maxima={}
    for raw,row,target,coarse,area in zip(np.column_stack([mass,a[:,3:]]),numerical,refs[1][0],refs[0][0],refs[1][1]):
        actual=centered_rates(raw,row);exact=centered_rates(raw,target);low=centered_rates(raw,coarse)
        for key in actual:
            # Domain length is1. Density covariance/h^2 and residual velocity
            # covariance/h^2 expose cell scale errors hidden by raw moments.
            scale=entry['amplitude']*(area if key.endswith('covariance') else 1.)
            norm=float(np.linalg.norm(actual[key]-exact[key])/scale)
            gap=float(np.linalg.norm(low[key]-exact[key])/scale)
            metrics[key]=metrics.get(key,0.)+area*norm;quadrature[key]=quadrature.get(key,0.)+area*gap;maxima[key]=max(maxima.get(key,0.),norm)
    result=dict(label=entry['label'],problem=entry['problem'],n=initial.get('n',entry['n']),t=t,pressure_degree=entry['degree'],
        centered_rate_L1_over_scale=metrics,centered_rate_maximum_over_scale=maxima,reference_quadrature_difference=quadrature,
        exact_reference_mass_rate_L1_over_amplitude=float(np.sum(abs(refs[1][0][:,0]))/entry['amplitude']),
        initial_sha256=hashlib.sha256(Path(entry['initial']).read_bytes()).hexdigest(),rhs_sha256=hashlib.sha256(Path(entry['rhs']).read_bytes()).hexdigest(),
        note='Initial smooth Euler consistency only. Derivative of moving centroid/central moments includes exact nonzero reference dm. Covariances normalized by cell area=h^2; other rates by benchmark amplitude. 12/16 reference quadratures are numerical cross-checks, not analytic enclosures. No state or derivative used to advance is modified.')
    atomic(entry['label']+'_moment_rate_audit.json',result);print(entry['label'],metrics,flush=True);return result


def run():
    install_domain();cfg=json.loads(Path('physical_moment_rate_audit_config.json').read_text());rows=[inspect(v) for v in cfg['cases']]
    groups={}
    for row in rows:groups.setdefault((row['problem'],row['pressure_degree']),[]).append(row)
    orders=[]
    for (problem,degree),items in groups.items():
        for a,b in zip(sorted(items,key=lambda v:v['n'])[:-1],sorted(items,key=lambda v:v['n'])[1:]):
            order={}
            for key,small in a['centered_rate_L1_over_scale'].items():
                large=b['centered_rate_L1_over_scale'][key];screen=20*max(a['reference_quadrature_difference'][key],b['reference_quadrature_difference'][key],1e-12)
                order[key]=None if min(small,large)<=screen else float(np.log(small/large)/np.log(b['n']/a['n']))
            orders.append(dict(problem=problem,pressure_degree=degree,n=[a['n'],b['n']],observed_orders=order))
    atomic('physical_moment_rate_audit_report.json',dict(completed=True,rows=rows,spatial_rate_orders=orders))


if __name__=='__main__':run()
