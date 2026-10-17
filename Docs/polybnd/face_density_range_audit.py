"""Near-optimal LP trace ranges and frozen-speed one-sided ALE transport.
These are independent per-face extrema, not a joint ALE/HLL solution.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy import sparse
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from subcell_density_audit import triangles,triangle_moments,normalized
from area_prior_density_audit import interface_averages

OPTIONS={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9}


def trace_weights(parts,a,b):
    """Interior trace of constant triangle densities on a convex cell face."""
    direction=b-a;length=np.linalg.norm(direction);tangent=direction/length
    weights=np.zeros(len(parts));midpoints=np.zeros((len(parts),2))
    for k,v in enumerate(parts):
        for p,q in zip(v,np.roll(v,-1,axis=0)):
            dp=p-a;dq=q-a
            crossp=abs(dp[0]*tangent[1]-dp[1]*tangent[0])
            crossq=abs(dq[0]*tangent[1]-dq[1]*tangent[0])
            if max(crossp,crossq)>1e-10*max(1.,length):continue
            lo=max(0.,min(dp@tangent,dq@tangent));hi=min(length,max(dp@tangent,dq@tangent))
            if hi-lo>1e-12*length:
                weights[k]+=(hi-lo)/length
                midpoints[k]=a+.5*(lo+hi)*tangent
    if abs(weights.sum()-1)>1e-8:raise RuntimeError('face trace does not partition face')
    return weights,midpoints


def problem(vertices,target,level,metadata,means,low,high,parts_override=None):
    parts=triangles(vertices,level) if parts_override is None else parts_override
    raw=np.array([triangle_moments(v) for v in parts]);mass=target[0]
    geo=raw.sum(0);center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);k=len(parts)
    A=np.column_stack([normalized(r,center,scale)/r[0] for r in raw]);b=normalized(target,center,scale)/mass
    if metadata['flat']:
        lv=hv=means;fraction=np.zeros(k)
    else:
        values=interface_averages(parts,None,metadata,low,high)
        lv=low;hv=high;fraction=np.clip((values-low)/(high-low),0,1)
    l=lv*raw[:,0]/mass;h=hv*raw[:,0]/mass;prior=(1-fraction)*l+fraction*h
    I=sparse.eye(k,format='csr')
    Aub=sparse.vstack([sparse.hstack([-I,-I]),sparse.hstack([sparse.diags(1-2*fraction),-I]),sparse.hstack([I,-I])],format='csr')
    bub=np.r_[-prior,(1-fraction)*l-fraction*h,prior]
    Aeq=sparse.hstack([sparse.csr_matrix(A),sparse.csr_matrix((6,k))],format='csr')
    lower=low*.99*raw[:,0]/mass;upper=high*1.01*raw[:,0]/mass
    bounds=list(zip(lower,upper))+[(0,None)]*k
    c=np.r_[np.zeros(k),np.ones(k)]
    base=linprog(c,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=b,bounds=bounds,method='highs',options=OPTIONS)
    if not base.success:raise RuntimeError('base LP failed: '+base.message)
    direct=lambda z:float(np.sum((1-fraction)*abs(z-l)+fraction*abs(z-h)))
    return dict(parts=parts,raw=raw,mass=mass,k=k,A=A,b=b,Aub=Aub,bub=bub,Aeq=Aeq,bounds=bounds,
        lower=lower,upper=upper,c=c,base=base,spatial_objective=direct)


def extremum(p,coefficient,budget,sign):
    k=p['k'];c=np.r_[sign*coefficient,np.zeros(k)]
    Aub=sparse.vstack([p['Aub'],sparse.csr_matrix(p['c'][None,:])],format='csr')
    cap=float(p['base'].fun)+budget
    sol=linprog(c,A_ub=Aub,b_ub=np.r_[p['bub'],cap],A_eq=p['Aeq'],b_eq=p['b'],bounds=p['bounds'],method='highs',options=OPTIONS)
    if not sol.success:return dict(success=False,reason=sol.message)
    z=sol.x[:k];objective=p['spatial_objective'](z)
    moment=float(np.max(abs(p['A']@z-p['b'])))
    bound=float(max(0.,np.max(p['lower']-z),np.max(z-p['upper'])))
    excess=max(0.,objective-cap)
    return dict(success=True,value=float(coefficient@z),moment_error=moment,mass_relative_error=float(abs(z.sum()-1)),
        fraction_bound_error=bound,objective_excess=excess,objective=objective,
        gate=moment<1e-8 and bound<1e-8 and excess<1e-8)


def audit(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    rebuilt,_=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    state=np.asarray(cp['state']);a=state[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
    previous=json.loads(Path('predecessor_completion.json').read_text())['cases'][name]
    # Freeze the face speed from the original smooth closure at this saved state.
    derivative,diag=model.rhs(cp['t'],state,mass);rates=derivative[:-4].reshape(-1,20);gen=rates[:,:2];wdot=rates[:,2]
    states=[model.decode(v,np.r_[mass[i],a[i,3:]]) for i,v in enumerate(cells)]
    faces=[[] for _ in mass]
    for edge_id,(i,j,va,vb,n,length) in enumerate(edges):
        faces[i].append((edge_id,j,va,vb,n,length,i,j))
        if j>=0:faces[j].append((edge_id,i,va,vb,-n,length,i,j))
    rows=[]
    for level in (2,3):
        pred=next(r for r in previous['comparisons'] if r['level']==level and r['angle_degrees']==0.)
        if not pred['all_gates']:raise RuntimeError('predecessor gate failure')
        for i,v in enumerate(cells):
            metadata=pred['cells'][i]['prior'];p=problem(v,np.r_[mass[i],a[i,3:9]],level,metadata,means[i],low,high)
            difference=abs(p['spatial_objective'](p['base'].x[:p['k']])-pred['cells'][i]['inferred_prior_L1_relative'])
            if difference>1e-7:raise RuntimeError('base optimum inconsistent with predecessor')
            for edge_id,neighbor,va,vb,normal,length,ei,ej in faces[i]:
                weights,midpoints=trace_weights(p['parts'],va,vb)
                density_coeff=weights*p['mass']/p['raw'][:,0]/high
                relative=np.zeros(p['k'])
                for q in np.flatnonzero(weights):
                    point=midpoints[q]
                    if ej>=0:
                        distance=np.linalg.norm(a[ej,:2]-a[ei,:2])
                        wn=((a[ej,:2]-a[ei,:2])@gen[ej]-(gen[ej]-gen[ei])@(point-a[ei,:2]))/distance+(wdot[ei]-wdot[ej])/(2*distance)
                        if i==ej:wn=-wn
                    else:wn=float(model.boundary_velocity(point)@normal)
                    relative[q]=float(model.velocity(states[i],point)@normal-wn)
                speed=float(np.max(abs(relative)))
                # Frozen one-sided advective mass transport, not re-solved HLL flux.
                flux_coeff=weights*p['mass']/p['raw'][:,0]*relative/(high*max(speed,1e-12))
                for budget in (1e-8,1e-5):
                    density=[extremum(p,density_coeff,budget,sign) for sign in (1.,-1.)]
                    flux=[extremum(p,flux_coeff,budget,sign) for sign in (1.,-1.)] if speed>1e-12 else []
                    all_results=density+flux;gate=all(x.get('gate',False) for x in all_results)
                    width=(density[1]['value']-density[0]['value'])*high if gate else None
                    flux_width=(flux[1]['value']-flux[0]['value'])*high*length*speed if flux and gate else 0. if gate else None
                    rows.append(dict(level=level,cell=i,neighbor=neighbor,edge=edge_id,budget=budget,all_gates=gate,
                        density_extrema=density,flux_extrema=flux,density_range_width=width,
                        density_width_over_contrast=width/(high-low) if width is not None else None,
                        frozen_transport_range_width=flux_width,face_length=length,max_relative_speed=speed,
                        base_objective=float(p['base'].fun),base_predecessor_difference=difference))
            atomic(name+'.partial.json',dict(rows=rows))
        print(name,'level',level,'completed',flush=True)
    result=dict(case=name,status='completed',all_gates=all(x['all_gates'] for x in rows),rows=rows,
        frozen_smooth_mass_flux_residual=diag['mass_flux'],inferred_densities=[low,high],
        note='Independent per-face one-sided extrema. Fixed numerical m,q,Q and 1% bounds; objective <= baseline optimum + absolute dimensionless budget. Frozen original smooth face speed/velocity only. No HLL or weight re-solve; no joint attainment or evolution/existence claim.')
    atomic(name+'.json',result);return result


def main():
    if not json.loads(Path('predecessor_completion.json').read_text())['completed']:raise RuntimeError('predecessor incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
        print(name,cases[name]['status'],flush=True)
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
