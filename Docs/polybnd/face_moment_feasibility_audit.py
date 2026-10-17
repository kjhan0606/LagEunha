"""Oracle face-moment feasibility on saved states; diagnostic use only.
All exact face constraints on a cell are imposed simultaneously. No evolution.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from subcell_density_audit import triangles,triangle_moments,normalized
from area_prior_density_audit import interface_averages
from moment_aligned_density_audit import aligned_parts
from face_density_range_audit import trace_weights
from face_density_selection_audit import one_sided_interval

OPTIONS={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10}


def reference_face_moments(refs,phases,a,b,normal):
    """M0=int_0^1 rho dt, M1=int_0^1 (2t-1)rho dt, interior trace."""
    m0=m1=coverage=0.
    for ref,rho in zip(refs,phases):
        interval=one_sided_interval(ref,a,b,-normal)
        if interval is None:continue
        lo,hi=interval;fraction=hi-lo
        m0+=rho*fraction;m1+=rho*fraction*(lo+hi-1);coverage+=fraction
    if abs(coverage-1)>1e-8:raise RuntimeError('exact reference coverage failed')
    return m0,m1,abs(coverage-1)


def face_matrix(parts,a,b):
    weights,mid=trace_weights(parts,a,b)
    tangent=b-a;s=(mid-a)@tangent/(tangent@tangent)
    return weights,weights*(2*s-1)


def build(vertices,target,parts,metadata,mean,low,high,geometry_moments=None,area_prior=None):
    raw=np.array([triangle_moments(v) for v in parts]) if geometry_moments is None else geometry_moments
    geo=raw.sum(0);mass=target[0];k=len(parts)
    center=geo[1:3]/geo[0];scale=np.sqrt(geo[0])
    # Variables u=rho/high and v=pointwise L1 epigraph/high. Every u has
    # ordinary-sized positive bounds, even for arbitrarily tiny subtriangles.
    if geometry_moments is None:A=high*np.column_stack([normalized(r,center,scale) for r in raw])/mass
    else:
        from vector_boundary_geometry import normalized_many
        A=high*normalized_many(raw,center,scale).T/mass
    b=normalized(target,center,scale)/mass
    if metadata['flat']:
        lv=hv=mean/high;fraction=np.zeros(k)
    else:
        means=interface_averages(parts,None,metadata,low,high) if area_prior is None else area_prior
        blend=metadata.get('blend',1.)
        lv=low/high if blend==1. else ((1-blend)*mean+blend*low)/high
        hv=1. if blend==1. else ((1-blend)*mean+blend*high)/high
        fraction=np.clip((means-low)/(high-low),0,1)
    prior=(1-fraction)*lv+fraction*hv;I=sparse.eye(k,format='csr')
    Aub=sparse.vstack([sparse.hstack([-I,-I]),sparse.hstack([sparse.diags(1-2*fraction),-I]),sparse.hstack([I,-I])],format='csr')
    bub=np.r_[-prior,(1-fraction)*lv-fraction*hv,prior]
    Aeq=sparse.hstack([sparse.csr_matrix(A),sparse.csr_matrix((6,k))],format='csr')
    area_scale=high*raw[:,0]/mass
    c=np.r_[np.zeros(k),area_scale]
    bounds=[(low*.99/high,1.01)]*k+[(0,None)]*k
    direct=lambda u:float(np.sum(area_scale*((1-fraction)*abs(u-lv)+fraction*abs(u-hv))))
    return dict(k=k,A=A,b=b,Aeq=Aeq,Aub=Aub,bub=bub,c=c,bounds=bounds,high=high,low=low,
        spatial_objective=direct,raw=raw,mass=mass,parts=parts,prior_density_units=prior)


def solve(p,F=None,face_target=None):
    k=p['k'];Aeq=p['Aeq'];beq=p['b']
    if F is not None:
        Aeq=sparse.vstack([Aeq,sparse.hstack([sparse.csr_matrix(F),sparse.csr_matrix(F.shape)])],format='csr')
        beq=np.r_[beq,face_target]
    row_scale=np.ones(len(beq));row_scale[:6]=p.get('moment_row_scale',1.)
    scaled_eq=sparse.diags(row_scale)@Aeq
    def optimize(method,options):
        return linprog(p['c'],A_ub=p['Aub'],b_ub=p['bub'],A_eq=scaled_eq,b_eq=beq*row_scale,
            bounds=p['bounds'],method=method,options={**OPTIONS,**options})
    sol=optimize('highs',{});attempts=[]
    if not sol.success and sol.status!=2 and p.get('retry_numerical_failure'):
        attempts.append(dict(method='highs',presolve=True,status=int(sol.status),message=sol.message))
        for method in ('highs-ds','highs-ipm'):
            sol=optimize(method,dict(presolve=False))
            attempts.append(dict(method=method,presolve=False,status=int(sol.status),message=sol.message))
            if sol.success:break
    if attempts and not sol.success:return dict(classification='numerical_failure',feasible=False,solver_status=int(sol.status),reason=sol.message,solver_attempts=attempts)
    if sol.status==2:return dict(classification='solver_infeasible',feasible=False,solver_status=2,reason=sol.message)
    if not sol.success:return dict(classification='numerical_failure',feasible=False,solver_status=int(sol.status),reason=sol.message)
    u=sol.x[:k];rho=u*p['high'];objective=p['spatial_objective'](u)
    moment=float(np.max(abs(p['A']@u-p['b'])))
    face=float(np.max(abs(F@u-face_target))) if F is not None else 0.
    bound=max(0.,p['low']*.99-float(rho.min()),float(rho.max())-p['high']*1.01)/max(1.,p['high']*1.01)
    epigraph_error=abs(objective-float(sol.fun))
    # HiGHS dual objective and stationarity residual, for numerical cap decisions.
    eqdual=sol.eqlin.marginals*row_scale;ubdual=sol.ineqlin.marginals
    lower=np.asarray(sol.lower.marginals);upper=np.asarray(sol.upper.marginals)
    lo=np.array([a for a,b in p['bounds']]);hi=np.array([np.inf if b is None else b for a,b in p['bounds']])
    dual=float(beq@eqdual+p['bub']@ubdual+lo@lower+hi[np.isfinite(hi)]@upper[np.isfinite(hi)])
    stationarity=float(np.max(abs(p['c']-Aeq.T@eqdual-p['Aub'].T@ubdual-lower-upper)))
    gap=abs(float(sol.fun)-dual)
    sign_error=float(max(0.,np.max(ubdual),-np.min(lower),np.max(upper)))
    gate=moment<1e-8 and face<1e-8 and bound<1e-8 and epigraph_error<1e-8 and stationarity<1e-8 and gap<1e-8 and sign_error<1e-8
    result=dict(classification='feasible_validated' if gate else 'numerical_failure',feasible=gate,solver_status=0,
        objective=objective,linear_objective=float(sol.fun),dual_objective=dual,duality_gap=gap,dual_stationarity_error=stationarity,dual_sign_error=sign_error,
        moment_error=moment,face_error=face,density_bound_error=bound,min_density=float(rho.min()),max_density=float(rho.max()),
        mass_relative_error=float(abs(p['A'][0]@u-1)),epigraph_error=epigraph_error,
        mass_fractions=(rho*p['raw'][:,0]/p['mass']).tolist())
    if p.get('keep_dual_certificate'):
        result['certificate']=dict(equality=eqdual,inequality=ubdual)
    if attempts:result['solver_attempts']=attempts
    return result


def capped(result,baseline,delta=1e-5):
    if result['classification']=='solver_infeasible':return dict(classification='infeasible_from_uncapped')
    if result['classification']!='feasible_validated':return dict(classification='indeterminate')
    cap=baseline['objective']+delta
    if result['objective']<=cap+1e-9:
        return dict(classification='feasible_with_witness',cap=cap,excess=max(0.,result['objective']-cap))
    if result['dual_objective']>cap+1e-8:
        return dict(classification='minimum_objective_exceeds_cap',cap=cap,minimum_excess=result['dual_objective']-cap)
    return dict(classification='indeterminate',cap=cap,objective=result['objective'])


def audit(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    rebuilt,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    a0=y0[:-4].reshape(-1,20);initial,_=model.geometry(a0[:,:2],a0[:,2],0.);refs=[];phases=np.where(a0[:,0]<.5,1.,cfg['amplitude'])
    for v in initial:
        ref=v.copy();ref[:,0]+=cp['t']*(.2+1.1*ref[:,1]);refs.append(ref)
    a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
    centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells]);neighbors=[set() for _ in mass];faces=[[] for _ in mass]
    for edge_id,(i,j,va,vb,n,length) in enumerate(edges):
        faces[i].append((edge_id,va,vb,n))
        if j>=0:
            faces[j].append((edge_id,va,vb,-n));neighbors[i].add(j);neighbors[j].add(i)
    previous=json.loads(Path('predecessor_completion.json').read_text())['cases'][name]['variants']
    neighbor_priors=[interface_prior(i,v,centers,means,neighbors)[1] for i,v in enumerate(cells)]
    moment_priors=previous['moment_uniform']['inferred_priors']
    reference=[]
    for i,v in enumerate(cells):
        reference.append([reference_face_moments(refs,phases,va,vb,normal) for edge_id,va,vb,normal in faces[i]])
    rows=[];logical_reuses=0
    for tag in ('neighbor_uniform','moment_uniform','neighbor_aligned','moment_aligned'):
        priors=moment_priors if tag.startswith('moment') else neighbor_priors
        for level in (2,3):
            predecessor={x['cell']:x for x in previous[tag]['records'] if x['level']==level}
            for i,v in enumerate(cells):
                parts=aligned_parts(v,level,priors[i]) if tag.endswith('aligned') else triangles(v,level)
                p=build(v,np.r_[mass[i],a[i,3:9]],parts,priors[i],means[i],low,high)
                baseline=solve(p)
                if not baseline['feasible']:raise RuntimeError('unconstrained baseline invalid: '+str(baseline))
                original=predecessor[i]['selected'][1]['volume_cap']-1e-5
                if abs(baseline['objective']-original)>1e-7:raise RuntimeError('baseline disagrees with prior batch')
                T0=[];T1=[]
                for edge_id,va,vb,normal in faces[i]:
                    w0,w1=face_matrix(parts,va,vb);T0.append(w0);T1.append(w1)
                T0=np.asarray(T0);T1=np.asarray(T1);target0=np.array([r[0]/high for r in reference[i]]);target1=np.array([r[1]/high for r in reference[i]])
                mean=solve(p,T0,target0)
                if mean['classification']=='solver_infeasible':
                    both=dict(classification='solver_infeasible',feasible=False,implied_by='mean_constraints_infeasible');logical_reuses+=1
                else:both=solve(p,np.vstack([T0,T1]),np.r_[target0,target1])
                rows.append(dict(variant=tag,level=level,cell=i,subcells=p['k'],baseline=baseline,
                    mean_only=mean,mean_and_first=both,mean_capped=capped(mean,baseline),joint_capped=capped(both,baseline),
                    max_reference_coverage_error=max(x[2] for x in reference[i]),face_ids=[x[0] for x in faces[i]],
                    exact_face_means=(target0*high).tolist(),exact_first_moments=(target1*high).tolist()))
            atomic(name+'.partial.json',dict(rows=rows))
            print(name,tag,'level',level,'completed',flush=True)
    failures=[(x['variant'],x['level'],x['cell']) for x in rows if x['mean_only']['classification']=='numerical_failure' or x['mean_and_first']['classification']=='numerical_failure'
        or x['mean_capped']['classification']=='indeterminate' or x['joint_capped']['classification']=='indeterminate']
    result=dict(case=name,status='completed' if not failures else 'numerical_failure',rows=rows,failures=failures,logical_reuses=logical_reuses,
        note='Exact face constraints are oracle diagnostics only. All faces imposed simultaneously. Fixed saved m,q,Q and 1% bounds. Finite-basis HiGHS feasibility; no continuum nonexistence proof. Cap feasibility derived from minimum prior-volume L1, without redundant capped solves. Physical density variables used for scaling.')
    atomic(name+'.json',result);return result


def main():
    previous=json.loads(Path('predecessor_completion.json').read_text())
    if not previous['completed'] or not previous['all_gates']:raise RuntimeError('predecessor failed')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
        print(name,cases[name]['status'],flush=True)
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
