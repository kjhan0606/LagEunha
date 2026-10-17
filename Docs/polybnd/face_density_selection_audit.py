"""Select face means by a secondary convex objective inside a volume-L1 cap.
Numerical inferred interfaces only enter selection; exact field only scores it.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog,minimize
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from face_density_range_audit import problem,trace_weights

OPTIONS={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10}


def line_prior_mean(a,b,metadata,mean,low,high):
    if metadata['flat']:return mean
    n=np.asarray(metadata['normal']);c=float(n@a-metadata['offset']);d=float(n@(b-a))
    if abs(d)<1e-14:fraction=float(c>=0)
    elif d>0:fraction=1-float(np.clip(-c/d,0,1))
    else:fraction=float(np.clip(-c/d,0,1))
    value=low+(high-low)*fraction;blend=metadata.get('blend',1.)
    return value if blend==1. else (1-blend)*mean+blend*value


def one_sided_interval(vertices,a,b,inward):
    """Clip a line parameter to a convex material polygon, taking interior trace."""
    orientation=np.sign(np.sum(vertices[:,0]*np.roll(vertices[:,1],-1)-vertices[:,1]*np.roll(vertices[:,0],-1)))
    lo,hi=0.,1.
    for p,q in zip(vertices,np.roll(vertices,-1,axis=0)):
        delta=q-p;length=np.linalg.norm(delta)
        if length<1e-14:continue
        n=orientation*np.array([delta[1],-delta[0]])/length
        c=float(n@(a-p));d=float(n@(b-a))
        if abs(d)<1e-13:
            if c>1e-12 or (abs(c)<=1e-12 and n@inward>1e-12):return None
        elif d>0:hi=min(hi,-c/d)
        else:lo=max(lo,-c/d)
        if hi<=lo:return None
    return (lo,hi) if hi>lo else None


def select(p,T,target,weights,budget,tolerance,start):
    """Fully corrective conditional gradient; LP oracle and small simplex QP.

    The reported Frank-Wolfe gap bounds objective suboptimality subject to
    numerical LP accuracy. Positive weights make the selected trace vector
    unique in the exact convex problem, not necessarily the interior field.
    """
    k=p['k'];cap=float(p['base'].fun)+budget
    Aub=sparse.vstack([p['Aub'],sparse.csr_matrix(p['c'][None,:])],format='csr')
    bub=np.r_[p['bub'],cap]
    def oracle(direction):
        scale=float(np.max(abs(direction)))
        if scale<1e-30:return p['base'].x.copy()
        sol=linprog(np.r_[direction/scale,np.zeros(k)],A_ub=Aub,b_ub=bub,
            A_eq=p['Aeq'],b_eq=p['b'],bounds=p['bounds'],method='highs',options=OPTIONS)
        if not sol.success:raise RuntimeError('selection oracle: '+sol.message)
        return sol.x
    x=p['base'].x.copy() if start==0 else oracle(np.sin(np.arange(len(target))+1)@T)
    support=[x];alpha=np.ones(1);fallbacks=0
    for iteration in range(201):
        y=T@x[:k];gradient=weights*(y-target);s=oracle(gradient@T);ys=T@s[:k]
        raw_gap=float(gradient@(y-ys));gap=max(0.,raw_gap)
        if gap<=tolerance:break
        if iteration==200:break
        support.append(s);V=np.array([T@v[:k] for v in support]).T
        def fun(a):
            diff=V@a-target;return .5*float(weights@(diff*diff))
        def jac(a):return V.T@(weights*(V@a-target))
        initial=np.r_[alpha,0.]
        objective_scale=1/max(fun(initial),1e-12)
        qp=minimize(lambda a:objective_scale*fun(a),initial,jac=lambda a:objective_scale*jac(a),method='SLSQP',bounds=[(0.,1.)]*len(support),
            constraints=[dict(type='eq',fun=lambda a:a.sum()-1,jac=lambda a:np.ones(len(a)))],
            options={'ftol':1e-15,'maxiter':200})
        candidate=np.maximum(qp.x,0);candidate/=candidate.sum()
        # Retain the exact Frank-Wolfe line step if the small QP stalls.
        d=ys-y;den=float(weights@(d*d))
        step=float(np.clip(raw_gap/den,0,1)) if den>0 else 0.
        line_candidate=np.r_[(1-step)*alpha,step]
        if not qp.success or fun(line_candidate)<fun(candidate):
            fallbacks+=1;candidate=line_candidate
        keep=candidate>1e-12
        alpha=candidate[keep];alpha/=alpha.sum();support=[v for v,take in zip(support,keep) if take]
        x=alpha@np.asarray(support)
    z=x[:k];y=T@z;objective=p['spatial_objective'](z)
    moment=float(np.max(abs(p['A']@z-p['b'])))
    bound=float(max(0.,np.max(p['lower']-z),np.max(z-p['upper'])))
    excess=max(0.,objective-cap)
    result=dict(iterations=iteration,gap=gap,raw_gap=raw_gap,converged=gap<=tolerance,
        secondary_objective=.5*float(weights@((y-target)**2)),moment_error=moment,
        mass_relative_error=float(abs(z.sum()-1)),fraction_bound_error=bound,
        volume_objective=objective,volume_cap=cap,objective_excess=excess,
        gate=gap<=tolerance and raw_gap>=-1e-8 and moment<1e-8 and bound<1e-8 and excess<1e-8,
        simplex_fallbacks=fallbacks,support=len(support))
    return result,z


def audit(name,prior_builder=None,parts_builder=None,levels=(2,3),budgets=(1e-8,1e-5),output_prefix=None):
    output_prefix=name if output_prefix is None else output_prefix
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    rebuilt,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    a0=y0[:-4].reshape(-1,20);initial,_=model.geometry(a0[:,:2],a0[:,2],0.)
    refs=[];phases=np.where(a0[:,0]<.5,1.,cfg['amplitude'])
    for v in initial:
        ref=v.copy();ref[:,0]+=cp['t']*(.2+1.1*ref[:,1]);refs.append(ref)
    state=np.asarray(cp['state']);a=state[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
    centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells]);neighbors=[set() for _ in mass]
    faces=[[] for _ in mass]
    for edge_id,(i,j,va,vb,n,length) in enumerate(edges):
        faces[i].append((edge_id,j,va,vb,n,length,i,j))
        if j>=0:
            neighbors[i].add(j);neighbors[j].add(i)
            faces[j].append((edge_id,i,va,vb,-n,length,i,j))
    priors=[interface_prior(i,v,centers,means,neighbors)[1] for i,v in enumerate(cells)]
    if prior_builder is not None:
        priors=[prior_builder(v,np.r_[mass[i],a[i,3:9]],means,priors[i]) for i,v in enumerate(cells)]
    derivative,_=model.rhs(cp['t'],state,mass);rates=derivative[:-4].reshape(-1,20);gen=rates[:,:2];wdot=rates[:,2]
    states=[model.decode(v,np.r_[mass[i],a[i,3:]]) for i,v in enumerate(cells)]
    def relative(point,normal,ei,ej,i):
        if ej>=0:
            distance=np.linalg.norm(a[ej,:2]-a[ei,:2])
            wn=((a[ej,:2]-a[ei,:2])@gen[ej]-(gen[ej]-gen[ei])@(point-a[ei,:2]))/distance+(wdot[ei]-wdot[ej])/(2*distance)
            if i==ej:wn=-wn
        else:wn=float(model.boundary_velocity(point)@normal)
        return float(model.velocity(states[i],point)@normal-wn)
    resume_path=Path('resume_records.json')
    saved=json.loads(resume_path.read_text()).get(name,[]) if resume_path.exists() and prior_builder is None and parts_builder is None else []
    reusable={(r['level'],r['cell'],r['budget']):r for r in saved if r['all_gates']}
    records=[]
    for level in levels:
        for i,v in enumerate(cells):
            if all((level,i,budget) in reusable for budget in budgets):
                records.extend(dict(reusable[(level,i,budget)],reused=True) for budget in budgets)
                continue
            parts=None if parts_builder is None else parts_builder(v,level,priors[i])
            p=problem(v,np.r_[mass[i],a[i,3:9]],level,priors[i],means[i],low,high,parts_override=parts)
            T=[];transport=[];target=[];exact=[];lengths=[];reference_flux=[];coverage=[]
            for edge_id,j,va,vb,normal,length,ei,ej in faces[i]:
                weights,mid=trace_weights(p['parts'],va,vb)
                T.append(weights*p['mass']/p['raw'][:,0]/high)
                rel=np.array([relative(point,normal,ei,ej,i) if weight else 0. for point,weight in zip(mid,weights)])
                transport.append(weights*p['mass']/p['raw'][:,0]*rel*length)
                target.append(line_prior_mean(va,vb,priors[i],means[i],low,high)/high)
                face_density=face_flux=covered=0.
                for ref,rho in zip(refs,phases):
                    interval=one_sided_interval(ref,va,vb,-normal)
                    if interval is None:continue
                    lo,hi=interval;fraction=hi-lo;point=va+.5*(lo+hi)*(vb-va)
                    face_density+=rho*fraction;face_flux+=rho*length*fraction*relative(point,normal,ei,ej,i);covered+=fraction
                exact.append(face_density);reference_flux.append(face_flux);coverage.append(abs(covered-1));lengths.append(length)
            T=np.asarray(T);transport=np.asarray(transport);target=np.asarray(target);lengths=np.asarray(lengths)
            if max(coverage)>1e-8:raise RuntimeError('exact trace coverage failed')
            weights=lengths/lengths.sum();base=p['base'].x[:p['k']]
            # Exact-field overlaps only score reconstructions; never enter select().
            intersections=[]
            for triangle in p['parts']:
                poly=Polygon(triangle);intersections.append([poly.intersection(Polygon(ref)).area for ref in refs])
            intersections=np.asarray(intersections)
            def score(z):
                rho=z*p['mass']/p['raw'][:,0]
                error=float(np.sum(intersections*abs(rho[:,None]-phases[None,:])))
                bound_error=max(0.,low*.99-float(rho.min()),float(rho.max())-high*1.01)/max(1.,high*1.01)
                return dict(face_density=(T@z*high).tolist(),transport=(transport@z).tolist(),density_L1_absolute=error,
                    min_density=float(rho.min()),max_density=float(rho.max()),density_bound_error_scaled=bound_error,
                    density_bound_gate=bound_error<1e-8,mass_fractions=z.tolist())
            baseline=score(base)
            for budget in budgets:
                if (level,i,budget) in reusable:
                    records.append(dict(reusable[(level,i,budget)],reused=True));continue
                results=[]
                for start,tolerance in ((0,1e-10),(1,1e-12)):
                    result,z=select(p,T,target,weights,budget,tolerance,start)
                    result.update(start=start,tolerance=tolerance,**score(z))
                    result['gate']=result['gate'] and result['density_bound_gate'];results.append(result)
                records.append(dict(level=level,cell=i,budget=budget,reused=False,selected=results,baseline=baseline,
                    subcells=p['k'],
                    inferred_target=(target*high).tolist(),exact_face_density=exact,exact_transport=reference_flux,
                    lengths=lengths.tolist(),face_ids=[f[0] for f in faces[i]],max_reference_coverage_error=max(coverage),
                    all_gates=baseline['density_bound_gate'] and all(x['gate'] for x in results),
                    start_face_difference=float(np.max(abs(np.array(results[0]['face_density'])-results[1]['face_density']))),
                    start_transport_difference=float(np.max(abs(np.array(results[0]['transport'])-results[1]['transport'])))))
            atomic(output_prefix+'.partial.json',dict(records=records))
        print(name,'level',level,'completed',flush=True)
    result=dict(case=name,status='completed',all_gates=all(x['all_gates'] for x in records),records=records,
        inferred_priors=priors,
        total_mass=float(mass.sum()),note='Positive face-length-weighted square error to inferred face means, under volume-L1 cap and fixed m,q,Q/1% bounds. Two starts/tolerances. Frozen original face speed; one-sided transport only. Exact field scores only. Unique face-mean vector in exact convex problem does not imply unique subface density or flux.')
    atomic(output_prefix+'.json',result);return result


def main():
    previous=json.loads(Path('predecessor_completion.json').read_text())
    if not previous['completed'] or not all(c.get('all_gates') for c in previous['cases'].values()):raise RuntimeError('predecessor gate failure')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
        print(name,cases[name]['status'],flush=True)
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
