"""Geometry-only boundary layers: oracle face-moment feasibility diagnostic."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from scipy.optimize import linprog
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from subcell_density_audit import triangles,triangle_moments
from moment_aligned_density_audit import polygon_moments
from face_moment_feasibility_audit import build,solve,capped,face_matrix

RADIAL=(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
SEGMENTS=(4,16,64)


def boundary_parts(vertices,radial,segments,core_level=2,moment_function=None,core_prior=None):
    if not 0<radial<1 or segments<1:raise ValueError('invalid boundary layer parameters')
    center=vertices.mean(0);inner=center+(1-radial)*(vertices-center)
    if core_prior is None:parts=triangles(inner,core_level)
    else:
        from moment_aligned_density_audit import aligned_parts
        parts=aligned_parts(inner,core_level,core_prior)
    for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
        for j in range(segments):
            p=a+(b-a)*(j/segments);q=a+(b-a)*((j+1)/segments)
            pi=center+(1-radial)*(p-center);qi=center+(1-radial)*(q-center)
            parts.extend([np.array([p,q,qi]),np.array([p,qi,pi])])
    raw=np.array([triangle_moments(v) for v in parts]) if moment_function is None else moment_function(parts)
    original=polygon_moments(vertices)
    if raw[:,0].min()<=0:raise RuntimeError('nonpositive boundary triangle')
    np.testing.assert_allclose(raw.sum(0),original,rtol=1e-10,atol=1e-14)
    return parts


def trace_interval(mean,lower,upper,segments):
    """Exact first-moment interval for equal constant bins at a fixed mean.

    Rearrangement maximizes moment by filling the rightmost bins first.
    All moments here are normalized by H.
    """
    if mean<lower-1e-10 or mean>upper+1e-10:return dict(mean_feasible=False)
    remaining=float(np.clip(mean-lower,0,upper-lower));maximum=0.
    for j in range(segments-1,-1,-1):
        amount=min(remaining,(upper-lower)/segments)
        maximum+=amount*(-1+(2*j+1)/segments);remaining-=amount
    return dict(mean_feasible=True,minimum_first=-maximum,maximum_first=maximum)


def boundary_face_matrix(vertices,parts_count,segments,a,b,core_level=2,core_count=None):
    # Boundary triangle ordering is explicit: only the first triangle of
    # each trapezoid touches its outer face over a nonzero interval.
    if core_count is None:core_count=len(vertices)*4**core_level
    if parts_count!=core_count+len(vertices)*2*segments:raise RuntimeError('unexpected part ordering')
    for edge,(p,q) in enumerate(zip(vertices,np.roll(vertices,-1,axis=0))):
        forward=np.linalg.norm(p-a)+np.linalg.norm(q-b)
        reverse=np.linalg.norm(p-b)+np.linalg.norm(q-a)
        if min(forward,reverse)>1e-9:continue
        sign=1. if forward<reverse else -1.
        ids=core_count+edge*2*segments+2*np.arange(segments)
        w0=np.zeros(parts_count);w1=np.zeros(parts_count)
        w0[ids]=1/segments;w1[ids]=sign*(-1+(2*np.arange(segments)+1)/segments)/segments
        return w0,w1
    raise RuntimeError('cell face not found')


def forget_witness(result):
    result=dict(result);result.pop('mass_fractions',None);return result


def audit(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    rebuilt,_=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
    centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells]);neighbors=[set() for _ in mass];faces=[[] for _ in mass]
    for edge_id,(i,j,va,vb,n,length) in enumerate(edges):
        faces[i].append((edge_id,va,vb))
        if j>=0:faces[j].append((edge_id,va,vb));neighbors[i].add(j);neighbors[j].add(i)
    priors=[interface_prior(i,v,centers,means,neighbors)[1] for i,v in enumerate(cells)]
    previous=json.loads(Path('predecessor_completion.json').read_text())['cases'][name]
    saved={r['cell']:r for r in previous['rows'] if r['variant']=='neighbor_uniform' and r['level']==2}
    targets=[]
    for i in range(len(mass)):
        if saved[i]['face_ids']!=[f[0] for f in faces[i]]:raise RuntimeError('saved face ordering mismatch')
        targets.append((np.array(saved[i]['exact_face_means'])/high,np.array(saved[i]['exact_first_moments'])/high))
    trace_checks={}
    for segments in SEGMENTS:
        checks=[]
        for i,(m0,m1) in enumerate(targets):
            intervals=[trace_interval(mean,low*.99/high,1.01,segments) for mean in m0]
            excess=max(max(0.,first-r['maximum_first'],r['minimum_first']-first) if r['mean_feasible'] else np.inf for first,r in zip(m1,intervals))
            checks.append(dict(cell=i,feasible=excess<=1e-10,max_first_excess=excess,intervals=intervals))
        trace_checks[str(segments)]=checks
    rows=[];logical_reuses=0
    resume_path=Path('resume_rows.json')
    saved=json.loads(resume_path.read_text()).get(name,[]) if resume_path.exists() else []
    reusable={(r['radial'],r['segments'],r['cell']):r for r in saved if r['numerical_valid'] and r['mean_and_first']['feasible']}
    for radial in RADIAL:
        for segments in SEGMENTS:
            config=[]
            for i,v in enumerate(cells):
                if (radial,segments,i) in reusable:
                    row=dict(reusable[(radial,segments,i)],reused=True)
                    config.append(row);rows.append(row);continue
                parts=boundary_parts(v,radial,segments);p=build(v,np.r_[mass[i],a[i,3:9]],parts,priors[i],means[i],low,high)
                p['moment_row_scale']=1e4
                baseline=solve(p);m0,m1=targets[i];T0=[];T1=[]
                for edge_id,va,vb in faces[i]:
                    w0,w1=boundary_face_matrix(v,len(parts),segments,va,vb);T0.append(w0);T1.append(w1)
                T0=np.asarray(T0);T1=np.asarray(T1)
                if baseline['classification']=='solver_infeasible':
                    mean=joint=dict(classification='solver_infeasible',feasible=False,implied_by='volume_moments_infeasible');logical_reuses+=2
                elif not baseline['feasible']:
                    raise RuntimeError('baseline numerical failure: '+str(forget_witness(baseline)))
                else:
                    mean=solve(p,T0,m0)
                    if mean['classification']=='solver_infeasible':
                        joint=dict(classification='solver_infeasible',feasible=False,implied_by='face_means_infeasible');logical_reuses+=1
                    elif not trace_checks[str(segments)][i]['feasible']:
                        joint=dict(classification='solver_infeasible',feasible=False,implied_by='tangential_trace_interval');logical_reuses+=1
                    else:joint=solve(p,np.vstack([T0,T1]),np.r_[m0,m1])
                if baseline['feasible']:
                    mc=capped(mean,baseline);jc=capped(joint,baseline)
                else:mc=jc=dict(classification='infeasible_from_volume_moments')
                bad=any(s['classification']=='numerical_failure' for s in (baseline,mean,joint)) or any(s['classification']=='indeterminate' for s in (mc,jc))
                row=dict(radial=radial,segments=segments,core_level=2,cell=i,subcells=p['k'],reused=False,moment_row_scale=1e4,
                    min_triangle_area=float(p['raw'][:,0].min()),baseline=forget_witness(baseline),
                    mean_only=forget_witness(mean),mean_and_first=joint,mean_capped=mc,joint_capped=jc,
                    numerical_valid=not bad)
                config.append(row);rows.append(row)
            atomic(name+'_r'+format(radial,'.0e')+'_n'+str(segments)+'.json',dict(rows=config))
            print(name,'radial',radial,'segments',segments,'joint',sum(x['mean_and_first']['feasible'] for x in config),'of',len(mass),flush=True)
    result=dict(case=name,status='completed' if all(x['numerical_valid'] for x in rows) else 'numerical_failure',
        rows=rows,trace_checks=trace_checks,logical_reuses=logical_reuses,
        note='Homothetic boundary layers independent of inferred/exact interfaces. Corelevel2; normal thickness radial*center-to-face distance, equal tangential bins. Neighbor prior only sets volume objective. Saved oracle M0/M1 targets diagnostic only. Fixed saved m,q,Q and1% bounds; cap delta1e-5 derived from minima. No evolution/HLL solve.')
    atomic(name+'.json',result);return result


def main():
    previous=json.loads(Path('predecessor_completion.json').read_text())
    if not previous['completed']:raise RuntimeError('predecessor incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:
            cases[name]=dict(status='failed',error=repr(exc));atomic(name+'.error.json',cases[name])
        print(name,cases[name]['status'],flush=True)
    report=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',report)
    raise SystemExit(0 if report['completed'] else 2)


if __name__=='__main__':main()
