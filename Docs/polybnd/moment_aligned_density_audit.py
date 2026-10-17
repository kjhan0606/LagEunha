"""Factorial comparison of moment-based normals and interface-aligned parts."""
import json,os
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from shapely.geometry import Polygon
from discontinuous_density_audit import clip
from subcell_density_audit import triangles,triangle_moments,normalized
from face_density_selection_audit import audit as selection_audit
from adaptive_coupled_mass import atomic


def polygon_moments(vertices):
    if len(vertices)<3:return np.zeros(7)
    return sum((triangle_moments(np.array([vertices[0],vertices[j],vertices[j+1]]))
        for j in range(1,len(vertices)-1)),np.zeros(7))


def cut_at_fraction(vertices,normal,fraction):
    if os.environ.get('LAGFORCE_ANALYTIC_CAP','0')=='1':
        from analytic_cap_inverse import cap_at_fraction
        return cap_at_fraction(vertices,normal,fraction)
    lo,hi=float((vertices@normal).min()),float((vertices@normal).max());area=Polygon(vertices).area
    for _ in range(48):
        mid=(lo+hi)/2
        if clip(vertices,normal,mid).area>fraction*area:lo=mid
        else:hi=mid
    offset=(lo+hi)/2;cut=clip(vertices,normal,offset)
    return offset,cut


def _fit_moment_prior(vertices,target,means,neighbor,refine=False):
    if os.environ.get('LAGFORCE_CENTERED_MOMENT_FIT','0')=='1':
        from centered_moment_prior import fit_centered_prior
        return fit_centered_prior(vertices,target,means,neighbor,refine=refine)
    low=float(means.min());high=float(means.max());geo=polygon_moments(vertices);mass=target[0]
    fraction=float(np.clip((mass/geo[0]-low)/(high-low),0,1)) if high>low else 0.
    if high-low<1e-12 or fraction<=1e-14 or fraction>=1-1e-14:
        return dict(flat=True,moment_fit_skipped=True)
    center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);goal=normalized(target,center,scale)/mass
    weight=np.array([1.,1.,1.,1.,np.sqrt(2.),1.])
    def residual(theta):
        normal=np.array([np.cos(theta),np.sin(theta)]);offset,cut=cut_at_fraction(vertices,normal,fraction)
        raw=low*geo
        if not cut.is_empty:raw=raw+(high-low)*polygon_moments(np.asarray(cut.exterior.coords[:-1]))
        diff=(normalized(raw,center,scale)/mass-goal)*weight
        return diff,normal,offset
    def candidate(theta):
        diff,normal,offset=residual(theta)
        return float(diff@diff),normal,offset
    baseline=float(np.arctan2(neighbor['normal'][1],neighbor['normal'][0])) if not neighbor['flat'] else 0.
    before=candidate(baseline)[0];grid=np.linspace(-np.pi,np.pi,48,endpoint=False)
    angles=list(grid)+[baseline];values=[candidate(x)[0] for x in angles]
    best=int(np.argmin(values));theta=angles[best];value=values[best]
    for index in np.argsort(values)[:3]:
        angle=angles[int(index)];step=2*np.pi/48
        result=minimize_scalar(lambda x:candidate(x)[0],bounds=(angle-step,angle+step),method='bounded',options={'xatol':1e-9})
        if result.success and result.fun<value:theta=float(result.x);value=float(result.fun)
    refinement={}
    if refine:
        before_polish=value;start_angle=theta;iterations=0
        # Bounded minimization has an angle-dependent sqrt(eps) stopping
        # scale. Polish its selected basin using the residual vector instead
        # of differences of squared objectives near machine precision.
        for _ in range(12):
            r,_,_=residual(theta);h=1e-4
            derivative=(residual(theta+h)[0]-residual(theta-h)[0])/(2*h)
            norm=float(derivative@derivative)
            if norm<=1e-30:break
            correction=float(np.clip((r@derivative)/norm,-np.pi/48,np.pi/48))
            if abs(correction)<=8*np.finfo(float).eps*max(1.,abs(theta)):break
            improved=False
            for damping in (1.,.5,.25,.125):
                trial=theta-damping*correction;score=candidate(trial)[0]
                if score<value:
                    theta=trial;value=score;iterations+=1;improved=True;break
            if not improved:break
        refinement=dict(moment_fit_before_polish=before_polish,moment_polish_iterations=iterations,
            moment_polish_angle_change=float(theta-start_angle))
    value,normal,offset=candidate(theta)
    return dict(flat=False,normal=normal.tolist(),offset=offset,fraction=fraction,
        moment_fit_before=before,moment_fit_after=value,
        direction_change_degrees=float(np.degrees(np.arctan2(np.sin(theta-baseline),np.cos(theta-baseline)))),
        note='Mass fraction fixed; normalized q,Q squared residual with Frobenius weight. Finite angle search, no global optimality claim.',**refinement)


_MOMENT_CACHE={}
_SAVED_DIRECTIONS=None


def moment_prior(vertices,target,means,neighbor,refine=False):
    global _SAVED_DIRECTIONS
    if _SAVED_DIRECTIONS is None:
        cache=Path('moment_direction_cache.json')
        _SAVED_DIRECTIONS=json.loads(cache.read_text()) if cache.exists() else {}
    centered=os.environ.get('LAGFORCE_CENTERED_MOMENT_FIT','0')=='1'
    analytic=os.environ.get('LAGFORCE_ANALYTIC_CAP','0')=='1'
    if not centered and not analytic and not refine and target.tobytes().hex() in _SAVED_DIRECTIONS:
        return dict(_SAVED_DIRECTIONS[target.tobytes().hex()])
    key=(vertices.tobytes(),target.tobytes(),means.tobytes(),json.dumps(neighbor,sort_keys=True),refine,centered,analytic)
    if key not in _MOMENT_CACHE:_MOMENT_CACHE[key]=_fit_moment_prior(vertices,target,means,neighbor,refine=refine)
    return dict(_MOMENT_CACHE[key])


def aligned_parts(vertices,level,metadata):
    base=triangles(vertices,level)
    if metadata['flat']:return base
    normal=np.asarray(metadata['normal']);offset=metadata['offset'];parts=[]
    for triangle in base:
        for n,c in ((normal,offset),(-normal,-offset)):
            poly=clip(triangle,n,c)
            if poly.is_empty:continue
            v=np.asarray(poly.exterior.coords[:-1])
            for j in range(1,len(v)-1):
                part=np.array([v[0],v[j],v[j+1]])
                if triangle_moments(part)[0]>1e-25:parts.append(part)
    raw=sum((triangle_moments(v) for v in parts),np.zeros(7));original=polygon_moments(vertices)
    np.testing.assert_allclose(raw,original,rtol=1e-10,atol=1e-14)
    return parts


def main():
    previous=json.loads(Path('predecessor_completion.json').read_text())
    if not previous['completed'] or not all(c.get('all_gates') for c in previous['cases'].values()):raise RuntimeError('predecessor gate failure')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        baseline=dict(previous['cases'][name]);baseline['records']=[r for r in baseline['records'] if r['budget']==1e-5]
        baseline.update(reused=True)
        variants={'neighbor_uniform':baseline}
        for tag,prior_builder,parts_builder in (
            ('moment_uniform',moment_prior,None),('neighbor_aligned',None,aligned_parts),('moment_aligned',moment_prior,aligned_parts)):
            try:
                variants[tag]=selection_audit(name,prior_builder=prior_builder,parts_builder=parts_builder,
                    budgets=(1e-5,),output_prefix=name+'_'+tag)
            except Exception as exc:variants[tag]=dict(status='failed',all_gates=False,error=repr(exc))
            print(name,tag,variants[tag]['status'],flush=True)
        cases[name]=dict(status='completed' if all(v['status']=='completed' for v in variants.values()) else 'failed',
            all_gates=all(v['all_gates'] for v in variants.values()),variants=variants)
        atomic('partial_report.json',dict(cases=cases))
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),
        all_gates=all(c['all_gates'] for c in cases.values()),cases=cases,
        note='2x2 neighbor/moment normals and uniform/aligned basis. Levels2/3, same1% bounds, volume cap delta1e-5. Previous neighbor_uniform results reused. Exact references only score results. No evolution/HLL re-solve.')
    atomic('completion_report.json',result)
    raise SystemExit(0 if result['completed'] and result['all_gates'] else 2)


if __name__=='__main__':main()
