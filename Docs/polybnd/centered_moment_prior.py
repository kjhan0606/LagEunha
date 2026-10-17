"""Experimental local-coordinate evaluation of the existing planar prior fit.

Only numerical conditioning changes: the six-moment least-squares objective,
phase values, fixed area fraction and angular search remain the same in exact
arithmetic. No true interface or neighbor steering enters this evaluation.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def fit_centered_prior(vertices,target,means,neighbor,refine=False):
    from moment_aligned_density_audit import polygon_moments,cut_at_fraction
    from subcell_density_audit import normalized
    low=float(means.min());high=float(means.max())
    geo=polygon_moments(vertices);mass=target[0]
    fraction=float(np.clip((mass/geo[0]-low)/(high-low),0,1)) if high>low else 0.
    if high-low<1e-12 or fraction<=1e-14 or fraction>=1-1e-14:
        return dict(flat=True,moment_fit_skipped=True)
    center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);local=(vertices-center)/scale
    local_geo=polygon_moments(local)[[0,1,2,3,4,6]]
    # Transform the stored physical moments once. Evaluating each trial cut
    # directly in this frame avoids repeatedly subtracting large raw moments.
    extended=np.longdouble
    goal=np.asarray(normalized(np.asarray(target,dtype=extended),
        np.asarray(center,dtype=extended),extended(scale))/extended(scale)**2,dtype=extended)
    cap_goal=np.asarray((goal-extended(low)*local_geo)/extended(high-low),dtype=float)
    weight=np.array([1.,1.,1.,1.,np.sqrt(2.),1.])
    amplitude=(high-low)*scale**2/mass

    def residual(theta):
        normal=np.array([np.cos(theta),np.sin(theta)])
        offset,cut=cut_at_fraction(local,normal,fraction)
        raw=np.zeros(7) if cut.is_empty else polygon_moments(np.asarray(cut.exterior.coords[:-1]))
        moments=np.r_[raw[:3],raw[3],raw[4],raw[6]]
        return (moments-cap_goal)*weight,normal,float(normal@center+scale*offset)

    def candidate(theta):
        diff,normal,offset=residual(theta)
        return float(diff@diff),normal,offset

    baseline=float(np.arctan2(neighbor['normal'][1],neighbor['normal'][0])) if not neighbor['flat'] else 0.
    before=candidate(baseline)[0]
    angles=list(np.linspace(-np.pi,np.pi,48,endpoint=False))+[baseline]
    values=[candidate(angle)[0] for angle in angles]
    best=int(np.argmin(values));theta=angles[best];value=values[best]
    for index in np.argsort(values)[:3]:
        angle=angles[int(index)];step=2*np.pi/48
        result=minimize_scalar(lambda x:candidate(x)[0],bounds=(angle-step,angle+step),
            method='bounded',options={'xatol':1e-9})
        if result.success and result.fun<value:theta=float(result.x);value=float(result.fun)
    refinement={}
    if refine:
        before_polish=value;start_angle=theta;iterations=0
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
        refinement=dict(moment_fit_before_polish=before_polish*amplitude**2,
            moment_polish_iterations=iterations,moment_polish_angle_change=float(theta-start_angle))
    value,normal,offset=candidate(theta)
    return dict(flat=False,normal=normal.tolist(),offset=offset,fraction=fraction,
        moment_fit_before=before*amplitude**2,moment_fit_after=value*amplitude**2,
        direction_change_degrees=float(np.degrees(np.arctan2(np.sin(theta-baseline),np.cos(theta-baseline)))),
        centered_moment_fit=True,
        note='Same six-moment objective evaluated in cell coordinates; finite angular search, no global optimality claim.',**refinement)
