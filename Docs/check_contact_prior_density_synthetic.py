"""Analytic planar-contact fixtures, independent of the Kepler benchmark."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_contact_prior_density_modules_20261018')
from contact_prior_density import ContactPriorDensity
from density_quadrature import rule
poly=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
rows=[];rng=np.random.default_rng(4271)
for kind in ('uniform','planar_contact'):
    low,high=(2.,2.) if kind=='uniform' else (1.,4.)
    cut=.37;m=low*cut+high*(1-cut)
    qx=(low*cut**2+high*(1-cut**2))/2;q=np.array([qx,m/2])
    Q=np.array([[(low*cut**3+high*(1-cut**3))/3,qx/2],[qx/2,m/3]])
    for angle,shift in ((0.,np.zeros(2)),(.73,np.array([2.1,-1.7]))):
        R=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        points=poly@R.T+shift;rq=R@q
        tq=rq+m*shift;tQ=R@Q@R.T+np.outer(rq,shift)+np.outer(shift,rq)+m*np.outer(shift,shift)
        normal=R[:,0];meta=dict(flat=False,low=low,high=high,normal=normal,offset=cut+normal@shift)
        if kind=='uniform':meta=dict(flat=True,value=2.)
        d=ContactPriorDensity(points,m,tq,tQ,meta,12)
        sample=rng.uniform(size=(100,2));exact=np.where(sample[:,0]>=cut,high,low)
        error=float(abs(d.value(sample@R.T+shift)-exact).max())
        allpoints=[];allweights=[]
        for part in d.parts:
            p,w=rule(part,24);allpoints.append(p);allweights.append(w)
        p=np.concatenate(allpoints);w=np.concatenate(allweights)*d.value(p)
        recovered=np.r_[w.sum(),w@p,(p.T@(w[:,None]*p)).ravel()]
        target=np.r_[m,tq,tQ.ravel()]
        moment_error=float(abs(recovered-target).max())
        assert error<1e-10 and moment_error<1e-10
        assert np.all(d.value(p)>0)
        rows.append(dict(kind=kind,angle=angle,value_max_error=error,independent_moment_max_error=moment_error,fit_residual=d.residual))
print(json.dumps(dict(passed=True,results=rows,scope='Density fit and prescribed-cut rotation/translation only. Inferred cuts and real-cell candidate not yet validated.'),indent=2))
