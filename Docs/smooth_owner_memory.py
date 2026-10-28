"""Experimental smooth, persistent-owner residual; no hydro evolution.

W_j=(1-d_j^2/R^2)_+^4, sum W=12 including owner. Quadratic
weighted least-squares prediction at owner, solved by SVD.
"""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def build_row(x, owner, tree=None, box=6., target=12.):
    if tree is None:tree=cKDTree(x % box,boxsize=box)
    k=min(128,len(x))
    while True:
        distances,ids=tree.query(x[owner] % box,k=k)
        distances=np.atleast_1d(distances);ids=np.atleast_1d(ids)
        assert distances[0]==0 and ids[0]==owner
        upper=min(float(distances[-1]),box/2)
        assert upper>0
        def weights(radius):return np.maximum(1-(distances/radius)**2,0)**4
        if weights(upper).sum()>target:break
        if k==len(x) or upper==box/2:raise RuntimeError('Uncovered or excessive support radius')
        k=min(2*k,len(x))
    lower=0.
    for _ in range(55):
        middle=(lower+upper)/2
        if weights(middle).sum()<target:lower=middle
        else:upper=middle
    radius=(lower+upper)/2;ww=weights(radius);keep=ww>0
    ids=ids[keep];ww=ww[keep]
    assert len(ids)>12
    dx=x[ids]-x[owner];dx-=box*np.rint(dx/box);dx/=radius
    P=np.column_stack((np.ones(len(ids)),dx,dx[:,0]**2,dx[:,0]*dx[:,1],dx[:,1]**2))
    A=np.sqrt(ww)[:,None]*P
    u,s,vt=np.linalg.svd(A,full_matrices=False)
    if s[-1]<=1e-10*s[0]:raise RuntimeError('Rank/conditioning failure')
    p0=np.array([1.,0,0,0,0,0])
    prediction=np.sqrt(ww)*(u@((vt@p0)/s))
    residual=-prediction;residual[np.where(ids==owner)[0][0]]+=1
    moments=P.T@residual
    return ids,residual,{'radius':radius,'support':len(ids),'condition':float(s[0]/s[-1]),
                        'moment_max':float(np.max(abs(moments))),
                        'weight_sum_error':float(abs(ww.sum()-target)),
                        'coverage_ratio':float(distances[-1]/radius)}


def audit():
    path=Path('/gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018/adaptive/raw_faces.dat.hgstate.rank0000')
    s=np.loadtxt(path);s=s[np.argsort(s[:,0])];x=s[:,3:5];tree=cKDTree(x%6,boxsize=6)
    owners=np.where(np.linalg.norm(x-3,axis=1)<.6)[0]
    stats=[]
    for owner in owners:
        _,_,meta=build_row(x,owner,tree);stats.append(meta)
    summary={key:{'min':min(row[key] for row in stats),'max':max(row[key] for row in stats)} for key in stats[0]}
    assert summary['moment_max']['max']<1e-12 and summary['weight_sum_error']['max']<1e-11
    crossings=[]
    for owner in [154,539,848]:
        ids,r,meta=build_row(x,owner,tree);R=meta['radius']
        dd,nn=tree.query(x[owner]%6,k=256)
        entering=int(nn[np.where(dd>R*1.01)[0][0]])
        direction=x[entering]-x[owner];direction-=6*np.rint(direction/6);direction/=np.linalg.norm(direction)
        at=x.copy();at[entering]=x[owner]+R*direction
        ids0,r0,m0=build_row(at,owner);dense=np.zeros(len(x));dense[ids0]=r0
        samples=[]
        for eps in [.02,.01,.005,.0025]:
            row={'relative_displacement':eps}
            for label,sign in [('inside',-1),('outside',1)]:
                moved=at.copy();moved[entering]=x[owner]+R*(1+sign*eps)*direction
                ii,rr,mm=build_row(moved,owner);dr=-dense.copy();dr[ii]+=rr
                coefficient=float(rr[np.where(ii==entering)[0][0]]) if entering in ii else 0.
                row[label]={'coefficient_change_L2':float(np.linalg.norm(dr)),
                            'entering_coefficient':coefficient,
                            'relative_radius_change':float((mm['radius']-m0['radius'])/R),
                            'moment_max':mm['moment_max']}
            samples.append(row)
        changes=[a['inside']['coefficient_change_L2'] for a in samples]
        assert all(changes[j+1]<changes[j]/8 for j in range(len(changes)-1))
        assert max(a['outside']['coefficient_change_L2'] for a in samples)<1e-12
        crossings.append({'owner':owner,'entering_particle':entering,'samples':samples})
    print(json.dumps({'scope':__doc__,'owners':len(owners),'static_summary':summary,
                      'crossings':crossings,
                      'limitations':'Static operator and three constructed crossings only. No full moving geometry, force-step integration, hydro, Laguerre, MPI or CUDA validation.'},indent=2))


if __name__=='__main__':audit()
