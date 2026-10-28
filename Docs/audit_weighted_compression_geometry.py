"""Static weighted-geometry response at a canonical C checkpoint.

Uses a regular-triangulation neighbor graph from the lower lifted hull, clips
only a bounded central patch, and checks the W=0 result against C geometry.
No fluid evolution, generator displacement, or particle mass change is made.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import lsmr


def clip(p, normal, bound):
    out=[]
    for a,b in zip(np.roll(p,1,axis=0),p):
        da=a@normal-bound;db=b@normal-bound
        if (da<=0)!=(db<=0):out.append(a+(b-a)*da/(da-db))
        if db<=0:out.append(b)
    return np.asarray(out).reshape(-1,2)


def geometry(x,w,v,ids,check_ids=()):
    n=len(x)
    lifted=np.column_stack((x,np.sum(x*x,axis=1)-w))
    hull=ConvexHull(lifted)
    tri=hull.simplices[hull.equations[:,2]<-1e-12]
    edges=np.concatenate((tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]]))
    edges=np.unique(np.sort(edges,axis=1),axis=0)
    neighbors=[[] for _ in range(n)]
    for i,j in edges:neighbors[i].append(j);neighbors[j].append(i)
    volumes=[];rates=[];row=[];col=[];val=[];cert=0.
    box=np.array([[-3.,-3.],[3.,-3.],[3.,3.],[-3.,3.]])
    for k,i in enumerate(ids):
        js=np.array(sorted(neighbors[i]),dtype=int)
        assert len(js)>0,('hidden cell',int(i))
        d=x[js]-x[i];bounds=np.sum(d*d,axis=1)+w[i]-w[js]
        p=box-x[i]
        for normal,bound in zip(2*d,bounds):
            p=clip(p,normal,bound)
            if len(p)<3:raise RuntimeError(f'empty cell {i}')
        q=np.roll(p,-1,axis=0)
        volume=.5*np.sum(p[:,0]*q[:,1]-p[:,1]*q[:,0])
        assert volume>0
        rate=0.;lapdiag=0.
        for aa,bb in zip(p,q):
            length=np.linalg.norm(bb-aa)
            if length<1e-13:continue
            mid=.5*(aa+bb)
            err=abs(2*d@mid-bounds)
            idx=int(np.argmin(err))
            assert err[idx]<1e-9,('unmatched or box face',int(i),float(err[idx]))
            j=js[idx];distance=np.linalg.norm(d[idx])
            speed=(d[idx]@v[j]-(v[j]-v[i])@mid)/distance
            rate+=length*speed
            c=length/(2*distance);lapdiag+=c
            row.append(k);col.append(j);val.append(-c)
        row.append(k);col.append(i);val.append(lapdiag)
        if i in check_ids:
            delta=x-x[i]
            bound=np.sum(delta*delta,axis=1)+w[i]-w
            violation=float(np.max(2*delta@p.T-bound[:,None]))
            cert=max(cert,violation)
            assert violation<1e-9,('all-site constraint violation',int(i),violation)
        volumes.append(volume);rates.append(rate)
    L=coo_matrix((val,(row,col)),shape=(len(ids),n)).tocsr()
    return np.array(volumes),np.array(rates),L,cert


def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args()
    args.output.mkdir(exist_ok=False)
    source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator')
    a=np.loadtxt(source/'rhs.dat');a=a[np.argsort(a[:,0])]
    x=a[:,1:3]-3;v=a[:,3:5];n=len(a);r=np.linalg.norm(x,axis=1)
    ids=np.flatnonzero(r<.70);unknown=np.flatnonzero(r<.50);w=np.zeros(n)
    check=set(ids[::max(1,len(ids)//16)].tolist())
    started=time.monotonic()
    V,D,L,cert=geometry(x,w,v,ids,check)
    parity=float(np.max(abs(V-a[ids,6])/a[ids,6]))
    assert parity<1e-7,('C volume mismatch',parity)
    f=np.loadtxt(source/'raw_faces.dat.common.rank0000')
    i,j=f[:,:2].astype(int).T;dist=np.linalg.norm(f[:,2:4],axis=1)
    mid=f[:,5:7]
    cflux=f[:,4]*(np.sum(f[:,2:4]*v[j],axis=1)-np.sum((v[j]-v[i])*mid,axis=1))/dist
    cD=np.bincount(i,weights=cflux,minlength=n)[ids]
    derivative_parity=float(np.max(abs(cD-D)/V))
    assert derivative_parity<1e-7,('C volume derivative mismatch',derivative_parity)
    eps=1e-5
    Vp,Dp,Lp,_=geometry(x+eps*v,w,v,ids)
    Vm,Dm,Lm,_=geometry(x-eps*v,w,v,ids)
    fd_error=float(np.max(abs((Vp-Vm)/(2*eps)-D)/V))
    assert fd_error<1e-5,('volume derivative finite difference',fd_error)
    # Mixed derivatives commute on a fixed topology: d_w(d_X V v)=d_X(V_w)v.
    J=(Lp-Lm)/(2*eps)
    Jscaled=diags(1/V)@J[:,unknown]@diags(a[unknown,6])
    target=-D/V
    reports=[]
    masks={'inner':r[ids]<.5,'guard':(r[ids]>=.5)&(r[ids]<.7),'all':np.ones(len(ids),bool)}
    def metrics(vol,rate):
        return {name:{'fractional_generator_volume_rate_rms':float(np.sqrt(np.average((rate[s]/vol[s])**2,weights=V[s]))),
                      'min_volume_over_baseline':float(np.min(vol[s]/V[s])),
                      'max_volume_over_baseline':float(np.max(vol[s]/V[s]))} for name,s in masks.items()}
    baseline=metrics(V,D)
    for damping in (.01,.1,1.):
        sol=lsmr(Jscaled,target,damp=damping,atol=1e-10,btol=1e-10,maxiter=2000)
        dw=np.zeros(n);dw[unknown]=a[unknown,6]*sol[0]
        maxsize=float(np.max(abs(sol[0])))
        cap=min(1.,.25/max(maxsize,1e-30))
        direction=dw*min(1.,.001/max(maxsize,1e-30))
        _,plus,_,_=geometry(x,direction,v,ids)
        _,minus,_,_=geometry(x,-direction,v,ids)
        predicted=J@direction
        mixed_relative=float(np.linalg.norm(.5*(plus-minus)-predicted)/max(np.linalg.norm(predicted),1e-30))
        for fraction in (.25,1.):
            trial=dw*cap*fraction
            vol,rate,_,certificate=geometry(x,trial,v,ids,check)
            name=f'damp{damping:g}_fraction{fraction:g}'
            np.savez(args.output/(name+'.npz'),weights=trial,ids=ids,volume=vol,generator_volume_rate=rate)
            report={'name':name,'damping':damping,'lsmr_stop':int(sol[1]),'iterations':int(sol[2]),
                    'linear_residual_norm':float(sol[3]),'mixed_derivative_directional_relative_error':mixed_relative,
                    'max_weight_over_initial_volume':float(np.max(abs(trial[unknown]/a[unknown,6]))),
                    'all_site_sample_max_violation':certificate,'regions':metrics(vol,rate)}
            reports.append(report)
            print(json.dumps(report),flush=True)
    result={'source':str(source),'source_sha256':hashlib.sha256((source/'rhs.dat').read_bytes()).hexdigest(),
            'scope':'Static geometry feasibility at original C positions and velocity. Only weights of current r<.5 cells varied; response measured through r<.7 including guard. Full site graph used. No mass assignment, hydrodynamics, time evolution or long-term claim.',
            'method':'Linearize generator-induced volume rate with respect to weights using mixed derivative of the geometric volume Jacobian. Regularized least squares then rebuild exact weighted geometry; bound trial |w|/V0 by .25.',
            'baseline_C_volume_max_relative_error':parity,'baseline_C_volume_rate_max_fractional_error':derivative_parity,
            'generator_volume_rate_fd_max_fractional_error':fd_error,'baseline_all_site_sample_max_violation':cert,
            'baseline':baseline,'trials':reports,'seconds':time.monotonic()-started,
            'limits':['Static optimization can exploit instantaneous topology and need not persist as rings slip.',
                      'No conservative weight/mass/energy coupling chosen; prior weight-only-remap obstruction still applies.',
                      'Mixed derivative loses smoothness at topology transitions; directional check reported rather than assumed.',
                      'All-site plane certificate is sampled; regular triangulation supplies constraints to every measured cell.']}
    (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'seconds':result['seconds'],'baseline':baseline,'output':str(args.output)}),flush=True)


if __name__=='__main__':main()
