"""Conservative quadratic velocity on actual discontinuous-density IC moments.
Initial face accuracy and independent moment integration, not time evolution.
"""
import importlib.util,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'))
import moment_ale_evolve as model
import quadratic_velocity_reconstruction as quad
from kinetic_factor_rhs import ALPHA
from euler_star_pressure import star
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
data=np.load(base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz');ids=data['ids'];q=data['moments'];cross=data['crosses_inner_edge']
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');x=np.column_stack([ic['x']-3,ic['y']-3]);P=float(ic['pressure'][0]);gamma=float(np.float32(1.666667));model.GAMMA=gamma
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
os.environ.update(LAGFORCE_DENSITY_QUADRATURE_ORDER='12',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_BUDGET='projected')
states=[];polys=[];factors=[];lookup={int(i):k for k,i in enumerate(ids)}
for i,a in zip(ids,q):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]];vv=ff[:,6:8];poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    raw=np.r_[a[:3],a[3],a[4],a[4],a[5],a[6:12],a[12],a[13],a[13],a[14],a[15]]
    s=model.decode(poly,raw);states.append(s);polys.append(poly);factors.append(np.linalg.cholesky(s['R']))
root=base/'Codex_continuous_quadratic_velocity_20261018';root.mkdir(exist_ok=False)
results=[]
for mode in ('affine','projected'):
    diagnostics={} if mode=='affine' else quad.prepare(states,polys,[],q[:,0],factors)
    errors=[];gen=[]
    for i,a,s,p,L in zip(ids,q,states,polys,factors):
        pts,w=geom.quadrature(p,24);dm=w*s['density'].value(pts)
        samples,_=quad.velocity_samples(s,L,a[0],pts,ALPHA)
        mean=samples.mean(axis=1);momentum=dm@mean;D=pts.T@(dm[:,None]*mean)
        Pi=np.einsum('n,nki,nkj->ij',dm/4,samples,samples)
        errors.append(np.r_[momentum-a[6:8],(D-a[8:12].reshape(2,2)).ravel(),(Pi-np.array([[a[12],a[13]],[a[13],a[14]]])).ravel()])
        gen.append(float(np.linalg.norm(model.velocity(s,x[i])-model.affine_velocity(s,x[i]))))
    sums={name:np.zeros(4) for name in ('interface','other')};counts={name:0 for name in sums};failures=[]
    z,w=np.polynomial.legendre.leggauss(8);z=(z+1)/2;w/=2
    for row in faces:
        ii,jj=int(row[0]),int(row[1])
        if ii>=jj or ii not in lookup or jj not in lookup:continue
        i,j=lookup[ii],lookup[jj];a=x[ii]+row[6:8];b=x[ii]+row[8:10]
        if np.linalg.norm((a+b)/2)>.6:continue
        n=(x[jj]-x[ii])/np.linalg.norm(x[jj]-x[ii]);length=np.linalg.norm(b-a)
        key='interface' if cross[i] or cross[j] else 'other'
        for g,weight in zip(z,w):
            point=a+g*(b-a);vl=model.velocity(states[i],point);vr=model.velocity(states[j],point)
            r=np.linalg.norm(point);exact=np.sqrt(gravity_coefficient(np.array([r])))[0]*np.array([-point[1],point[0]])
            rl=float(states[i]['density'].value(point)[0]);rr=float(states[j]['density'].value(point)[0])
            try:ps,_=star(rl,vl,P,rr,vr,P,n,gamma)
            except Exception as exc:failures.append(dict(i=ii,j=jj,error=repr(exc)));continue
            sums[key]+=length*weight*np.array([1.,(ps/P-1)**2,.5*(np.sum((vl-exact)**2)+np.sum((vr-exact)**2)),((vr-vl)@n)**2]);counts[key]+=1
    e=np.array(errors)
    results.append(dict(mode=mode,diagnostics=diagnostics,maximum_independent_p_D_Pi_errors=np.max(abs(e),axis=0).tolist(),
        maximum_generator_velocity_change=max(gen),failures=failures,
        regions={name:dict(samples=counts[name],pressure_relative_RMS=float(np.sqrt(a[1]/a[0])),velocity_RMS=float(np.sqrt(a[2]/a[0])),normal_jump_RMS=float(np.sqrt(a[3]/a[0]))) for name,a in sums.items()}))
out=dict(scope=__doc__,results=results,limits='Existing smooth positive density with spatial moments, not the added-potential density. Independent volume quadrature24 versus density fit12. No mass solve, moving geometry, gravity, or stability claim. Analytic circular velocity only measures error.')
(root/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
