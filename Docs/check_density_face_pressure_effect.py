"""Isolate density-closure effects on t0 face pressure at identical velocities.
Analytic density/velocity are diagnostic controls only, never evolution inputs.
"""
import importlib.util,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
import moment_ale_evolve as model
from potential_moment_density import PotentialDensity
from euler_star_pressure import star
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
data=np.load(base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz');ids=data['ids'];q=data['moments'];cross=data['crosses_inner_edge']
co=np.load(base/'Codex_potential_moment_density_20261018/coefficients.npz')['coefficients'];assert np.array_equal(co[:,0],ids)
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');x=np.column_stack([ic['x']-3,ic['y']-3]);P=float(ic['pressure'][0]);gamma=float(np.float32(1.666667));model.GAMMA=gamma
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']='12';states={};newdensity={};crossing={}
for i,a,c,cut in zip(ids,q,co,cross):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]];vv=ff[:,6:8];poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    raw=np.r_[a[:3],a[3],a[4],a[4],a[5],a[6:12],a[12],a[13],a[13],a[14],a[15]]
    states[int(i)]=model.decode(poly,raw);crossing[int(i)]=bool(cut)
    d=PotentialDensity.__new__(PotentialDensity);d.origin=c[1:3];d.scale=c[3];d.beta=c[4:10];d.potential_scale=c[10];d.theta=c[11:17];d.logZ=c[17];d.mass=a[0];d.potential=geom.potential;newdensity[int(i)]=d
result=[]
for nq in (4,8):
    z,w=np.polynomial.legendre.leggauss(nq);z=(z+1)/2;w=w/2
    sums={key:np.zeros(7) for key in ('interface','other')};count={key:0 for key in sums};fail=[];nullmax=0.
    for row in faces:
        i,j=int(row[0]),int(row[1])
        if i>=j or i not in states or j not in states:continue
        a=x[i]+row[6:8];b=x[i]+row[8:10]
        if np.linalg.norm((a+b)/2)>.6:continue
        edge=b-a;n=(x[j]-x[i])/np.linalg.norm(x[j]-x[i]);cuts=[0.,1.]
        A=edge@edge;B=2*a@edge;C=a@a-.25;disc=B*B-4*A*C
        if disc>0:
            cuts.extend(float(t) for t in ((-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)) if 0<t<1)
        key='interface' if crossing[i] or crossing[j] else 'other'
        cuts=sorted(cuts)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            for g,weight in zip(z,w):
                point=a+(lo+(hi-lo)*g)*edge;scale=np.sqrt(A)*(hi-lo)*weight
                vl=model.velocity(states[i],point);vr=model.velocity(states[j],point)
                radius=np.linalg.norm(point);truth=.01 if radius<.5 else 1.
                densities=[(float(states[i]['density'].value(point)[0]),float(states[j]['density'].value(point)[0])),(float(newdensity[i].value(point)[0]),float(newdensity[j].value(point)[0])),(truth,truth)]
                try:
                    ps=np.array([star(rl,vl,P,rr,vr,P,n,gamma)[0] for rl,rr in densities]);error=ps/P-1
                    # Shared exact circular velocity removes velocity reconstruction error.
                    vc=np.sqrt(gravity_coefficient(np.array([radius])))[0]*np.array([-point[1],point[0]])
                    for rl,rr in densities[:2]:nullmax=max(nullmax,abs(star(rl,vc,P,rr,vc,P,n,gamma)[0]/P-1))
                    sums[key]+=scale*np.r_[1.,error*error,((ps[:2]-ps[2])/P)**2,abs((ps[1]-ps[0])/P)]
                    count[key]+=1
                except Exception as exc:fail.append(dict(i=i,j=j,error=repr(exc)))
    summary={}
    for key,s in sums.items():
        summary[key]=dict(samples=count[key],length=s[0],pressure_relative_rms=dict(zip(['old','potential','analytic_density_control'],np.sqrt(s[1:4]/s[0]).tolist())),density_only_pressure_difference_rms=dict(zip(['old','potential'],np.sqrt(s[4:6]/s[0]).tolist())),mean_absolute_old_new_pressure_difference=float(s[6]/s[0]))
    result.append(dict(order=nq,regions=summary,failures=fail,shared_velocity_pressure_null_error=nullmax))
print(json.dumps(dict(scope=__doc__,results=result,limits='Face pressure only; no mass-root solve, gravity or time evolution. Analytic density is a diagnostic control and not a proposed solver. Smaller initial pressure error alone cannot establish stability.'),indent=2))
