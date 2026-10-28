"""Decode continuous-IC moments and inspect face traces; no evolution."""
import argparse, hashlib, importlib.util, json, os, sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--density-order',type=int,default=6)
args=parser.parse_args()
sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'))
import moment_ale_evolve as model
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
source=base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz'
data=np.load(source);ids=data['ids'];q=data['moments'];cross=data['crosses_inner_edge']
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
x=np.column_stack([ic['x']-3,ic['y']-3]);faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0];lookup={int(k):n for n,k in enumerate(ids)}
suffix='' if args.density_order==6 else f'_q{args.density_order}'
root=base/f'Codex_continuous_initial_reconstruction{suffix}_20261018';root.mkdir(exist_ok=False)
states={};failures=[];minR=np.inf;fitmax=0.;densitycheck=[];faceerrors=[]
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']=str(args.density_order)
for local,i in enumerate(ids):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vv=ff[:,6:8];poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    a=q[local];raw=np.r_[a[:3],a[3],a[4],a[4],a[5],a[6:12],a[12],a[13],a[13],a[14],a[15]]
    try:s=model.decode(poly,raw)
    except Exception as exc:failures.append(dict(id=int(i),error=repr(exc)));continue
    states[int(i)]=s;minR=min(minR,float(np.linalg.eigvalsh(s['R']).min()));fitmax=max(fitmax,s['density'].residual)
    if cross[local]:
        pts,dv=geom.quadrature(poly,max(24,args.density_order+8));rho=s['density'].value(pts);dm=dv*rho
        c=a[1:3]/a[0];scale=np.sqrt(dv.sum());z=(pts-c)/scale
        # Independent higher-order integration of fitted density, not its fit rule.
        F=np.column_stack([np.ones(len(z)),z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])
        target=np.r_[a[0],0.,0.,(a[3]-a[1]**2/a[0])/scale**2,(a[4]-a[1]*a[2]/a[0])/scale**2,(a[5]-a[2]**2/a[0])/scale**2]
        densitycheck.append(float(np.max(abs(dm@F-target))/a[0]))
    if np.linalg.norm(x[i])<.6:
        for row in ff:
            point=x[i]+.5*(row[6:8]+row[8:10]);r=np.linalg.norm(point)
            exact=np.sqrt(gravity_coefficient(np.array([r])))[0]*np.array([-point[1],point[0]])
            err=np.linalg.norm(model.velocity(s,point)-exact)
            faceerrors.append((bool(cross[local]),float(err)))
jumps=[];vacuum=[]
for row in faces:
    i,j=int(row[0]),int(row[1])
    if i>=j or i not in states or j not in states:continue
    point=x[i]+.5*(row[6:8]+row[8:10])
    if np.linalg.norm(point)>.6:continue
    n=(x[j]-x[i])/np.linalg.norm(x[j]-x[i]);left,right=states[i],states[j]
    jump=float((model.velocity(right,point)-model.velocity(left,point))@n)
    rl=float(left['density'].value(point)[0]);rr=float(right['density'].value(point)[0])
    threshold=2*(np.sqrt(model.GAMMA*left['P']/rl)+np.sqrt(model.GAMMA*right['P']/rr))/(model.GAMMA-1)
    jumps.append(abs(jump))
    if jump>=threshold:vacuum.append(dict(i=i,j=j,jump=jump,threshold=float(threshold)))
out=dict(scope=__doc__,source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    density_order=args.density_order,independent_check_order=max(24,args.density_order+8),
    cells=len(ids),decoded=len(states),failures=failures,minimum_raw_residual_eigenvalue=minR,maximum_density_fit_residual=fitmax,
    crossing_density_independent_moment_relative_error_max=max(densitycheck),
    face_midpoint_velocity_error={str(c):dict(count=sum(k==c for k,e in faceerrors),maximum=max(e for k,e in faceerrors if k==c),rms=float(np.sqrt(np.mean([e*e for k,e in faceerrors if k==c])))) for c in (False,True)},
    shared_faces=len(jumps),maximum_absolute_normal_jump=max(jumps),vacuum_faces=vacuum,
    pressure_gamma=float(model.GAMMA),initial_gamma=float(np.float32(1.666667)),
    limits='Affine traces only, midpoint samples only, no Riemann flux/RHS/time evolution. Smooth positive density fit does not resolve the actual discontinuous contact inside a cell. R checked without clipping; no stable covariance factors built here.')
(root/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
