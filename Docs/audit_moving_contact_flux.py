"""Manufactured exact contacts on all certified segments; not general Riemann states."""
from pathlib import Path
import json,importlib.util,numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_checkpoint_weight_remap_20261018'
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load('rd',base/'run_canonical_av_branches_20261018.py');a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');faces=json.loads((root/'certified_face_segments.json').read_text())['faces'];gamma=5/3
# Primitive contact: equal normal velocity/P; arbitrary density and tangent.
def flux(rhol,rhor,vl,vr,P,n,w0,w1,L,low,high,quadrature=False):
 un=float(vl@n);assert abs(vr@n-un)<1e-12;cuts=[low,high]
 if abs(w1)>1e-30:
  crossing=(un-w0)/w1
  if low<crossing<high:cuts.insert(1,crossing)
 F=np.zeros(4)
 for lo,hi in zip(cuts[:-1],cuts[1:]):
  mid=(lo+hi)/2;rho,v=(rhol,vl) if un-w0-w1*mid>=0 else (rhor,vr)
  E=P/(gamma-1)+.5*rho*(v@v)
  if quadrature:
   z,g=np.polynomial.legendre.leggauss(2);u=lo+(z+1)*.5*(hi-lo);ds=g*.5*(hi-lo)*L;rel=un-w0-w1*u
   F+=np.r_[rho*(ds@rel),rho*v*(ds@rel)+P*n*ds.sum(),E*(ds@rel)+P*un*ds.sum()]
  else:
   A=L*((un-w0)*(hi-lo)-.5*w1*(hi*hi-lo*lo));length=L*(hi-lo)
   F+=np.r_[rho*A,rho*v*A+P*n*length,E*A+P*un*length]
 return F,len(cuts)-1
rows=[]
for fraction in (.25,1.):
 maxquad=maxswap=maxboost=maxstationary=0.;count=crossings=0;tot=np.zeros(4);wrongwork=0.;b=np.array([.7,-.4])
 for face in faces:
  if face['fraction']!=fraction:continue
  n=np.array(face['normal']);tangent=np.array([-n[1],n[0]]);L=np.linalg.norm(np.array(face['b'])-face['a']);w0,w1=face['normal_speed_coefficients']
  for low,high,left,right in face['segments']:
   vl0=np.array([a['vx'][left],a['vy'][left]]);vr0=np.array([a['vx'][right],a['vy'][right]]);un=.5*((vl0+vr0)@n);vl=un*n+(vl0@tangent)*tangent;vr=un*n+(vr0@tangent)*tangent;P=.5*(a['pressure'][left]+a['pressure'][right]);rhol=float(a['den'][left]);rhor=float(a['den'][right])
   F,pieces=flux(rhol,rhor,vl,vr,P,n,w0,w1,L,low,high);Q,_=flux(rhol,rhor,vl,vr,P,n,w0,w1,L,low,high,True)
   reverse,_=flux(rhor,rhol,vr,vl,P,-n,-w0-w1,w1,L,1-high,1-low)
   boosted,_=flux(rhol,rhor,vl+b,vr+b,P,n,w0+b@n,w1,L,low,high)
   expected=np.r_[F[0],F[1:3]+b*F[0],F[3]+b@F[1:3]+.5*(b@b)*F[0]]
   stationary,_=flux(rhol,rhor,vl,vr,P,n,un,0.,L,low,high);exact=np.r_[0.,P*n*L*(high-low),P*un*L*(high-low)]
   maxquad=max(maxquad,float(abs(F-Q).max()));maxswap=max(maxswap,float(abs(F+reverse).max()));maxboost=max(maxboost,float(abs(boosted-expected).max()));maxstationary=max(maxstationary,float(abs(stationary-exact).max()));count+=1;crossings+=pieces-1
   # Omitting pressure work leaves this energy flux discrepancy in magnitude.
   wrongwork+=abs(P*un*L*(high-low));tot+=abs(F)
  # Opposite cell updates use the same owner F; orientation comparison above
  # independently evaluates reversed parameterization rather than just -F.
 rows.append({'fraction':fraction,'segments':count,'relative_speed_zero_crossings':crossings,'analytic_vs_split_Gauss_max_abs':maxquad,'orientation_reverse_max_abs':maxswap,'Galilean_flux_transform_max_abs':maxboost,'comoving_contact_max_abs':maxstationary,'sum_abs_contact_pressure_work':wrongwork,'sum_abs_flux':tot.tolist()})
result={'scope':__doc__,'gamma':gamma,'rows':rows,'limits':['Original normal velocities/pressures replaced by common contact values; density/tangent drawn from actual donor states.','Constant segment primitives; this does not integrate the previous fitted density or solve general jumps.','Exact contact only; no shock/rarefaction, full moment flux, force update, MPI or CUDA parity.','Common contact speed is not generally the geometric face speed; zero mass flux is tested only when wn=un.']}
with (root/'moving_contact_flux.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
