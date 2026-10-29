"""Compare actual affine force with an analytic Kepler-gradient reference."""
import hashlib
import json
import socket
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK,ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces
from audit_kepler_t3_support_20261029 import orbit
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_MODE_TORQUE_2026-10-29.json').read_text())
initial=Path(late['initial']);root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
if sha(initial)!=late['initial_sha256'] or ma['t']!=0:raise ValueError('initial mismatch')
xy0=np.column_stack((a['x']-3,a['y']-3));r0=np.linalg.norm(xy0,axis=1)
active=np.flatnonzero(r0<2.35);inner=(r0>=.5)&(r0<.75)
flag=np.zeros(len(a),bool);flag[active]=True;lookup=np.full(len(a),-1,int);lookup[active]=np.arange(len(active))
mass=a['mass'].astype(float);zero=np.zeros(len(a));phi=np.arctan2(xy0[:,1],xy0[:,0]);rad=np.round(r0,10)
L0=xy0[:,0]*a['vy']-xy0[:,1]*a['vx']
def geometry(xy):
    area,g,tri,centers=power_areas(xy,zero,active,lookup,return_geometry=True)
    if g['hidden_active'] or np.any(area<=0):raise ValueError('invalid geometry')
    owner,u,v,offset,counts=oriented_faces(tri,centers,xy,active,lookup)
    area2,cl,qxx,qxy,qyy=tensor_from_faces(u,v,owner,offset,counts,mass[active])
    if np.max(abs(area2/area-1))>1e-8 or np.any(qxx*qyy-qxy*qxy<=0):raise ValueError('invalid moments')
    S=np.empty((len(active),2,2));S[:,0,0]=qxx;S[:,0,1]=S[:,1,0]=qxy;S[:,1,1]=qyy
    return area,cl,S,tri,centers
def force(xy,geo,A):
    area,cl,S,tri,centers=geo;G=A@np.linalg.inv(S)
    kin=np.einsum('nij,nij->n',G,A);center=np.zeros_like(xy);center[active]=xy[active]+cl
    rho=np.zeros(len(a));rho[active]=mass[active]/area
    Gfull=np.zeros((len(a),2,2));Gfull[active]=G
    kpm=np.zeros(len(a));kpm[active]=kin/mass[active]
    i,j,c,length=finite_faces_touching_active(tri,centers,flag)
    if np.any((inner[i]&~flag[j])|(inner[j]&~flag[i])):raise ValueError('missing inner neighbor')
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1);n=delta/dist[:,None]
    tang=np.column_stack((-n[:,1],n[:,0]));mid=.5*(xy[i]+xy[j]);F=np.zeros_like(xy)
    for sign in (-1,1):
        q=c+sign*(length/(2*np.sqrt(3)))[:,None]*tang
        ui=np.einsum('nij,nj->ni',Gfull[i],q-center[i]);uj=np.einsum('nij,nj->ni',Gfull[j],q-center[j])
        df=.5*rho[i]*(np.sum(ui*ui,axis=1)-kpm[i])-.5*rho[j]*(np.sum(uj*uj,axis=1)-kpm[j])
        off=(q-mid)/dist[:,None];np.add.at(F,i,.5*length[:,None]*df[:,None]*(.5*n+off));np.add.at(F,j,.5*length[:,None]*df[:,None]*(.5*n-off))
    return F
def projection(err,torque):
    power=0.;cross=0.
    for radius in np.unique(rad[inner]):
        ix=np.flatnonzero(inner&(rad==radius));ix=ix[np.argsort(phi[ix])];n=len(ix)
        U=np.fft.rfft(err[ix]-np.mean(err[ix]))/n;V=np.fft.rfft(torque[ix]-np.mean(torque[ix]))/n
        factor=np.full(len(U),2.);factor[0]=1.
        if n%2==0:factor[-1]=1.
        freq=np.fft.rfftfreq(n);band=(freq>=.2)&(freq<.3);weight=float(np.sum(mass[ix]))
        power+=weight*np.sum(factor[band]*abs(U[band])**2);cross+=weight*np.sum(factor[band]*np.real(np.conj(U[band])*V[band]))
    return float(cross/power)
cases=[]
for path in (Path(late['cases'][0]['endpoint']),root/'khout.011058.dat'):
    p,meta=ordered_global(path,ids);ref=next(v for v in prior['cases'] if v['time']==meta['t'])
    if sha(path)!=ref['source_sha256'] or np.any(p['w2']!=0) or not np.array_equal(p['mass'],mass):raise ValueError('changed state')
    xy=np.column_stack((p['x']-3,p['y']-3));exact=orbit(xy0+3,meta['t'])-3
    actual_geo=geometry(xy);G=np.stack((np.stack((p['stress.gUxx'][active],p['stress.gUxy'][active]),axis=1),np.stack((p['stress.gUyx'][active],p['stress.gUyy'][active]),axis=1)),axis=1)
    A=G@actual_geo[2]
    Fa=force(xy,actual_geo,A)
    exact_geo=geometry(exact);x=exact[active,0];y=exact[active,1]
    r2=x*x+y*y+.01175**2;omega=r2**(-.75);domega=-1.5*r2**(-1.75)
    Ge=np.empty_like(G);Ge[:,0,0]=-x*y*domega;Ge[:,0,1]=-omega-y*y*domega
    Ge[:,1,0]=omega+x*x*domega;Ge[:,1,1]=x*y*domega
    # Independent centered difference on a spread of inner-disk positions.
    sample=np.flatnonzero(inner[active])[::max(1,int(np.sum(inner)/32))][:32]
    def vel(q):
        om=(np.sum(q*q,axis=1)+.01175**2)**(-.75)
        return np.column_stack((-q[:,1]*om,q[:,0]*om))
    max_gradient_error=0.
    for axis in (0,1):
        shift=np.zeros((len(sample),2));shift[:,axis]=1e-6
        numerical=(vel(exact[active[sample]]+shift)-vel(exact[active[sample]]-shift))/(2e-6)
        max_gradient_error=max(max_gradient_error,float(np.max(abs(numerical-Ge[sample,:,axis]))))
    if max_gradient_error>1e-8:raise ValueError('analytic gradient mismatch')
    Fe=force(exact,exact_geo,Ge@exact_geo[2])
    ta=(xy[:,0]*Fa[:,1]-xy[:,1]*Fa[:,0])/mass
    te=(exact[:,0]*Fe[:,1]-exact[:,1]*Fe[:,0])/mass
    err=xy[:,0]*p['vy']-xy[:,1]*p['vx']-L0
    omega_gravity=np.linalg.norm(xy,axis=1)/(np.sum(xy*xy,axis=1)+.01175**2)**1.5
    m=mass[inner];ga=np.linalg.norm(Fa[inner],axis=1)/m;ge=np.linalg.norm(Fe[inner],axis=1)/m
    cases.append({'time':meta['t'],'snapshot_sha256':ref['source_sha256'],
                  'actual_torque_projection':projection(err,ta),
                  'exact_geometry_analytic_gradient_torque_projection':projection(err,te),
                  'actual_minus_reference_projection':projection(err,ta-te),
                  'actual_accel_over_gravity_mass_mean':float(np.average(ga/omega_gravity[inner],weights=m)),
                  'reference_accel_over_actual_gravity_mass_mean':float(np.average(ge/omega_gravity[inner],weights=m)),
                  'analytic_gradient_finite_difference_max_abs_error':max_gradient_error,
                  'reference_energy_model':'A is constructed from analytic softened-Kepler gradient and exact-position cell moment'})
result={'cases':cases,'initial_sha256':sha(initial),'script_sha256':sha(Path(__file__)),
        'scope':'Frozen counterfactual: exact circular positions and analytically reconstructed Kepler velocity gradients with uniform-density polygon moments. Not a time-evolved affine state; exterior exact motion does not reproduce the actual outer-boundary treatment. Inner-site neighbors are included. This is not a stability, native Euler, or conservation proof.'}
out=docs/'KEPLER_Q2_AFFINE_REFERENCE_ANALYTIC_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
