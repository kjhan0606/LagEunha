"""Frozen-Q2 affine mass, momentum, angular-momentum and energy geometry audit."""
import ctypes
import hashlib
import json
import socket
import sys
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ID_MASK, ordered_global
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_affine_spin_energy_20261027 import tensor_from_faces
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
late=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
prior=json.loads((docs/'KEPLER_Q2_NORMAL_TRACTION_2026-10-29.json').read_text())
root=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1')
initial=Path(late['initial']);path=root/'khout.011058.dat'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
a,ma=ordered_global(initial);p,mt=ordered_global(path,a['u4if']&ID_MASK)
ref=next(v for v in prior['cases'] if v['time']==mt['t'])
if sha(initial)!=late['initial_sha256'] or sha(path)!=ref['source_sha256'] or ma['t']!=0 or np.any(p['w2']!=0):raise ValueError('changed input')
xy0=np.column_stack((a['x']-3,a['y']-3));r0=np.linalg.norm(xy0,axis=1)
selected=np.flatnonzero((r0>=.5)&(r0<.75));flag=np.zeros(len(a),bool);flag[selected]=True
lookup=np.full(len(a),-1,int);lookup[selected]=np.arange(len(selected))
xy=np.column_stack((p['x']-3,p['y']-3));base=np.random.default_rng(20261030).normal(0,2e-6,len(a))
rates=np.random.default_rng(20261029).normal(0,1e-4,len(a))
def moments(weight,with_faces=False):
    area,info,tri,centers=power_areas(xy,weight,selected,lookup,return_geometry=True)
    if info['hidden_active'] or np.any(area<=0):raise ValueError('hidden or degenerate')
    owner,u,v,offset,counts=oriented_faces(tri,centers,xy,selected,lookup)
    geom,cl,qxx,qxy,qyy=tensor_from_faces(u,v,owner,offset,counts,area)
    if np.max(abs(geom/area-1))>1e-8:raise ValueError('polygon moment inconsistency')
    cx=xy[selected,0]+cl[:,0];cy=xy[selected,1]+cl[:,1]
    result=np.column_stack((area,area*cx,area*cy,qxx+area*cx*cx,qxy+area*cx*cy,qyy+area*cy*cy))
    if with_faces:return result,tri,centers
    return result
_,tri,centers=moments(base,True)
i,j,c,length=finite_faces_touching_active(tri,centers,flag)
dr=xy[j]-xy[i];distance=np.linalg.norm(dr,axis=1);tangent=np.column_stack((-dr[:,1],dr[:,0]))/distance[:,None]
start=c-.5*length[:,None]*tangent;end=c+.5*length[:,None]*tangent
conv=lambda v:np.ascontiguousarray(v,dtype=np.float64)
arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib=ctypes.CDLL(str(Path(sys.argv[1]).resolve()));fn=lib.affine_batch
fn.argtypes=[ctypes.c_int]+[arr]*9;fn.restype=ctypes.c_int
oi=np.empty(len(i)*13);oj=np.empty(len(i)*13)
if not fn(len(i),conv(rates[i]),conv(rates[j]),conv(distance),conv(start[:,0]),conv(start[:,1]),conv(end[:,0]),conv(end[:,1]),oi,oj):raise ValueError('candidate rejected')
oi=oi.reshape(-1,13);oj=oj.reshape(-1,13);deriv=np.zeros((len(selected),13))
pair_residual=float(np.max(np.abs(oi+oj)))
if pair_residual>1e-14:raise ValueError('shared-face conservation failed')
ii=flag[i];jj=flag[j];np.add.at(deriv,lookup[i[ii]],oi[ii]);np.add.at(deriv,lookup[j[jj]],oj[jj])
def affine_integrals(v):
    m,x,y,xx,xy,yy=v.T
    bx,by=.7,-.4
    a,b,c,d=.31,-.72,.45,-.12
    px=bx*m+a*x+b*y; py=by*m+c*x+d*y
    dxx=bx*x+a*xx+b*xy; dxy=by*x+c*xx+d*xy
    dyx=bx*y+a*xy+b*yy; dyy=by*y+c*xy+d*yy
    energy=2*m+.5*((bx*bx+by*by)*m+2*(bx*a+by*c)*x+
        2*(bx*b+by*d)*y+(a*a+c*c)*xx+2*(a*b+c*d)*xy+(b*b+d*d)*yy)
    return np.column_stack((v,px,py,dxx,dxy,dyx,dyy,energy))
rows=[]
for h in (1e-4,5e-5):
    numerical=(affine_integrals(moments(base+h*rates))-affine_integrals(moments(base-h*rates)))/(2*h)
    errors=numerical-deriv
    relative=[float(np.linalg.norm(errors[:,k])/np.linalg.norm(numerical[:,k])) for k in range(13)]
    if max(relative)>1e-7:raise ValueError('affine geometric derivative mismatch')
    rows.append({'step':h,'relative_rms_by_moment':relative,'max_abs_by_moment':np.max(abs(errors),axis=0).tolist()})
result={'passed':True,'pair_max_abs_residual':pair_residual,'snapshot_sha256':ref['source_sha256'],'snapshot_time':mt['t'],'sites':len(selected),'faces':len(i),
        'moments':['mass','x','y','xx','xy','yy','px','py','xpx','xpy','ypx','ypy','energy'],'steps':rows,'bridge_sha256':sha(Path(sys.argv[1])),
        'script_sha256':sha(Path(__file__)),
        'scope':'Frozen bounded Q2 geometry, fixed generators, synthetic signed weights, uniform density, shared affine velocity and constant specific internal energy. Centered differences of analytic polygon integrals against actual C two-point affine weight-advection operator; no native negative-weight support, moving-generator flux, Riemann pressure/work, physical closure, finite-time evolution, MPI or stability proof.'}
out=docs/'KEPLER_WEIGHT_AFFINE_GEOMETRY_2026-10-29.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
