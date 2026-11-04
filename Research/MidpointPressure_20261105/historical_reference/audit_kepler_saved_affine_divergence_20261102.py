#!/usr/bin/env python3
"""Offline consistency probe of a reconstructed flux D and its exact adjoint.

This is NOT a native evolution operator, a selected method, or a release gate.
Both velocity-divergence and pressure-gradient consistency are measured to avoid
accepting a volume correction that merely transfers error to the force.
"""
from pathlib import Path
import hashlib,importlib.util,json,time
import numpy as np
from scipy.sparse import coo_matrix

t0=time.perf_counter()
spec=importlib.util.spec_from_file_location('m',Path(__file__).with_name('audit_kepler_weight_physical_comparison_20261101.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
R=m.ROOT;source=R/'native_u64_split_weight_screen_20261102_v1'
out=R/'native_u64_saved_affine_divergence_20261102_v1'
if out.exists():raise FileExistsError(out)
out.mkdir()
reg=json.loads((source/'registration.json').read_text())
raw,ids,_=m.load(Path(reg['source']));initial,iids,_=m.load(R/'native_u64_w0_initial_20261101_v1/initial696.bin')
assert np.array_equal(ids,iids)
ct=np.dtype([('id','<i8')]+[(k,'<f8')for k in ['W','rate','mu','Q','Edot','Hrow','Hcol','U','Udot','V']])
ft=np.dtype([('i','<i4'),('j','<i4')]+[(k,'<f8')for k in ['ax','ay','bx','by']])
c=np.fromfile(source/'run/weight_native_cells.bin',dtype=ct)
f=np.fromfile(source/'run/weight_native_faces.bin',dtype=ft)
order=np.searchsorted(ids,c['id']);assert np.array_equal(ids[order],c['id'])
x=m.column(raw,24)[order]-3;y=m.column(raw,32)[order]-3
vx=m.column(raw,56)[order];vy=m.column(raw,64)[order]
mass=m.column(raw,48)[order];V=c['V'];n=len(c)
r0=np.hypot(m.column(initial,24)[order]-3,m.column(initial,32)[order]-3)
real=f['j']>=0;i=f['i'][real].astype(int);j=f['j'][real].astype(int)
dx=x[j]-x[i];dy=y[j]-y[i];d2=dx*dx+dy*dy
count=lambda values:np.bincount(i,values,minlength=n)
Mxx=count(dx*dx/d2);Mxy=count(dx*dy/d2);Myy=count(dy*dy/d2)
det=Mxx*Myy-Mxy*Mxy
if np.any(det<=0):raise ValueError('Unregularized least-squares stencil singular')
cx=(Myy[i]*dx/d2-Mxy[i]*dy/d2)/det[i]
cy=(Mxx[i]*dy/d2-Mxy[i]*dx/d2)/det[i]
rows=np.concatenate((i,np.arange(n)));cols=np.concatenate((j,np.arange(n)))
Cx=coo_matrix((np.concatenate((cx,-count(cx))),(rows,cols)),shape=(n,n)).tocsr()
Cy=coo_matrix((np.concatenate((cy,-count(cy))),(rows,cols)),shape=(n,n)).tocsr()
ex=f['bx'][real]-f['ax'][real];ey=f['by'][real]-f['ay'][real]
fx=.5*(f['ax'][real]+f['bx'][real]);fy=.5*(f['ay'][real]+f['by'][real])
gx=fx-dx;gy=fy-dy
length=np.hypot(ex,ey);distance=np.sqrt(d2)
bi_x=.5*ey+length*(fx-dx/2)/distance
bi_y=-.5*ex+length*(fy-dy/2)/distance
bj_x=.5*ey-length*(fx-dx/2)/distance
bj_y=-.5*ex-length*(fy-dy/2)/distance
def B(ux,uy):return count(bi_x*ux[i]+bi_y*uy[i]+bj_x*ux[j]+bj_y*uy[j])
def BT(p):
    return(count(p[i]*bi_x)+np.bincount(j,p[i]*bj_x,minlength=n),count(p[i]*bi_y)+np.bincount(j,p[i]*bj_y,minlength=n))
def D(ux,uy):
    xx=Cx@ux;xy=Cy@ux;yx=Cx@uy;yy=Cy@uy
    face_x=.5*(ux[i]+xx[i]*fx+xy[i]*fy+ux[j]+xx[j]*gx+xy[j]*gy)
    face_y=.5*(uy[i]+yx[i]*fx+yy[i]*fy+uy[j]+yx[j]*gx+yy[j]*gy)
    return count(ey*face_x-ex*face_y)
def DT(p):
    force=[]
    for normal in [ey,-ex]:
        a=p[i]*normal/2
        plain=count(a)+np.bincount(j,a,minlength=n)
        dualx=count(a*fx)+np.bincount(j,a*gx,minlength=n)
        dualy=count(a*fy)+np.bincount(j,a*gy,minlength=n)
        force.append(plain+Cx.T@dualx+Cy.T@dualy)
    return tuple(force)
wall=np.zeros(n,dtype=bool);wall[f['i'][~real]]=True
interior=(r0<=2)&(~wall)
one=np.ones(n);zero=np.zeros(n)
patches={}
for name,ux,uy,expected in [('translation_x',one,zero,zero),('translation_y',zero,one,zero),('rotation',-y,x,zero),('pure_strain',x,-y,zero),('simple_shear',y,zero,zero),('x_dilation',x,zero,one),('isotropic_dilation',x,y,2*one)]:
    b=B(ux,uy)/V;d=D(ux,uy)/V
    patches[name]=dict(B_max_abs_divergence_error=float(abs(b[interior]-expected[interior]).max()),D_max_abs_divergence_error=float(abs(d[interior]-expected[interior]).max()))
omega=(x*x+y*y+.01175**2)**(-.75)
bcold=B(-omega*y,omega*x);dcold=D(-omega*y,omega*x)
bactual=B(vx,vy);dactual=D(vx,vy)
pressure=(m.GAMMA-1)*c['U']/V
force=DT(pressure)
work1=np.sum(vx*force[0]+vy*force[1],dtype=np.longdouble)
work2=np.sum(pressure*dactual,dtype=np.longdouble)
workscale=np.sum(abs(vx*force[0])+abs(vy*force[1])+abs(pressure*dactual),dtype=np.longdouble)
constant_force=DT(one)
# Exact polygon centroid of the stored piecewise-constant cell, for a second
# linear-pressure probe using a cell average rather than a site sample.
cross=f['ax']*f['by']-f['ay']*f['bx']
area=np.bincount(f['i'],cross/2,minlength=n)
center_x=x+np.bincount(f['i'],(f['ax']+f['bx'])*cross,minlength=n)/(6*area)
center_y=y+np.bincount(f['i'],(f['ay']+f['by'])*cross,minlength=n)/(6*area)
gradient={}
for label,px in [('site',x),('cell_average',center_x)]:
    p=1+.1*px
    gradient[label]={'B':BT(p),'D':DT(p)}
groups={}
for name,mask in [('core',r0<.1),('inner',(r0>=.1)&(r0<.5)),('main',(r0>=.5)&(r0<=2))]:
    weights=mass[mask].astype(np.longdouble);weights/=weights.sum()
    mean=lambda a:float(np.sum(weights*a[mask],dtype=np.longdouble))
    Bmean=mean(abs(bcold/V));Dmean=mean(abs(dcold/V))
    g={}
    for label,operators in gradient.items():
        g[label]={}
        for op,(ax,ay)in operators.items():
            g[label][op]=dict(relative_linear_pressure_force_error_L1=mean(np.hypot(ax+.1*V,ay)/(.1*V)),relative_force_x_error_L1=mean(abs(ax/(.1*V)+1)),relative_force_y_error_L1=mean(abs(ay/(.1*V))))
    groups[name]=dict(cells=int(mask.sum()),B_circular_abs_divergence_mean=Bmean,D_circular_abs_divergence_mean=Dmean,D_to_B_divergence_ratio=Dmean/Bmean,B_actual_abs_divergence_mean=mean(abs(bactual/V)),D_actual_abs_divergence_mean=mean(abs(dactual/V)),centroid_offset_over_sqrtV_mean=mean(np.hypot(center_x-x,center_y-y)/np.sqrt(V)),linear_pressure_patch=g)
checks=dict(affine_divergence_patches=all(v['D_max_abs_divergence_error']<1e-9 for v in patches.values()),work_exact_adjoint=bool(abs(work1-work2)<512*np.finfo(float).eps*max(workscale,1e-300)),constant_pressure_null=bool(max(abs(constant_force[0]).max(),abs(constant_force[1]).max())<1e-12),global_flux_conservative=bool(abs(np.sum(dactual,dtype=np.longdouble))<1e-12),polygon_area_matches_native=bool(abs(area/V-1).max()<1e-10))
report=dict(algebra_checks_passed=all(checks.values()),checks=checks,selected_for_evolution=False,saved_time=1.0000486862349771,patches=patches,cohorts=groups,work_duality_absolute=float(abs(work1-work2)),work_duality_relative=float(abs(work1-work2)/max(workscale,1e-300)),constant_pressure_force_max=float(max(abs(constant_force[0]).max(),abs(constant_force[1]).max())),global_material_volume_rate=float(np.sum(dactual,dtype=np.longdouble)),global_kinematic_volume_rate=float(np.sum(bactual,dtype=np.longdouble)),gradient_stencil_min_det=float(det.min()),elapsed_CPU_diagnostic_seconds=time.perf_counter()-t0,evidence=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in [Path(reg['source']),source/'run/weight_native_cells.bin',source/'run/weight_native_faces.bin']],interpretation='D reconstructs a common physical face velocity from unregularized affine least-squares gradients. DT is its exact algebraic transpose. Passing affine divergence or work identities alone does NOT establish pressure-gradient accuracy, Euler consistency, stability, positive volumes, viable weighted geometry or native performance. Both site-sampled and exact cell-average linear pressure tests are reported.',limitations=['Read-only CPU evaluation of saved native data; no Python or native time integration, new mesh, force modification or live-state change.','No weight solve, kinematic closure or spin model has been implemented or selected.','Linear pressure probes are consistency measurements, not replacements for accepted entropy/EOS labels.','CPU diagnostic elapsed time is not an evolution speed benchmark.'])
(out/'affine_divergence_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if not report['algebra_checks_passed']:raise SystemExit('Offline operator algebra failed; do not use as evidence of a valid replacement')
