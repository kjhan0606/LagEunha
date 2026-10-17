"""Static translation/scale conditioning of the two-phase normal inversion.

The same physical corner cap is initialized in three coordinate frames.
Long-double analytic polynomial transforms are quantized once to float64;
the reference normal is used only to initialize and score these new fixtures.
"""
import json,os,time
import numpy as np
from adaptive_coupled_mass import atomic
from moment_aligned_density_audit import _fit_moment_prior


def triangle_raw(v):
    v=np.asarray(v,dtype=np.longdouble);edge1=v[1]-v[0];edge2=v[2]-v[0]
    area=abs(edge1[0]*edge2[1]-edge1[1]*edge2[0])/2;s=v.sum(0)
    return np.r_[area,area*s/3,(area*(v.T@v+np.outer(s,s))/12).ravel()]


def translated_raw(raw,origin,scale):
    raw=np.asarray(raw,dtype=np.longdouble);origin=np.asarray(origin,dtype=np.longdouble);scale=np.longdouble(scale)
    m=raw[0];q=raw[1:3];Q=raw[3:].reshape(2,2)
    return np.r_[scale**2*m,scale**2*(scale*q+m*origin),
        (scale**2*(scale**2*Q+scale*(np.outer(origin,q)+np.outer(q,origin))+m*np.outer(origin,origin))).ravel()]


def run():
    started=time.monotonic();theta=np.longdouble('.73');normal=np.array([np.cos(theta),np.sin(theta)])
    square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]],dtype=np.longdouble)
    geometry=triangle_raw(square[[0,1,2]])+triangle_raw(square[[0,2,3]])
    rows=[]
    for fraction in (1e-2,1e-5,1e-7,1e-9):
        delta=np.sqrt(2*np.longdouble(str(fraction))*normal.prod())
        triangle=np.array([[1,1],[1-delta/normal[0],1],[1,1-delta/normal[1]]])
        physical=geometry+99*triangle_raw(triangle)
        for frame,origin,scale in [('unit',[0.,0.],1.),('n8_cell',[.45,.6],.125),('shift_scale_stress',[10.,-7.],.01)]:
            vertices=np.asarray(square*np.longdouble(scale)+np.asarray(origin,dtype=np.longdouble),dtype=float)
            exact=translated_raw(physical,origin,scale);target=np.asarray(exact,dtype=float)
            for mode in ('global','centered'):
                os.environ['LAGFORCE_CENTERED_MOMENT_FIT']='1' if mode=='centered' else '0'
                result=_fit_moment_prior(vertices,target,np.array([1.,100.]),dict(flat=False,normal=[1.,0.]),refine=True)
                inferred=np.asarray(result.get('normal',[np.nan,np.nan]))
                angle=np.arctan2(inferred[1],inferred[0]);error=abs(np.arctan2(np.sin(angle-float(theta)),np.cos(angle-float(theta))))
                row=dict(fraction=fraction,frame=frame,mode=mode,normal_angle_error_radians=float(error),
                    objective=result.get('moment_fit_after'),flat=result['flat'],
                    max_raw_quantization_relative=float(np.max(abs(exact-target)/np.maximum(abs(exact),np.longdouble('1e-30')))))
                rows.append(row);atomic('conditioning_progress.json',dict(rows=rows))
                print('conditioning',fraction,frame,mode,error,flush=True)
    atomic('moment_storage_conditioning_report.json',dict(completed=True,rows=rows,seconds=time.monotonic()-started,
        note='New static synthetic two-phase fixtures only. Same shape is represented in unit, physical-cell and translated small-cell frames before float64 quantization. No evolved state is reinitialized, corrected or supplied with exact interface data. Long double is a diagnostic reference, not an interval proof; angular accuracy is measured, not assumed.'))


if __name__=='__main__':run()
