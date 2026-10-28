"""Split the circular density interface explicitly; preserve prior quadratures."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec

base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
old=base/'Codex_continuous_kepler_initial_moments_20261018'
report=json.loads((old/'report.json').read_text()); data=np.load(old/'moments_q32.npz')
ids=data['ids'];cross=data['crosses_inner_edge']; original=data['moments']
ic,t=read(Path(report['source']['source']));assert t==0
x=np.column_stack([ic['x']-3,ic['y']-3]);P=float(ic['pressure'][0]);gamma=float(np.float32(1.666667))
env=report['environment'];rin=float(env['EUNHA_KEPLER_RIN']);floor=float(env['EUNHA_KEPLER_FLOOR'])
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
root=base/'Codex_continuous_kepler_interface_split_20261018';root.mkdir(exist_ok=False)
results=[]; previous=None
for degree in (8,16):
    z,w=np.polynomial.legendre.leggauss(degree);moments=original.copy(); errors=[]
    for local in np.flatnonzero(cross):
        i=ids[local]
        rows=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
        v=rows[:,6:8];poly=v[np.argsort(np.arctan2(v[:,1],v[:,0]))]+x[i]
        aa=poly;dd=np.roll(poly,-1,axis=0)-poly
        knots=list(poly[:,0])
        for a,d in zip(aa,dd):
            A=d@d;B=2*a@d;C=a@a-rin*rin;disc=B*B-4*A*C
            if disc>=0:
                for s in ((-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)):
                    if 0<s<1:knots.append(float(a[0]+s*d[0]))
        lo=float(np.min(poly[:,0]));hi=float(np.max(poly[:,0]))
        knots.extend(k for k in (-rin,rin) if lo<k<hi)
        knots=np.unique(knots)
        def vertical(xx):
            valid=(abs(dd[:,0])>1e-16)
            s=np.zeros(len(dd));s[valid]=(xx-aa[valid,0])/dd[valid,0]
            valid &= (s>=-1e-12)&(s<=1+1e-12)
            yy=aa[valid,1]+s[valid]*dd[valid,1]
            if len(yy)<2:return np.zeros(17)
            bottom=float(yy.min());top=float(yy.max());cuts=[bottom,top]
            if abs(xx)<rin:
                h=np.sqrt(rin*rin-xx*xx)
                cuts.extend(y for y in (-h,h) if bottom<y<top)
            cuts=sorted(cuts);out=np.zeros(17)
            for yl,yh in zip(cuts[:-1],cuts[1:]):
                y=.5*(yh+yl)+.5*(yh-yl)*z;dv=.5*(yh-yl)*w
                pts=np.column_stack([np.full(len(y),xx),y]);r=np.linalg.norm(pts,axis=1)
                density=floor if xx*xx+(.5*(yl+yh))**2<rin*rin else 1.
                dm=density*dv;vel=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack([-y,np.full(len(y),xx)])
                q=dm@pts;Q=pts.T@(dm[:,None]*pts);p=dm@vel;D=pts.T@(dm[:,None]*vel);Pi=vel.T@(dm[:,None]*vel)
                out+=np.r_[dm.sum(),q,Q[0,0],Q[0,1],Q[1,1],p,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],P*dv.sum()/(gamma-1)+.5*np.trace(Pi),dm@geom.potential(pts)]
            return out
        val,err=quad_vec(vertical,lo,hi,points=knots[1:-1],epsabs=1e-13,epsrel=1e-11,limit=1000)
        assert np.isfinite(val).all();moments[local]=val;errors.append(float(err))
    thermal=moments[:,15]-.5*(moments[:,12]+moments[:,14]);assert np.all(thermal>0)
    np.savez(root/f'moments_split_q{degree}.npz',ids=ids,moments=moments,crosses_inner_edge=cross)
    entry=dict(degree=degree,crossing_cells=int(cross.sum()),maximum_outer_error_estimate=max(errors),
        minimum_thermal=float(thermal.min()),mass_difference_from_point=float(np.sum(moments[:,0]-ic['mass'][ids])),
        max_relative_mass_change_from_unsplit_q32=float(np.max(abs(moments[cross,0]/original[cross,0]-1))))
    if previous is not None:
        entry['successive_max_absolute_moment_difference']=float(np.max(abs(moments-previous)))
        entry['successive_max_relative_mass_difference']=float(np.max(abs(moments[:,0]/previous[:,0]-1)))
    previous=moments;results.append(entry)
out=dict(scope=__doc__,source=report['source'],results=results,
    limits='Only the 101 crossing cells are reintegrated. Other cells retain q32. This is a new continuous-field projection, not the original checkpoint state. Density reconstruction across internal interface, force evolution, and long-time stability remain unverified.')
(root/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
