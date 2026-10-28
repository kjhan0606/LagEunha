"""Compact piecewise donor face representation: segment count and line quadrature."""
import json, importlib.util, hashlib
from pathlib import Path
import numpy as np

def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_checkpoint_weight_remap_20261018'
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py')
reader=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=reader.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']))
lift=np.load(root/'lift_coefficients.npz');ids=lift['ids'];rho=a['mass'][ids]/lift['volume']
lam=lift['coefficient'];coeff=lam/lift['feature_scale'];means=lift['cell_feature_mean']@lam
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old])
normals=[np.column_stack(((np.roll(p,-1,axis=0)-p)[:,1],-(np.roll(p,-1,axis=0)-p)[:,0])) for p in old]
target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
result={'source_time':t,'scope':__doc__,'frame':'Each receiver mean velocity is the fixed diagnostic face velocity. Not the true evolving Laguerre face speed.','cases':{}}
for fraction in (.25,1.):
    saved=np.load(root/f'state_w{fraction:g}_q8.npz');q=saved['state'];assert np.array_equal(ids,saved['ids'])
    mass=q[:,0];c=q[:,1:3]/mass[:,None];v=q[:,6:8]/mass[:,None]
    S=np.array([[q[:,3],q[:,4]],[q[:,4],q[:,5]]]).transpose(2,0,1)-mass[:,None,None]*c[:,:,None]*c[:,None,:]
    C=q[:,8:12].reshape(-1,2,2)-mass[:,None,None]*c[:,:,None]*v[:,None,:]
    G=np.linalg.solve(S,C).transpose(0,2,1)
    new=mod.polygons(x,fraction*target,ids,geo.clip);records=[];skipped=0;eligible=0;minrho=float('inf')
    for k,p in enumerate(new):
        if np.linalg.norm(x[ids[k]])>=.6:continue
        for aa,bb in zip(p,np.roll(p,-1,axis=0)):
            d=bb-aa;length=np.linalg.norm(d)
            if length<1e-10:continue
            eligible+=1;n=np.array([d[1],-d[0]])/length
            candidates=np.flatnonzero(np.all(hi>=np.minimum(aa,bb)-1e-12,axis=1)&np.all(lo<=np.maximum(aa,bb)+1e-12,axis=1))
            intervals=[]
            for src in candidates:
                low,high=0.,1.
                for normal,point in zip(normals[src],old[src]):
                    den=normal@d;rhs=normal@(point-aa)
                    if abs(den)<1e-18:
                        if rhs < -1e-14: high=-1.;break
                    elif den>0:high=min(high,rhs/den)
                    else:low=max(low,rhs/den)
                if high-low>1e-10:intervals.append((low,high,src))
            # Reject coincident donor-boundary traces or incomplete partition.
            intervals.sort();cursor=0.;ok=True
            for low,high,src in intervals:
                if abs(low-cursor)>1e-8:ok=False;break
                cursor=high
            if not ok or abs(cursor-1)>1e-8:skipped+=1;continue
            flux=[]
            for order in (2,4,8,12):
                z,w=np.polynomial.legendre.leggauss(order);pair=np.zeros((2,4))
                for low,high,src in intervals:
                    u=low+(z+1)*.5*(high-low);xy=aa+u[:,None]*d;weight=w*.5*(high-low)*length
                    rel=xy-x[ids[src]];donor=vel[ids[src]]
                    f=rel@coeff[:2]+(rel[:,0]*donor[1]-rel[:,1]*donor[0])*coeff[2]+(mod.potential(xy)-mod.potential(x[ids[src]][None,:])[0])*coeff[3]-means[src]
                    density=rho[src]*(1+rho[src]*f);minrho=min(minrho,float(density.min()));assert density.min()>0
                    affine=v[k]+(xy-c[k])@G[k].T
                    for j,vv in enumerate((np.broadcast_to(donor,xy.shape),affine)):
                        relative=vv-v[k];mf=density*(relative@n)*weight
                        # Relative momentum and kinetic flux in the same frozen local frame.
                        pair[j]+=np.r_[mf.sum(),mf@relative,.5*np.sum(mf*np.sum(relative**2,axis=1))]
                flux.append(pair)
            records.append((k,length,flux[0],flux[1],flux[2],flux[3],len(intervals)))
    row={'eligible_faces':eligible,'retained_faces':len(records),'excluded_coincident_or_uncovered_faces':skipped,'minimum_sampled_density':minrho,'regions':{}}
    for label,test in [('inner',lambda r:r<.5),('edge',lambda r:.5<=r<.6)]:
        rec=[r for r in records if test(np.linalg.norm(x[ids[r[0]]]))]
        actual=np.array([r[5][0] for r in rec]);counts=np.array([r[6] for r in rec]);checks={}
        for order,idx in ((2,2),(4,3),(8,4)):
            approx=np.array([r[idx][0] for r in rec]);diff=approx-actual
            checks[str(order)]={'max_absolute_difference':float(abs(diff).max()),'relative_L2_difference':(np.linalg.norm(diff,axis=0)/np.maximum(np.linalg.norm(actual,axis=0),1e-300)).tolist()}
        row['regions'][label]={'faces':len(rec),'segments':int(counts.sum()),'segments_mean':float(counts.mean()),'segments_max':int(counts.max()),'segment_count_histogram':{str(i):int(np.sum(counts==i)) for i in np.unique(counts)},'quadrature_checks_vs_order12':checks,'components':['mass','relative_px','relative_py','relative_K']}
    result['cases'][str(fraction)]=row
result['limitations']=['One-sided receiver reconstruction; no common Riemann flux, pressure work, gravity, time update or stability claim.','Same donor density used for both velocity traces to isolate velocity reconstruction.','Coincident donor-boundary traces excluded; not a full closed-cell balance.','Relative kinetic flux is measured in receiver-mean frame, not lab-frame total energy.','Constant donor velocity and fixed diagnostic face speed allow one density line integral per segment to generate all four flux components.','Segment discovery and evolving provenance are not eliminated by this quadrature comparison.']
result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
with (root/'compact_face_trace_audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
