"""Adaptive whole-polygon positivity bound for the saved L2 density lift."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    root=Path(__file__).parent;base=Path('/gpfs/kjhan/LagForce')
    ap=argparse.ArgumentParser()
    ap.add_argument('--candidate',type=Path,default=base/'Codex_weight_remap_initialization_check_v3_20261018/l2_candidate.npz')
    ap.add_argument('--checkpoint',type=Path)
    args=ap.parse_args()
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    candidate=args.candidate
    z=np.load(candidate);ids=z['ids'];lam=z['coefficient'];scale=z['feature_scale'];means=z['cell_feature_mean']
    if args.checkpoint:
        source=args.checkpoint
        reader=load('reader',base/'Codex_paired_response_20261017/run.py')
        layout=json.loads((base/'Codex_paired_response_20261017/layout.json').read_text())
        _,a,order,_=reader.read(source,layout);a=a[order]
        x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']));mass=a['mass']
    else:
        source=base/'Codex_gradient_operator_20261017/operator/rhs.dat'
        a=np.loadtxt(source);a=a[np.argsort(a[:,0])]
        x=a[:,1:3]-3;vel=a[:,3:5];mass=a[:,7]
    rho=mass[ids]/z['volume'];coeff=lam/scale
    polygons=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);phi0=mod.potential(x[ids])
    rows=[]
    for k,p in enumerate(polygons):
        v=vel[ids[k]];linear=np.array([coeff[0]+coeff[2]*v[1],coeff[1]-coeff[2]*v[0]])
        offset=1-rho[k]*(means[k]@lam)
        center=p.mean(axis=0);stack=[(np.array([center,aa,bb]),0) for aa,bb in zip(p,np.roll(p,-1,axis=0))]
        lo=np.inf;hi=-np.inf;maxdepth=0;leaves=0;unresolved=0
        while stack:
            t,depth=stack.pop();values=(t-x[ids[k]])@linear
            inside=True;distance=[]
            for aa,bb in zip(t,np.roll(t,-1,axis=0)):
                edge=bb-aa;inside &= (edge[0]*(-aa[1])-edge[1]*(-aa[0])>=0)
                s=np.clip(-(aa@edge)/(edge@edge),0,1) if edge@edge else 0.
                distance.append(np.linalg.norm(aa+s*edge))
            rmin=0. if inside else min(distance);rmax=np.max(np.linalg.norm(t,axis=1))
            phis=coeff[3]*(mod.potential(np.array([[rmin,0.],[rmax,0.]]))-phi0[k])
            lower=offset+rho[k]*(values.min()+phis.min())-1e-11
            upper=offset+rho[k]*(values.max()+phis.max())+1e-11
            if lower>0 or depth>=8:
                lo=min(lo,float(lower));hi=max(hi,float(upper));maxdepth=max(depth,maxdepth);leaves+=1
                unresolved+=int(lower<=0)
            else:
                aa,bb,cc=t;ab=(aa+bb)/2;bc=(bb+cc)/2;ca=(cc+aa)/2
                stack.extend((np.array(q),depth+1) for q in ((aa,ab,ca),(ab,bb,bc),(ca,bc,cc),(ab,bc,ca)))
        rows.append({'id':int(ids[k]),'density_factor_lower_bound':lo,'density_factor_upper_bound':hi,
                     'subdivision_depth':maxdepth,'triangle_leaves':leaves,'unresolved_leaves':unresolved})
    result={'scope':'Continuous-density positivity bound for the specified L2 fit on these polygons, not a general moment closure or evolution proof.',
            'method':'Centroid triangle fan, adaptive subdivision; exact linear extrema at vertices plus monotone radial spline-potential range. 1e-11 floating margin; not arbitrary-precision interval arithmetic.',
            'candidate':str(candidate),'source':str(source),'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (candidate,source,Path(__file__))},
            'cells':len(rows),'all_cells_positive':all(r['unresolved_leaves']==0 for r in rows),
            'minimum_density_factor_lower_bound':min(r['density_factor_lower_bound'] for r in rows),
            'maximum_density_factor_upper_bound':max(r['density_factor_upper_bound'] for r in rows),
            'maximum_subdivision_depth':max(r['subdivision_depth'] for r in rows),
            'total_triangle_leaves':sum(r['triangle_leaves'] for r in rows),'cells_detail':rows}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cells_detail'},indent=2))


if __name__=='__main__':main()
