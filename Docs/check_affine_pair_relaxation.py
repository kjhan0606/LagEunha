"""Finite affine pair impulse check on saved remap states; synthetic contacts."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from affine_pair_relaxation import relax,kinetic,angular
path=Path('/gpfs/kjhan/LagForce/Codex_checkpoint_weight_remap_20261018/state_w1_q8.npz')
with np.load(path) as data:q=data['state']
def unpack(q):
    m=q[0];c=q[1:3]/m;v=q[6:8]/m
    S=np.array([[q[3],q[4]],[q[4],q[5]]])-m*np.outer(c,c)
    C=q[8:12].reshape(2,2)-m*np.outer(c,v)
    T=np.array([[q[12],q[13]],[q[13],q[14]]])-m*np.outer(v,v)
    G=np.linalg.solve(S,C).T
    return dict(m=m,c=c,v=v,S=S,G=G,R=T-G@S@G.T,U=q[15]-.5*(q[12]+q[14]))
states=[unpack(x) for x in q];centers=np.array([s['c'] for s in states])
errors=dict(energy=0.,angular=0.,momentum=0.,slip=0.,semigroup=0.,null_impulse=0.)
heat_total=0.;count=0
for i in range(64):
    dist=np.linalg.norm(centers-centers[i],axis=1);dist[i]=np.inf;j=int(np.argmin(dist))
    a,b=states[i],states[j];delta=b['c']-a['c'];normal=delta/np.linalg.norm(delta)
    tangent=np.array([-normal[1],normal[0]])
    point=.5*(a['c']+b['c'])+.2*np.linalg.norm(delta)*tangent
    E0=sum(kinetic(s)+s['U'] for s in [a,b]);L0=angular(a)+angular(b);p0=a['m']*a['v']+b['m']*b['v']
    for dt in [0.,1e-5,.01,1.,1e6]:
        aa,bb,d=relax(a,b,point,tangent,.5,dt);count+=1
        errors['energy']=max(errors['energy'],abs(sum(kinetic(s)+s['U'] for s in [aa,bb])-E0))
        errors['angular']=max(errors['angular'],abs(angular(aa)+angular(bb)-L0))
        errors['momentum']=max(errors['momentum'],float(np.max(abs(aa['m']*aa['v']+bb['m']*bb['v']-p0))))
        slip=((aa['v']+aa['G']@(point-aa['c']))-(bb['v']+bb['G']@(point-bb['c'])))@tangent
        errors['slip']=max(errors['slip'],abs(slip-d['slip_after_expected']))
        assert d['generated_heat']>=0 and aa['U']>=a['U'] and bb['U']>=b['U']
        assert np.array_equal(aa['R'],a['R']) and np.array_equal(bb['R'],b['R'])
        ha,hb,_=relax(a,b,point,tangent,.5,dt/2)
        ha,hb,_=relax(ha,hb,point,tangent,.5,dt/2)
        for end,half in [(aa,ha),(bb,hb)]:
            for key in ['v','G','U']:errors['semigroup']=max(errors['semigroup'],float(np.max(abs(end[key]-half[key]))))
        heat_total+=d['generated_heat']
    # Full affine shear null, using actual mass-center/inertia geometry.
    G=np.array([[0.,1.1],[0.,0.]])
    na=dict(a,G=G,v=G@a['c']+np.array([.2,-.1]))
    nb=dict(b,G=G,v=G@b['c']+np.array([.2,-.1]))
    _,_,d=relax(na,nb,point,tangent,.5,1.)
    errors['null_impulse']=max(errors['null_impulse'],float(np.linalg.norm(d['impulse'])))
assert max(errors.values())<1e-10,errors
out=dict(scope=__doc__,pairs=64,cases=count,max_errors=errors,
         aggregate_heat_over_independent_tests=heat_total,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
         limits=['Geometry frozen within each pair update; synthetic contact points, not actual shared faces.',
                 'No production timestep, force, spin or heat state changed.',
                 'Residual covariance retained; its transport closure remains unresolved.',
                 'Sequential face updates need ordering/time-convergence and curved-flow checks.'])
print(json.dumps(out,indent=2))
