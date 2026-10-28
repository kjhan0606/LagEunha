"""Existing covariance-budget reconstruction: smooth-field face convergence.
No evolution, gravity, or production-backend claim.
"""
import argparse,json,os,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_coupled_relaxation_evolution_20261018')
from adaptive_moment_batch import setup
import moment_ale_evolve as model
import quadratic_velocity_reconstruction as quad
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--samples',action='store_true')
parser.add_argument('--riemann',action='store_true',help='Also check shared HLLC of resolved means at stationary faces')
args=parser.parse_args()
alpha=np.sqrt(2)*np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.]])

def exact(p):
    p=np.atleast_2d(p)
    return np.column_stack([.2+.1*np.sin(2*np.pi*p[:,1]),.1*np.cos(2*np.pi*p[:,0])])

results=[]
for n in (4,8,16):
    mass,y=setup('shear',n,0.);a=y[:-4].reshape(-1,20)
    polys,edges=model.geometry(a[:,:2],a[:,2],0.)
    states=[];factors=[];factor_discrepancy=0.
    for i,poly in enumerate(polys):
        raw=np.r_[mass[i],a[i,3:]].copy();s=model.decode(poly,raw)
        pts=s['density'].points;dm=s['density'].weighted_mass;v=exact(pts)
        raw[7:9]=dm@v;raw[9:13]=(pts.T@(dm[:,None]*v)).ravel()
        Pi=v.T@(dm[:,None]*v);raw[13:17]=Pi.ravel()
        raw[17]=.002777777777777778*s['density'].area/(model.GAMMA-1)+.5*np.trace(Pi)
        s=model.decode(poly,raw)
        residual=v-s['v']-(pts-s['c'])@s['B'].T
        _,R=np.linalg.qr(np.sqrt(dm)[:,None]*residual,mode='reduced');L=R.T
        factor_discrepancy=max(factor_discrepancy,float(np.max(abs(L@L.T-s['R']))))
        states.append(s);factors.append(L)
    for mode in ('affine','quadratic'):
        for s in states:s.pop('velocity_quadratic',None)
        diag={}
        if mode=='quadratic':diag=quad.prepare(states,polys,edges,mass,factors)
        sums=np.zeros(4);lengthsum=0.;contacts=0;sample_identity_error=0.;hsums=np.zeros(3);hlength=0.
        gx,gw=np.polynomial.legendre.leggauss(4)
        for i,j,p0,p1,normal,length in edges:
            if j<0:continue
            for z,w in zip(gx,gw):
                point=p0+(z+1)/2*(p1-p0);truth=exact(point)[0];weight=length*w/2
                if args.riemann:
                    pressure=.002777777777777778
                    ul=model.state(1.,model.velocity(states[i],point),pressure)
                    ur=model.state(1.,model.velocity(states[j],point),pressure)
                    flux,_,_=model.sample(model.fan(ul,ur,normal),0.)
                    un=truth@normal
                    reference=np.r_[un,truth*un+pressure*normal,(.5*(truth@truth)+model.GAMMA/(model.GAMMA-1)*pressure)*un]
                    delta=flux-reference
                    hsums+=weight*np.array([delta[0]**2,sum(delta[1:3]**2),delta[3]**2]);hlength+=weight
                for owner in (i,j):
                    v=model.velocity(states[owner],point);vn=v@normal;un=truth@normal
                    momentum=np.outer(v,v)@normal;reference=np.outer(truth,truth)@normal
                    energy=(.5*(v@v)+model.GAMMA/(model.GAMMA-1)*.002777777777777778)*vn
                    if args.samples:
                        velocities,_=quad.velocity_samples(states[owner],factors[owner],mass[owner],point,alpha)
                        samples=velocities[0];normal_v=samples@normal
                        fluct=samples-v;cov=fluct.T@fluct/len(samples)
                        sample_momentum=np.mean(samples*normal_v[:,None],axis=0)
                        sample_energy=float(np.mean((.5*np.sum(samples*samples,axis=1)+model.GAMMA/(model.GAMMA-1)*.002777777777777778)*normal_v))
                        expected_p=momentum+cov@normal
                        expected_e=energy+.5*np.trace(cov)*vn+v@cov@normal
                        sample_identity_error=max(sample_identity_error,float(np.max(abs(samples.mean(0)-v))),float(np.max(abs(sample_momentum-expected_p))),abs(sample_energy-expected_e))
                        momentum=sample_momentum;energy=sample_energy
                    exactE=(.5*(truth@truth)+model.GAMMA/(model.GAMMA-1)*.002777777777777778)*un
                    sums+=weight*np.array([sum((v-truth)**2),(vn-un)**2,sum((momentum-reference)**2),(energy-exactE)**2])
                    lengthsum+=weight;contacts+=1
        errors=np.sqrt(sums/lengthsum)
        results.append(dict(n=n,mode=mode,velocity_RMS=errors[0],mass_flux_RMS=errors[1],
            momentum_flux_RMS=errors[2],energy_flux_RMS=errors[3],contacts=contacts,
            covariance_factor_discrepancy=factor_discrepancy,sample_identity_error=sample_identity_error,
            shared_mean_HLLC_RMS=(np.sqrt(hsums/hlength).tolist() if args.riemann else None),diagnostic=diag))
        assert sample_identity_error<1e-12
ratios={}
for mode in ('affine','quadratic'):
    rows=[r for r in results if r['mode']==mode]
    ratios[mode]={key:[rows[i][key]/rows[i+1][key] for i in (0,1)]
                  for key in ('velocity_RMS','mass_flux_RMS','momentum_flux_RMS','energy_flux_RMS')}
print(json.dumps(dict(scope=__doc__,closure=quad.closure_name(),residual_samples=args.samples,shared_mean_hllc=args.riemann,results=results,refinement_ratios=ratios,
    limits='Manufactured smooth sinusoidal velocity, rho1, low uniform pressure. Optional residual-sample traces and separate mean-state stationary HLLC diagnostic. No residual-sample Riemann coupling, weight-rate solve, or evolution. Different seeded irregular meshes are not nested.'),indent=2))
