"""Fixed coefficient factor maps: common-input accuracy and stage propagation."""
import hashlib,json,math
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_factor_update_backends_v2_20261018'
manifest=json.loads((root/'manifest.json').read_text());N=manifest['cells'];eps=np.finfo(float).eps;rt=np.sqrt(2.)
offset=np.fromfile(root/'offsets.bin',dtype='<i8');report=dict(scope=__doc__,heights=manifest['heights'],comparisons={})
assert manifest['inputs']['cpu1']==manifest['inputs']['omp4']
def check(label,actual,target,tolerance):
    delta=abs(actual-target);bad=delta>tolerance
    if np.any(bad):
        where=np.unravel_index(np.argmax(np.where(bad,delta,0)),delta.shape)
        raise AssertionError((label,where,float(delta[where]),float(np.broadcast_to(tolerance,delta.shape)[where])))
    ratio=np.divide(delta,tolerance,out=np.zeros_like(delta),where=tolerance>0)
    return dict(maximum_absolute_difference=float(delta.max()),maximum_tolerance_fraction=float(ratio.max()))
def symmetric_operator(R):
    return np.array([[2*R[0,0],rt*R[0,1],0.],[rt*R[0,1],R[0,0]+R[1,1],rt*R[0,1]],[0.,rt*R[0,1],2*R[1,1]]])
def make_bounds(deck,expected,Q):
    tol=np.zeros_like(expected);metrics=[]
    for index,(a,target) in enumerate(zip(deck,expected)):
        if a[35]:metrics.append(None);continue
        i=index%N;baseL=a[:4].reshape(2,2);linear=a[4:8].reshape(2,2);A=a[8:12].reshape(2,2);K=a[12:16].reshape(2,2);raw=a[16:34];h=a[34]
        d=target[22:62];loss=d[:4].reshape(2,2);Ac=d[4:8].reshape(2,2);eh=d[8:12].reshape(2,2);ef=d[12:16].reshape(2,2)
        R=linear@linear.T;J=linear@K@linear.T;M=symmetric_operator(R)
        dR=32*eps*(abs(linear)@abs(linear.T));dJ=128*eps*(abs(linear)@abs(K)@abs(linear.T))
        dM=symmetric_operator(dR);invm=abs(np.linalg.inv(M));perturb=invm@dM;eta=np.linalg.norm(perturb,np.inf)
        assert eta<1
        z=np.array([loss[0,0],rt*loss[0,1],loss[1,1]]);rhs=np.array([J[0,0],rt*J[0,1],J[1,1]])
        drhs=np.array([dJ[0,0],rt*dJ[0,1],dJ[1,1]])
        roundoff=256*eps*(np.linalg.norm(M,np.inf)*np.max(abs(z))+np.max(abs(rhs)))
        dz=np.linalg.solve(np.eye(3)-perturb,invm@(drhs+dM@abs(z)+roundoff))
        dl=np.array([[dz[0],dz[1]/rt],[dz[1]/rt,dz[2]]])+32*eps*np.linalg.norm(loss,np.inf)
        dac=dl+16*eps*(abs(A)+abs(loss));dx=.5*h*np.linalg.norm(dac,np.inf)
        Xnorm=.5*h*np.linalg.norm(Ac,np.inf);scaled=Xnorm;s=0
        while scaled>.5:scaled*=.5;s+=1
        trunc=math.exp(scaled)*scaled**19/math.factorial(19)
        exponential_error=trunc+256*eps*math.exp(scaled)
        matrix_norm=math.exp(scaled)
        for _ in range(s):
            exponential_error=2*matrix_norm*exponential_error+exponential_error**2+32*eps*matrix_norm**2
            matrix_norm*=matrix_norm
        deh=exponential_error+math.exp(Xnorm+dx)*dx+128*eps*math.exp(Xnorm)
        defull=2*np.linalg.norm(eh,np.inf)*deh+deh**2+32*eps*np.linalg.norm(eh,np.inf)**2+128*eps*math.exp(2*Xnorm)
        qi=Q[offset[i]:offset[i+1]].T
        U=np.column_stack([ef@baseL,np.sqrt(h)*(eh@qi)])
        du0=defull*np.ones((2,2))@abs(baseL)+32*eps*(abs(ef)@abs(baseL))
        duq=np.sqrt(h)*(deh*np.ones((2,2))@abs(qi)+32*eps*(abs(eh)@abs(qi)))
        du=np.column_stack([du0,duq]);cols=U.shape[1]
        gamma=64*cols*eps/(1-64*cols*eps);condition=float(np.linalg.cond(U))
        dL=16*condition*(np.linalg.norm(du)+gamma*np.linalg.norm(U))
        newL=target[:4].reshape(2,2);newR=d[16:20].reshape(2,2)
        dr=2*dL*(np.ones((2,2))@abs(newL.T)+abs(newL)@np.ones((2,2)))+2*dL*dL+32*eps*(abs(newL)@abs(newL.T))
        m=raw[0];q=raw[1:3];p=raw[7:9];S=raw[3:7].reshape(2,2)-np.outer(q,q)/m;C=raw[9:13].reshape(2,2)-np.outer(q,p)/m
        B=d[20:24].reshape(2,2)
        ds=16*eps*(abs(raw[3:7].reshape(2,2))+abs(np.outer(q,q)/m));dc=16*eps*(abs(raw[9:13].reshape(2,2))+abs(np.outer(q,p)/m))
        invs=abs(np.linalg.inv(S));ps=invs@ds;etas=np.linalg.norm(ps,np.inf);assert etas<1
        solve_roundoff=128*eps*(np.linalg.norm(S,np.inf)*np.max(abs(B.T),axis=0)+np.max(abs(C),axis=0))
        db=np.linalg.solve(np.eye(2)-ps,invs@(dc+ds@abs(B.T)+solve_roundoff[None,:])).T
        mix=db@abs(S)@abs(B.T)+abs(B)@abs(S)@db.T+abs(B)@ds@abs(B.T)
        mix+=db@ds@abs(B.T)+abs(B)@ds@db.T+db@abs(S)@db.T+db@ds@db.T
        pi_scale=abs(np.outer(p,p)/m)+abs(B)@abs(S)@abs(B.T)+abs(newR)
        dpi=mix+dr+128*eps*pi_scale
        thermal_tol=.5*np.trace(dpi)+16*eps*(abs(raw[17])+.5*(abs(target[17])+abs(target[20])))
        tol[index,:4]=dL;tol[index,4+13:4+17]=dpi.ravel()
        td=tol[index,22:62];td[:4]=dl.ravel();td[4:8]=dac.ravel();td[8:12]=deh;td[12:16]=defull;td[16:20]=dr.ravel();td[20:24]=db.ravel()
        td[24:28]=dpi.ravel()+16*eps*(abs(target[17:21])+abs(raw[13:17]));td[28:31]=thermal_tol
        metrics.append(dict(condition_U=condition,qr_gamma=gamma,exp_bound=deh,exp_norm=Xnorm,maximum_input_factor_perturbation=float(np.linalg.norm(du)),qr_factor_bound=float(dL),U=U,condition_S=float(np.linalg.cond(S,np.inf))))
    return tol,metrics
refs={};decks={};bounds={};factors={}
for deck_name in ('cpu1','cuda'):
    decks[deck_name]=np.fromfile(root/(deck_name+'_input.bin'),dtype='<f8').reshape(-1,36)
    refs[deck_name]=np.load(root/(deck_name+'_reference.npz'))['expected'].reshape(-1,63)
    Q=np.fromfile(root/(deck_name+'_Q.bin'),dtype='<f8').reshape(-1,2)
    bounds[deck_name]=make_bounds(decks[deck_name],refs[deck_name],Q)
for name,deck_name in (('common_cpu1','cpu1'),('common_omp4','cpu1'),('common_cuda','cpu1'),('linked_cuda','cuda')):
    out=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(-1,63);target=refs[deck_name];deck=decks[deck_name];tol,metrics=bounds[deck_name]
    assert np.isfinite(out).all() and not np.any(out[:,62])
    frozen=deck[:,35].astype(bool);active=~frozen;r={}
    for label,sl in [('Lnew',slice(0,4)),('rebuilt_raw',slice(4,22)),('loss',slice(22,26)),('Ac',slice(26,30)),('Ehalf',slice(30,34)),('Efull',slice(34,38)),('Rnew',slice(38,42)),('Bprovisional',slice(42,46)),('deltaPi',slice(46,50)),('thermal_partition',slice(50,53))]:
        r[label]=check(label,out[:,sl],target[:,sl],tol[:,sl])
    assert np.array_equal(out[frozen,:4].view(np.uint64),deck[frozen,:4].view(np.uint64))
    assert np.array_equal(out[frozen,4:22].view(np.uint64),deck[frozen,16:34].view(np.uint64))
    assert not np.any(out[frozen,22:62])
    for j in list(range(13))+[17]:assert np.array_equal(out[:,4+j].view(np.uint64),deck[:,16+j].view(np.uint64))
    assert np.array_equal(out[:,18],out[:,19])
    assert np.all(out[active,0]>0) and np.all(out[active,3]>0) and not np.any(out[active,1])
    assert np.all(out[active,51]>0) and np.max(out[:,53])<128*eps
    for index in np.flatnonzero(active):
        metric=metrics[index];U=metric['U'];Lnew=out[index,:4].reshape(2,2)
        thin=np.linalg.solve(Lnew,U).T
        orth=np.linalg.norm(thin.T@thin-np.eye(2),np.inf)
        assert orth<64*metric['qr_gamma']*metric['condition_U']**2+64*metric['exp_bound']
        assert out[index,58]<64*metric['qr_gamma']*metric['condition_U']**2
        assert out[index,59]<128*eps*np.max(abs(U))
        assert out[index,60]<64*metric['qr_gamma']*out[index,61]
        check('independent UUT',Lnew@Lnew.T,U@U.T,2*tol[index,38:42].reshape(2,2))
        assert out[index,56]<=.5 and out[index,55]<=60
    r['maximum_loss_backward_residual']=float(np.max(out[:,53]));r['maximum_QR_orthogonality_residual']=float(np.max(out[:,58]));r['maximum_QR_gram_residual']=float(np.max(out[:,60]));r['minimum_final_thermal']=float(np.min(out[active,51]))
    r['input_deck']=deck_name;r['sha256']=hashlib.sha256((root/(name+'.bin')).read_bytes()).hexdigest();report['comparisons'][name]=r;factors[name]=out
report['common_cpu_omp_bitwise_equal']=np.array_equal(factors['common_cpu1'].view(np.uint64),factors['common_omp4'].view(np.uint64))
report['linked_cpu_omp_reused_common']='Input/Q/reference SHAs identical; common CPU/OpenMP runs are also their native-stage runs. No duplicate calculation.'
effect=[]
for k,h in enumerate(manifest['heights']):
    sl=slice(k*N,(k+1)*N);cpu=factors['common_cpu1'][sl];gpu=factors['linked_cuda'][sl];same=factors['common_cuda'][sl]
    stage_cpu=decks['cpu1'][sl];stage_gpu=decks['cuda'][sl]
    effect.append(dict(h=h,h_times_maximum_covA_input_difference=float(h*abs(stage_gpu[:,8:12]-stage_cpu[:,8:12]).max()),common_input_factor_maximum=float(abs(same[:,:4]-cpu[:,:4]).max()),linked_factor_maximum=float(abs(gpu[:,:4]-cpu[:,:4]).max()),linked_R_maximum=float(abs(gpu[:,38:42]-cpu[:,38:42]).max()),linked_Pi_maximum=float(abs(gpu[:,17:21]-cpu[:,17:21]).max()),linked_thermal_maximum=float(abs(gpu[:,51]-cpu[:,51]).max()),python_predicted_linked_R_maximum=float(abs(refs['cuda'][sl,38:42]-refs['cpu1'][sl,38:42]).max())))
report['stage_error_propagation']=effect
report['tolerance']='Separate S/R Sylvester inverse-resolvent input roundoff, normwise exponential perturbation and Taylor/squaring bound, Givens QR normwise backward error with column count and condition, provisional S/C solve and Pi/thermal propagation. No absolute epsilon floor; finite normal-range audit, not universal interval certification.'
report['limits']='Two frozen-coefficient maps h5e-5/1e-4 only. Not actual midpoint evolution or temporal-order proof; rank-deficient Llinear unsupported; moving geometry, massroot, MPI and long Kepler accuracy remain unresolved.'
report['provenance']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('moment_factor_update.h','moment_factor_update_harness.cu','factor_harness')}
(root/'comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
