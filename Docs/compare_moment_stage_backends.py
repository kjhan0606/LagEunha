"""Fixed captured stage: actual backend transport, gravity and covariance."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_moment_stage_capture_20261018'
sys.path.insert(0,str(root))
import hydro_cell_error_bounds as hb
ref=np.load(root/'reference.npz');inputs=np.load(root/'cell_inputs.npz')
expected=ref['expected'];raw=inputs['raw'];L=inputs['factors'];frozen=inputs['frozen'];N=len(raw)
eps=np.finfo(float).eps;check=hb.check
prep=base/'Codex_hydro_cell_backends_20261018'
off=np.fromfile(prep/'offsets.bin',dtype='<i8');inc=np.fromfile(prep/'incidence.bin',dtype='<i8')
goff=np.fromfile(root/'gravity_offsets.bin',dtype='<i8')
gn=np.fromfile(root/'gravity_nodes.bin',dtype='<f8').reshape(-1,7)
dm=gn[:,2];x=abs(gn[:,:2]);v=abs(gn[:,3:5]);acc=abs(gn[:,5:7])
primitive=np.zeros((len(gn),18));primitive[:,7:9]=dm[:,None]*acc
primitive[:,9:13]=(dm[:,None,None]*x[:,:,None]*acc[:,None,:]).reshape(-1,4)
primitive[:,13:17]=(dm[:,None,None]*(v[:,:,None]*acc[:,None,:]+acc[:,:,None]*v[:,None,:])).reshape(-1,4)
primitive[:,17]=dm*(v*acc).sum(1)
gabs=np.add.reduceat(primitive,goff[:-1],axis=0)
gc=np.diff(goff);gamma=gc*eps/(1-gc*eps)
gtol=(64*eps+2*gamma[:,None])*gabs;grav=expected[:,30:48]
q=raw[:,1:3];p=raw[:,7:9];m=raw[:,0];B=expected[:,18:22].reshape(-1,2,2)
S=hb.S;invS=abs(np.linalg.inv(S));dS=hb.dS
C=raw[:,9:13].reshape(-1,2,2)-q[:,:,None]*p[:,None,:]/m[:,None,None]
dCstate=8*eps*(abs(raw[:,9:13].reshape(-1,2,2))+abs(q[:,:,None]*p[:,None,:]/m[:,None,None]))
BT=abs(B.transpose(0,2,1));snorm=np.max(abs(S).sum(2),axis=1)
bnorm=np.max(BT,axis=1);cnorm=np.max(abs(C),axis=1)
roundoff=128*eps*(snorm[:,None]*bnorm+cnorm)
perturb=invS@dS;etaS=np.max(perturb.sum(2),axis=1);assert etaS.max()<1
dB=np.linalg.solve(np.eye(2)-perturb,invS@(dCstate+dS@BT+roundoff[:,None,:])).transpose(0,2,1)
gC=grav[:,9:13].reshape(-1,2,2)-q[:,:,None]*grav[:,None,7:9]/m[:,None,None]
gT=grav[:,13:17].reshape(-1,2,2)-(grav[:,7:9,None]*p[:,None,:]+p[:,:,None]*grav[:,None,7:9])/m[:,None,None]
dgC=gtol[:,9:13].reshape(-1,2,2)+abs(q[:,:,None])*gtol[:,None,7:9]/m[:,None,None]+16*eps*(abs(grav[:,9:13].reshape(-1,2,2))+abs(q[:,:,None]*grav[:,None,7:9]/m[:,None,None]))
dgT=gtol[:,13:17].reshape(-1,2,2)+(gtol[:,7:9,None]*abs(p[:,None,:])+abs(p[:,:,None])*gtol[:,None,7:9])/m[:,None,None]+16*eps*(abs(grav[:,13:17].reshape(-1,2,2))+(abs(grav[:,7:9,None]*p[:,None,:])+abs(p[:,:,None]*grav[:,None,7:9]))/m[:,None,None])
mix=dB@abs(gC)+abs(B)@dgC+dB@dgC
dG=dgT+mix+mix.transpose(0,2,1)+64*eps*(abs(gT)+abs(B)@abs(gC)+abs(gC.transpose(0,2,1))@abs(B.transpose(0,2,1)))
dG=.5*(dG+dG.transpose(0,2,1))
R=L@L.transpose(0,2,1);dR=16*eps*(abs(L)@abs(L.transpose(0,2,1)));rt=np.sqrt(2.)
def operator_matrix(r):
    t=np.zeros((N,3,3));t[:,0,0]=2*r[:,0,0];t[:,1,1]=r[:,0,0]+r[:,1,1];t[:,2,2]=2*r[:,1,1]
    t[:,0,1]=t[:,1,0]=rt*r[:,0,1];t[:,1,2]=t[:,2,1]=rt*r[:,0,1]
    return t
T=operator_matrix(R);dT=operator_matrix(dR);absinv=abs(np.linalg.inv(T))
Ag=expected[:,52:56].reshape(-1,2,2);G=expected[:,48:52].reshape(-1,2,2)
z=np.column_stack([Ag[:,0,0],rt*Ag[:,0,1],Ag[:,1,1]])
rhs=np.column_stack([G[:,0,0],rt*G[:,0,1],G[:,1,1]])
drhs=np.column_stack([dG[:,0,0],rt*dG[:,0,1],dG[:,1,1]])
perturbR=absinv@dT;etaR=np.max(perturbR.sum(2),axis=1);assert etaR.max()<1
t_norm=np.max(abs(T).sum(2),axis=1)
solve_roundoff=128*eps*(t_norm*np.max(abs(z),axis=1)+np.max(abs(rhs),axis=1))
g=np.einsum('nij,nj->ni',absinv,drhs+np.einsum('nij,nj->ni',dT,abs(z))+solve_roundoff[:,None])
dz=np.linalg.solve(np.eye(3)-perturbR,g[:,:,None])[:,:,0]
dAg=np.stack([dz[:,0],dz[:,1]/rt,dz[:,1]/rt,dz[:,2]],axis=1).reshape(-1,2,2)
Qref=ref['Q'];Qun=Qref[:,:8].reshape(-1,4,2);qtol=32*eps*ref['Qabsolute']
qcombinedtol=np.concatenate([qtol.reshape(-1,8),qtol.reshape(-1,8)],axis=1)
incomingtol=np.zeros((N,2,2));kab=ref['Kabsolute'];count=4*np.diff(off)
kgamma=count*eps/(1-count*eps);ktol=(64*eps+2*kgamma[:,None,None])*kab
for i in range(N):
    sl=slice(off[i],off[i+1]);qq=abs(Qun[sl]);dq=qtol[sl]
    value=np.einsum('nsi,nsj->ij',dq,qq)+np.einsum('nsi,nsj->ij',qq,dq)+np.einsum('nsi,nsj->ij',dq,dq)
    value+=(32*eps+2*kgamma[i])*np.einsum('nsi,nsj->ij',qq,qq)
    incomingtol[i]=value
    if frozen[i]:qcombinedtol[sl,8:]=0.
total_tol=hb.maskedtol+gtol+8*eps*(abs(expected[:,:18])+abs(grav));total_tol[frozen]=0.
unmasktol=hb.hydrotol+gtol+8*eps*(abs(expected[:,75:93])+abs(grav))
covtol=dB+dAg+8*eps*(abs(B)+abs(Ag));maskedcovtol=covtol.copy();maskedcovtol[frozen]=0.
report=dict(cells=N,nodes=len(hb.a),gravity_nodes=len(gn),maximum_S_condition=float(np.max(np.linalg.cond(S))),maximum_R_condition=float(np.max(np.linalg.cond(R))),maximum_sylvester_condition=float(np.max(np.linalg.cond(T,np.inf))),maximum_eta_S=float(etaS.max()),maximum_eta_R=float(etaR.max()),backends={})
arrays={}
for name in ('cpu1','omp4','cuda'):
    path=root/(name+'.bin');out=np.fromfile(path,dtype='<f8').reshape(N,94)
    Q=np.fromfile(root/(name+'_Q.bin'),dtype='<f8').reshape(-1,16)
    led=np.fromfile(root/(name+'_gravity_ledger.bin'),dtype='<f8')
    assert np.isfinite(out).all() and np.isfinite(Q).all() and not np.any(out[:,93])
    r={}
    for label,sl,tol in [('raw_total',slice(0,18),total_tol),('raw_unmasked',slice(75,93),unmasktol),('gravity',slice(30,48),gtol),('B',slice(18,22),dB.reshape(-1,4)),('G',slice(48,52),dG.reshape(-1,4)),('Ag',slice(52,56),dAg.reshape(-1,4)),('cov_A',slice(56,60),maskedcovtol.reshape(-1,4)),('K_unmasked',slice(63,67),ktol.reshape(-1,4)),('incoming_unmasked',slice(67,71),incomingtol.reshape(-1,4))]:
        r[label]=check(label,out[:,sl],expected[:,sl],tol)
    r['Qfactor']=check('Qfactor',Q,Qref,qcombinedtol)
    maskedktol=ktol.copy();maskedktol[frozen]=0.
    maskedintol=incomingtol.copy();maskedintol[frozen]=0.
    r['K_masked']=check('K_masked',out[:,22:26],expected[:,22:26],maskedktol.reshape(-1,4))
    r['incoming_masked']=check('incoming_masked',out[:,26:30],expected[:,26:30],maskedintol.reshape(-1,4))
    assert not np.any(out[frozen,:18]) and not np.any(out[frozen,22:30]) and not np.any(out[frozen,56:60])
    for i in range(N):
        qq=Q[off[i]:off[i+1],:8].reshape(-1,2)
        check('QQT',qq.T@qq,out[i,67:71].reshape(2,2),(incomingtol[i]+64*eps*(abs(qq).T@abs(qq))))
    assert np.max(out[:,61])<128*eps
    r['maximum_sylvester_backward_residual']=float(np.max(out[:,61]))
    r['maximum_pre_symmetry_G_difference']=float(np.max(abs(out[:,60])))
    thermal=out[:,47]-.5*(out[:,43]+out[:,46]);torque=out[:,40]-out[:,41]
    r['gravity_thermal_identity']=check('gravity thermal',thermal,np.zeros(N),gtol[:,17]+.5*(gtol[:,13]+gtol[:,16])+8*eps*gabs[:,17])
    r['gravity_origin_torque']=check('gravity torque',torque,np.zeros(N),gtol[:,10]+gtol[:,11])
    ltol=np.r_[gtol[~frozen,7:9].sum(0),gtol[~frozen,17].sum(),(gtol[~frozen,10]+gtol[~frozen,11]).sum()]
    labs=np.r_[gabs[~frozen,7:9].sum(0),gabs[~frozen,17].sum(),(gabs[~frozen,10]+gabs[~frozen,11]).sum()]
    r['gravity_ledger']=check('gravity ledger',led,np.load(root/'gravity_ledger.npy'),ltol+2*N*eps*labs)
    # Independent stage identity catches a common sign error in both compiled
    # gravity ledger and its Python reference. Hydro ledger is already audited.
    hledger=np.fromfile(prep/(name+'_ledger.bin'),dtype='<f8').reshape(2,4).sum(0)
    total=np.r_[out[:,7:9].sum(0),out[:,17].sum(),(out[:,10]-out[:,11]).sum()]
    absolute=np.r_[abs(out[:,7:9]).sum(0),abs(out[:,17]).sum(),(abs(out[:,10])+abs(out[:,11])).sum()]
    primitive_bound=np.r_[(hb.facebound[:,7:9]+gabs[:,7:9]).sum(0),(hb.facebound[:,17]+gabs[:,17]).sum(),(hb.facebound[:,10]+hb.facebound[:,11]+gabs[:,10]+gabs[:,11]).sum()]
    r['stage_global_ledger_identity']=check('global stage ledger',total,hledger+led,(128*eps+2*N*eps)*(absolute+primitive_bound))
    r['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();report['backends'][name]=r;arrays[name]=out
report['cpu_omp_bitwise_equal']=np.array_equal(arrays['cpu1'].view(np.uint64),arrays['omp4'].view(np.uint64))
report['cpu_cuda_bitwise_equal']=np.array_equal(arrays['cpu1'].view(np.uint64),arrays['cuda'].view(np.uint64))
report['tolerance']='Primitive gravity/transport absolute sums and gamma_n reductions; component-mixing S/C and symmetric Sylvester perturbation inverse resolvents with normwise solve roundoff; no absolute epsilon floor. Residuals checked separately.'
report['scope']='Fixed captured reconstruction, face motion and q12 gravity quadrature stage RHS; no factor time update, geometry/reconstruction/massroot implementation, MPI, evolving stage or long-time accuracy claim.'
report['provenance']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('moment_covariance_stage.h','moment_stage_harness.cu','stage_harness')}
(root/'comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
