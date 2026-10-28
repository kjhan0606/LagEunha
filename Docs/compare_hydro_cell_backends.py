"""Backend CSR/source comparison with assembly and solve error kept separate."""
import hashlib,json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_hydro_cell_backends_20261018'
capture=base/'Codex_kinetic_face_capture_20261018'
ref=np.load(root/'reference.npz');expected=ref['expected'];raw=ref['raw'];pv=ref['pv']
S=ref['spatial_covariance'];N=len(raw);eps=np.finfo(float).eps
counts=ref['incidence_counts'];gamma=counts*eps/(1-counts*eps)
owners=np.fromfile(root/'owners.bin',dtype='<i8').reshape(-1,2)
frozen=np.fromfile(root/'frozen.bin',dtype='<i8').astype(bool)
a=np.fromfile(capture/'inputs.bin',dtype='<f8').reshape(-1,33)
# Independent absolute elementary-term bound before cancellation at nodes.
x=abs(a[:,:2]);vl=a[:,5:13].reshape(-1,4,2);vr=a[:,13:21].reshape(-1,4,2)
wl=abs(a[:,21:25]);wr=abs(a[:,25:29]);fm=wl.sum(1)+wr.sum(1)
fp=np.einsum('nk,nkj->nj',wl,abs(vl))+np.einsum('nk,nkj->nj',wr,abs(vr))+abs(a[:,31,None])*abs(a[:,2:4])
pi=np.einsum('nk,nki,nkj->nij',wl,abs(vl),abs(vl))+np.einsum('nk,nki,nkj->nij',wr,abs(vr),abs(vr))
energy=(wl*(.5*(vl*vl).sum(2)+abs(a[:,29,None]))).sum(1)+(wr*(.5*(vr*vr).sum(2)+abs(a[:,30,None]))).sum(1)+abs(a[:,31]*a[:,32])
nb=np.column_stack([fm,x*fm[:,None],(x[:,:,None]*x[:,None,:]*fm[:,None,None]).reshape(-1,4),fp,(x[:,:,None]*fp[:,None,:]).reshape(-1,4),pi.reshape(-1,4),energy])*a[:,4,None]
c=raw[:,1:3]/raw[:,0,None]
facebound=np.zeros((N,18));workbound=np.zeros(N);coordinate_H=np.zeros((N,2,2))
for (i,j),b,v in zip(owners,nb,a):
    for k in (i,j):
        if k>=0:
            facebound[k]+=b;workbound[k]+=abs(v[4]*v[31]*v[32])
            coordinate_H[k]+=np.outer(8*eps*(abs(v[:2])+abs(c[k])),abs(v[4]*v[31]*v[2:4]))
ftol=(64*eps+2*gamma[:,None])*facebound
Ftol=(32*eps+2*gamma[:,None])*ref['absolute_F']
Htol=(32*eps+2*gamma[:,None,None])*ref['absolute_H']+coordinate_H
c=raw[:,1:3]/raw[:,0,None];p=raw[:,7:9];D=raw[:,9:13].reshape(-1,2,2)
A=expected[:,96:100].reshape(-1,2,2);a0=expected[:,100:102]
qq=np.einsum('ni,nj->nij',raw[:,1:3],raw[:,1:3])/raw[:,0,None,None]
dS=8*eps*(abs(raw[:,3:7].reshape(-1,2,2))+abs(qq))
dH=Htol+8*eps*(abs(expected[:,74:78].reshape(-1,2,2))+abs(pv[:,None,None])*np.eye(2))
# First-order perturbation bound including covariance subtraction. Reject if
# perturbation is not small enough for the inverse resolvent bound (<1).
inv=abs(np.linalg.inv(S));perturb=inv@dS
eta=np.max(perturb.sum(2),axis=1)
assert eta.max()<1
AT=abs(A.transpose(0,2,1));H=abs(expected[:,74:78].reshape(-1,2,2)+pv[:,None,None]*np.eye(2))
snorm=np.max(abs(S).sum(2),axis=1)
# The solve error is a column norm bound distributed across both rows;
# componentwise relative error would miss mixing into small/zero entries.
solve_roundoff=64*eps*(snorm[:,None]*np.max(AT,axis=1)+np.max(H,axis=1))
g=inv@(dH+dS@AT+solve_roundoff[:,None,:])
deltaA=np.linalg.solve(np.eye(2)-perturb,g).transpose(0,2,1)
deltaa0=Ftol/raw[:,0,None]+np.einsum('nij,nj->ni',deltaA,abs(c))+16*eps*(abs(expected[:,72:74])/raw[:,0,None]+np.einsum('nij,nj->ni',abs(A),abs(c)))
pi_scale=np.einsum('ni,nj->nij',abs(a0),abs(p))+np.einsum('ni,nj->nij',abs(p),abs(a0))+abs(A)@abs(D)+abs(D.transpose(0,2,1))@abs(A.transpose(0,2,1))
pitol=np.einsum('ni,nj->nij',deltaa0,abs(p))+np.einsum('ni,nj->nij',abs(p),deltaa0)+deltaA@abs(D)+abs(D.transpose(0,2,1))@deltaA.transpose(0,2,1)+64*eps*pi_scale
stol=np.zeros((N,18));stol[:,1:3]=8*eps*abs(p)
stol[:,3:7]=8*eps*(abs(D)+abs(D.transpose(0,2,1))).reshape(-1,4)
stol[:,9:13]=8*eps*(abs(raw[:,13:17].reshape(-1,2,2))+abs(pv[:,None,None])*np.eye(2)).reshape(-1,4)
stol[:,13:17]=pitol.reshape(-1,4)
hydrotol=ftol+stol+8*eps*(abs(expected[:,:18])+facebound+abs(expected[:,78:96]))
ptol=np.zeros((N,18));ptol[:,7:9]=Ftol
ptol[:,9:13]=(Htol+abs(c[:,:,None])*Ftol[:,None,:]+16*eps*(abs(expected[:,74:78].reshape(-1,2,2))+abs(pv[:,None,None])*np.eye(2)+abs(c[:,:,None])*abs(expected[:,None,72:74]))).reshape(-1,4)
ptol[:,13:17]=pitol.reshape(-1,4);ptol[:,17]=(16*eps+2*gamma)*workbound
ptol[frozen]=0.;maskedtol=hydrotol.copy();maskedtol[frozen]=0.
def check(name,actual,target,tolerance):
    delta=abs(actual-target);bad=delta>tolerance
    if np.any(bad):
        where=np.unravel_index(np.argmax(np.where(bad,delta,0)),delta.shape)
        raise AssertionError((name,where,float(delta[where]),float(tolerance[where])))
    ratio=np.divide(delta,tolerance,out=np.zeros_like(delta),where=tolerance>0)
    return dict(maximum_absolute_difference=float(delta.max()),maximum_tolerance_fraction=float(ratio.max()))
report=dict(cells=N,nodes=len(a),maximum_S_condition=float(expected[:,103].max()),maximum_covariance_roundoff_perturbation=float(eta.max()),backends={})
arrays={}
for name in ('cpu1','omp4','cuda'):
    out=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(N,106)
    led=np.fromfile(root/(name+'_ledger.bin'),dtype='<f8')
    assert np.isfinite(out).all() and np.isfinite(led).all() and not np.any(out[:,105])
    r={}
    r['F']=check('F',out[:,72:74],expected[:,72:74],Ftol)
    r['Hface']=check('Hface',out[:,74:78].reshape(-1,2,2),expected[:,74:78].reshape(-1,2,2),Htol)
    r['source']=check('source',out[:,78:96],expected[:,78:96],stol)
    r['unmasked_hydro']=check('unmasked',out[:,:18],expected[:,:18],hydrotol)
    r['masked_hydro']=check('masked',out[:,18:36],expected[:,18:36],maskedtol)
    r['pressure']=check('pressure',out[:,36:54],expected[:,36:54],ptol)
    r['transport_remainder']=check('remainder',out[:,54:72],expected[:,54:72],maskedtol+ptol+8*eps*(abs(expected[:,18:36])+abs(expected[:,36:54])))
    assert np.max(out[:,102])<64*eps
    r['maximum_solve_backward_residual']=float(np.max(out[:,102]))
    r['A']=check('A',out[:,96:100].reshape(-1,2,2),A,deltaA)
    r['a0']=check('a0',out[:,100:102],a0,deltaa0)
    assert np.array_equal(out[:,78+4],out[:,78+5]) and np.array_equal(out[:,78+14],out[:,78+15])
    assert not np.any(out[frozen,18:72])
    reconstructed=out[:,36:54]+out[:,54:72]
    assert np.all(abs(reconstructed-out[:,18:36])<=8*eps*(abs(out[:,36:54])+abs(out[:,54:72])))
    # Ledger includes backend-ordered reduction; scale by primitive face and
    # cell terms, with reduction length explicit (no absolute epsilon floor).
    eb=nb[owners[:,1]<0][:,[7,8,17]].sum(0)
    angular=(nb[owners[:,1]<0,10]+nb[owners[:,1]<0,11]).sum()
    extbound=np.r_[eb,angular]
    frozenbound=np.r_[hydrotol[frozen,7:9].sum(0),hydrotol[frozen,17].sum(),(hydrotol[frozen,10]+hydrotol[frozen,11]).sum()]
    frozenabs=np.r_[abs(expected[frozen,7:9]).sum(0),abs(expected[frozen,17]).sum(),(abs(expected[frozen,10])+abs(expected[frozen,11])).sum()]
    ltol=np.r_[(128*eps+2*len(a)*eps)*extbound,frozenbound+2*N*eps*frozenabs]
    r['ledger']=check('ledger',led,ref['ledger'],ltol)
    r['sha256']=hashlib.sha256((root/(name+'.bin')).read_bytes()).hexdigest()
    report['backends'][name]=r;arrays[name]=out
report['cpu_omp_bitwise_equal']=np.array_equal(arrays['cpu1'].view(np.uint64),arrays['omp4'].view(np.uint64))
report['cpu_cuda_bitwise_equal']=np.array_equal(arrays['cpu1'].view(np.uint64),arrays['cuda'].view(np.uint64))
report['tolerance']='Node absolute elementary terms and ordered sum gamma_n; F/H absolute traction sums plus centroid/subtraction roundoff; source S/H perturbation bound solve(I-abs(invS)@dS,g) with column-norm solve roundoff. Solve backward residual <64eps checked separately.'
report['scope']='Actual backend node flux, CSR cell assembly, compatible pressure/volume source and exterior/frozen ledger on fixed captured state; excludes factor K/Qfactor, gravity, massroot, geometry/reconstruction, MPI, evolution and Psi.'
report['provenance']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('moment_hydro_cell_source.h','moment_hydro_cell_harness.cu','hydro_harness')}
(root/'comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
