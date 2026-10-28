"""Instantaneous departure-kinetic-energy generator on matched t6 states."""
import argparse, json, struct, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
from audit_canonical_pair_history import gravity_coefficient
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_t6_entropy_force_probes_20261018'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=root)
root=parser.parse_args().root
pressure=json.loads((root/'entropy_pressure_work.json').read_text())
ic,_=read(b/'Codex_canonical_initial_20261017/run/khout.000000.dat')
r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'inner_background':(r0<.5)&(ic['den']<.5),'inner_shear_band':(r0>=.3)&(r0<.5)}
out={'scope':'Instantaneous generator of half sum(m*|v-vcirc|^2). HG is its zero-step modal generator, not an integrated split-map budget. Hydro remainder includes TCF and any other hydro RHS terms. No eigenmode or long-time stability claim.', 'cases':{}}
for name in ('adaptive','legacy'):
    d=root/name;a,t=read(d/'khout.001783.dat');n=len(a)
    assert sha(d/'khout.001783.dat')==pressure['cases'][name]['source_sha256']
    state=np.loadtxt(d/'raw_faces.dat.hgstate.rank0000');state=state[np.argsort(state[:,0])]
    rhs=np.loadtxt(d/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
    assert state.shape==(n,13) and rhs.shape==(n,25)
    assert np.array_equal(state[:,0],np.arange(n))
    for key,col in [('x',3),('y',4),('mass',5),('vx',6),('vy',7),('ie',8)]:
        assert np.max(abs(state[:,col]-a[key]))<1e-12
    path=d/'raw_faces.dat.hgm.rank0000'
    with path.open('rb') as stream: magic,nm,ns,ncol=struct.unpack('<4q',stream.read(32))
    assert magic==20261017 and ncol==9 and path.stat().st_size==32+nm*72+ns*24
    modes=np.memmap(path,mode='r',dtype='<f8',offset=32,shape=(nm,9))
    support=np.memmap(path,mode='r',dtype='<f8',offset=32+nm*72,shape=(ns,3))
    owner=modes[:,0].astype(int);start=modes[:,2].astype(int);length=modes[:,3].astype(int);ids=support[:,0].astype(int)
    assert np.all(modes[:,4:6]==0) and np.all(length>0)
    assert np.array_equal(start,np.r_[0,np.cumsum(length[:-1])]) and length.sum()==ns
    mass=state[:,5];U=state[:,8];v=state[:,6:8];x=state[:,3:5]-3
    r=np.linalg.norm(x,axis=1);assert np.all(r>0)
    frozen=r>=2.5;free=~frozen[ids]
    alpha=np.add.reduceat(np.sum(support[:,1:3]**2,axis=1)/mass[ids]*free,start)
    dot=np.sum(support[:,1:3]*v[ids],axis=1);res=np.add.reduceat(dot,start)
    C=modes[:,6]*(modes[:,7]*abs(res)+np.sqrt((modes[:,7]*res)**2+modes[:,8]*U[owner]));C=np.where(alpha>0,C,0.)
    assert np.isfinite(C).all() and np.all(C>=0)
    hg=np.zeros((n,2));np.add.at(hg,ids,-np.repeat(C*res,length)[:,None]*support[:,1:3]*free[:,None])
    heat=float(np.sum(C[~frozen[owner]]*res[~frozen[owner]]**2))
    boundary_work=float(np.sum(C*res*np.add.reduceat(dot*(~free),start)))
    boundary_heat=float(np.sum(C[frozen[owner]]*res[frozen[owner]]**2))
    kinetic=float(np.sum(v*hg));closure=kinetic+heat-boundary_work+boundary_heat
    assert abs(closure)<1e-9*max(heat,abs(kinetic),1e-30)
    k=gravity_coefficient(r);omega=np.sqrt(k);z=r/.5
    low=z<=.5;mid=(z>.5)&(z<1);dM=np.zeros_like(z)
    dM[low]=32*z[low]**2-192*z[low]**4+192*z[low]**5
    dM[mid]=64*z[mid]**2-192*z[mid]**3+192*z[mid]**4-64*z[mid]**5
    slope=.5*(z*dM/(k*r**3)-3)
    jx=np.column_stack((-x[:,1],x[:,0]));jv=np.column_stack((-v[:,1],v[:,0]))
    reference=omega[:,None]*(jv+(slope*np.sum(x*v,axis=1)/r**2)[:,None]*jx)
    def circular(pos):
        return np.sqrt(gravity_coefficient(np.linalg.norm(pos,axis=1)))[:,None]*np.column_stack((-pos[:,1],pos[:,0]))
    errors=[]
    for h in (1e-5,5e-6):
        fd=(circular(x+h*v)-circular(x-h*v))/(2*h)
        errors.append(float(np.linalg.norm(fd-reference)/np.linalg.norm(reference)))
    assert errors[-1]<1e-7
    dv=v-circular(x);gravity=-k[:,None]*x
    row={'time':t,'HG_energy_closure':closure,'reference_derivative_checks':errors,'regions':{}}
    for label,mask in masks.items():
        def power(force):return float(np.sum(dv[mask]*force[mask]))
        p=pressure['cases'][name]['regions'][label]['actual']
        hyd=power(rhs[:,8:10]);hp=power(hg);gp=power(mass[:,None]*gravity);rp=-power(mass[:,None]*reference)
        row['regions'][label]={'pressure':p,'hydro_minus_pressure':hyd-p,'HG':hp,'gravity':gp,'reference_motion':rp,'gravity_plus_reference':gp+rp,'total':hyd+hp+gp+rp}
    out['cases'][name]=row
path=root/'complete_departure_work.json'
if path.exists():assert json.loads(path.read_text())==out
else:
    with path.open('x') as stream:json.dump(out,stream,indent=2)
print(json.dumps(out,indent=2))
