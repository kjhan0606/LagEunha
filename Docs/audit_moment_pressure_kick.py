"""Frozen-geometry shared pressure work and finite affine velocity kick.

Actual remapped t2 moments. Constant cell pressure, arithmetic face pressure,
affine reconstructed velocity; patch exterior is prescribed checkpoint data.
This is an experimental pressure substep, not a full Euler discretization:
there is no geometric drift, advection, gravity, shock solver or weight ODE.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def decode(q):
    m=q[:,0];c=q[:,1:3]/m[:,None];v=q[:,6:8]/m[:,None]
    Q=np.array([[q[:,3],q[:,4]],[q[:,4],q[:,5]]]).transpose(2,0,1)
    S=Q-m[:,None,None]*c[:,:,None]*c[:,None,:]
    C=q[:,8:12].reshape(-1,2,2)-m[:,None,None]*c[:,:,None]*v[:,None,:]
    B=np.linalg.solve(S,C).transpose(0,2,1)
    T=np.array([[q[:,12],q[:,13]],[q[:,13],q[:,14]]]).transpose(2,0,1)-m[:,None,None]*v[:,:,None]*v[:,None,:]
    R=T-B@C;U=q[:,15]-.5*(q[:,12]+q[:,14]);return c,v,S,C,B,R,U

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();base=Path('/gpfs/kjhan/LagForce')
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py';geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py';readpath=base/'Codex_paired_response_20261017/run.py'
    mod=load('moments',modpath);geo=load('geo',geopath);reader=load('reader',readpath)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat';layout=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,t=reader.read(source,json.loads(layout.read_text()));a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));vel=np.column_stack((a['vx'],a['vy']))
    targetpath=base/'Codex_evolving_weight_transport_20261017/weights_0.npz';target=np.load(targetpath)['weights']
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.mk_batch,gpu.mk_cuda_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    def call(f,data):
        out=np.empty((len(data),17));assert f(len(data),np.ascontiguousarray(data),out)==0;assert np.isfinite(out).all();return out
    rows=[];paths=[source,layout,modpath,geopath,readpath,targetpath,Path(__file__),root/'moment_affine_recovery.h',root/'moment_affine_kick.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so']
    for fraction in (.25,1.):
        statepath=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(statepath)
        saved=np.load(statepath);ids=saved['ids'];q=saved['state'];V=saved['volume'];w=fraction*target
        c,v,S,C,B,R,U=decode(q);gamma=1.6666669845581055;P=(gamma-1)*U/V
        polygons=mod.polygons(x,w,ids,geo.clip);lookup=np.full(len(a),-1,int);lookup[ids]=np.arange(len(ids))
        hull=ConvexHull(np.column_stack((x,np.sum(x*x,axis=1)-w)));tri=hull.simplices[hull.equations[:,2]<-1e-12]
        edges=np.unique(np.sort(np.concatenate((tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]])),axis=1),axis=0)
        neighbors=[[] for _ in x]
        for i,j in edges:neighbors[i].append(j);neighbors[j].append(i)
        faces=[]
        for k,(i,p) in enumerate(zip(ids,polygons)):
            js=np.array(neighbors[i]);delta=x[js]-x[i];bound=np.sum(delta*delta,axis=1)+w[i]-w[js]
            for aa,bb in zip(p,np.roll(p,-1,axis=0)):
                d=bb-aa;length=np.linalg.norm(d)
                if length<1e-13:continue
                mid=(aa+bb)/2;error=abs(2*delta@(mid-x[i])-bound);best=np.argmin(error);assert error[best]<1e-9
                j=js[best];other=lookup[j]
                if other>=0 and i>j:continue
                normal=np.array([d[1],-d[0]])/length
                faces.append((k,other,int(j),mid,normal,length))
        def forces(pressure,ghost_pressure):
            F=np.zeros((len(q),2));J=np.zeros((len(q),2,2));boundary=np.zeros(3)
            for i,j,g,point,n,length in faces:
                pf=.5*(pressure[i]+(pressure[j] if j>=0 else ghost_pressure[g]));force=-pf*length*n
                F[i]+=force;J[i]+=np.outer(point,force)
                if j>=0:F[j]-=force;J[j]-=np.outer(point,force)
                else:boundary+=np.r_[force,point[0]*force[1]-point[1]*force[0]]
            J+=pressure[:,None,None]*V[:,None,None]*np.eye(2)
            return F,J,boundary
        F,J,boundary=forces(P,a['pressure'])
        uniformF,uniformJ,_=forces(np.ones(len(q)),np.ones(len(a)))
        assert np.max(abs(uniformF))<1e-10 and np.max(abs(uniformJ))<1e-10
        # Acceleration inferred from force and its first spatial moment.
        centered=J-c[:,:,None]*F[:,None,:];A=np.linalg.solve(S,centered).transpose(0,2,1);acc=F/q[:,0,None]
        for dt in (.001,.0005,.00025):
            dE=np.zeros(len(q));boundaryE=0.
            for i,j,g,point,n,length in faces:
                vi=v[i]+B[i]@(point-c[i])+.5*dt*(acc[i]+A[i]@(point-c[i]))
                vj=v[j]+B[j]@(point-c[j])+.5*dt*(acc[j]+A[j]@(point-c[j])) if j>=0 else vel[g]
                pf=.5*(P[i]+(P[j] if j>=0 else a['pressure'][g]));work=-dt*pf*length*np.dot(.5*(vi+vj),n)
                dE[i]+=work
                if j>=0:dE[j]-=work
                else:boundaryE+=work
            data=np.column_stack((q,F,J.reshape(-1,4),np.full(len(q),dt),dE))
            assert cpu.mk_threads(1)==1;one=call(cpu.mk_batch,data);assert cpu.mk_threads(4)==4;four=call(cpu.mk_batch,data);device=call(gpu.mk_cuda_batch,data)
            assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
            _,vn,_,Cn,Bn,Rn,Un=decode(one);change=one-q
            # Independent covariance congruence, including the quadratic term.
            T=np.array([[q[:,12],q[:,13]],[q[:,13],q[:,14]]]).transpose(2,0,1)-q[:,0,None,None]*v[:,:,None]*v[:,None,:]
            Tnew=T+dt*(A@C+(A@C).transpose(0,2,1))+dt*dt*A@S@A.transpose(0,2,1)
            Pinew=Tnew+q[:,0,None,None]*vn[:,:,None]*vn[:,None,:]
            oracle=np.column_stack((Pinew[:,0,0],Pinew[:,0,1],Pinew[:,1,1]));oracle_error=float(np.max(abs(oracle-one[:,12:15])))
            energy_error=float(change[:,15].sum()-boundaryE)
            angular_error=float((change[:,9]-change[:,10]).sum()-dt*boundary[2])
            residual_error=float(np.max(abs(Rn-R)/q[:,0,None,None]))
            assert oracle_error<1e-12 and abs(energy_error)<1e-12 and abs(angular_error)<1e-12 and np.all(Un>0)
            row={'weight_fraction':fraction,'dt':dt,'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'energy_boundary_ledger_error':energy_error,'angular_boundary_ledger_error':angular_error,
                 'momentum_boundary_ledger_error':(change[:,6:8].sum(axis=0)-dt*boundary[:2]).tolist(),
                 'finite_covariance_oracle_max_error':oracle_error,'R_change_per_mass_max_absolute':residual_error,
                 'minimum_thermal_ratio':float(np.min(Un/U)),'thermal_change_over_initial_U':float((Un-U).sum()/U.sum()),
                 'uniform_pressure_force_max':float(np.max(abs(uniformF))),'uniform_pressure_first_moment_force_max':float(np.max(abs(uniformJ))),
                 'mass_position_potential_unchanged':bool(np.array_equal(one[:,:6],q[:,:6]) and np.array_equal(one[:,16],q[:,16]))}
            rows.append(row);print(json.dumps(row),flush=True)
    result={'scope':'Finite pressure kick on frozen actual remapped t2 geometry,1590 cells; shared face pressure work plus affine acceleration preserving residual covariance algebraically. Exterior checkpoint pressure/velocity prescribed.',
            'source_time':t,'rows':rows,'limits':['Not the old production pressure force; no adoption or zero-weight reduction claim.', 'No drift, mass remap, gravity, shock treatment, physical-time weight law or MPI in this diagnostic.', 'Face pressures frozen over each kick. Geometric/transport coupling and full time-order remain unverified.', 'R inherited from rounded raw moments can have tiny negative eigenvalues; the kick preserves R algebraically, not a general numerical positivity guarantee.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
