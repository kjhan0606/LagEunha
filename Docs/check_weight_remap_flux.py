"""Audit a shared C/CUDA weight-remap pair kernel against algebra and real faces."""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import numpy as np


def main():
    root=Path(__file__).parent
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.wr_flux_batch,gpu.wr_cuda_batch):
        f.argtypes=[ctypes.c_int,ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    cpu.wr_set_threads.argtypes=[ctypes.c_int];cpu.wr_set_threads.restype=ctypes.c_int
    def call(f,mode,data):
        data=np.ascontiguousarray(data,dtype=np.float64);out=np.empty((len(data),4))
        status=f(mode,len(data),data,out);assert status==0,('backend error',status)
        assert np.isfinite(out).all()
        return out
    source=Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator/rhs.dat')
    a=np.loadtxt(source);a=a[np.argsort(a[:,0])];n=len(a)
    x=a[:,1:3]-3;v=a[:,3:5];r=np.linalg.norm(x,axis=1)
    omega=(x[:,0]*v[:,1]-x[:,1]*v[:,0])/r**2;omega[r>=2.5]=0
    stamp=.001;c,s=np.cos(omega*stamp),np.sin(omega*stamp)
    pos=np.column_stack((c*x[:,0]-s*x[:,1],s*x[:,0]+c*x[:,1]));vel=np.column_stack((-omega*pos[:,1],omega*pos[:,0]))
    prior=Path('/gpfs/kjhan/LagForce/Codex_evolving_weight_transport_20261017')
    left=np.load(prior/'weights_0.npz');right=np.load(prior/'weights_0.001.npz')
    ids=right['ids'];w=right['weights'];wd=(w-left['weights'])/stamp
    code=Path('/gpfs/kjhan/LagForce/Codex_weighted_compression_geometry_20261017/audit.py')
    spec=importlib.util.spec_from_file_location('geometry_audit',code);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    volume,_,lap,_=module.geometry(pos,w,vel,ids)
    sparse=lap.tocoo();i=ids[sparse.row];j=sparse.col;take=i<j
    i,j,conductance=i[take],j[take],-sparse.data[take]
    rho=.01;pressure=2.7777796321446218e-5;gamma=1.6666669845581055
    rate=rho*conductance*(wd[i]-wd[j]);e=pressure/((gamma-1)*rho)
    actual=np.column_stack((rate,vel[i],vel[j],np.full(len(i),e),np.full(len(i),e)))
    rng=np.random.default_rng(94271)
    synthetic=np.column_stack((rng.uniform(-1,1,2048),rng.normal(size=(2048,4)),np.exp(rng.normal(size=(2048,2)))))
    reports={}
    for name,data in (('actual_weight_faces',actual),('synthetic_states',synthetic)):
        reports[name]={}
        for mode,label in ((0,'donor'),(1,'centered_kinetic_compatible'),(2,'kinetic_compatible_thermal_donor')):
            assert cpu.wr_set_threads(1)==1
            one=call(cpu.wr_flux_batch,mode,data)
            assert cpu.wr_set_threads(4)==4
            four=call(cpu.wr_flux_batch,mode,data)
            device=call(gpu.wr_cuda_batch,mode,data)
            assert np.array_equal(one.view('u8'),four.view('u8'))
            assert np.array_equal(one.view('u8'),device.view('u8'))
            fm=data[:,0];vi=data[:,1:3];vj=data[:,3:5]
            if mode==0:
                donor_v=np.where((fm>=0)[:,None],vj,vi);donor_e=np.where(fm>=0,data[:,6],data[:,5])
                oracle=np.column_stack((fm,fm[:,None]*donor_v,fm*(donor_e+.5*np.sum(donor_v**2,axis=1))))
            else:
                thermal=.5*(data[:,5]+data[:,6]) if mode==1 else np.where(fm>=0,data[:,6],data[:,5])
                oracle=np.column_stack((fm,fm[:,None]*.5*(vi+vj),fm*(thermal+.5*np.sum(vi*vj,axis=1))))
            assert np.max(abs(one-oracle)/(1+abs(oracle)))<2e-15
            reverse=data.copy();reverse[:,0]*=-1;reverse[:,1:3]=vj;reverse[:,3:5]=vi
            reverse[:,5]=data[:,6];reverse[:,6]=data[:,5]
            assert np.max(abs(call(cpu.wr_flux_batch,mode,reverse)+one))<1e-13
            boost=np.array([17.,-9.]);shifted=data.copy();shifted[:,1:3]+=boost;shifted[:,3:5]+=boost
            moved=call(cpu.wr_flux_batch,mode,shifted)
            expected=one.copy();expected[:,1:3]+=fm[:,None]*boost
            expected[:,3]+=one[:,1:3]@boost+.5*np.dot(boost,boost)*fm
            boost_error=float(np.max(abs(moved-expected)/(1+abs(expected))))
            assert boost_error<1e-13
            zero=data.copy();zero[:,0]=0
            assert np.all(call(gpu.wr_cuda_batch,mode,zero)==0)
            row={'faces':len(data),'openmp_1_4_bitwise_equal':True,'cpu_cuda_bitwise_equal':True,
                 'oracle_max_abs_difference':float(np.max(abs(one-oracle))),
                 'galilean_max_scaled_difference':boost_error,'reverse_flux_check':True,'zero_mass_rate_zero_flux':True}
            if name=='actual_weight_faces':
                total=np.zeros((n,4));np.add.at(total,i,one);np.add.at(total,j,-one)
                mdot=total[:,0];pdot=total[:,1:3];Edot=total[:,3]
                Kdot=np.sum(vel*pdot,axis=1)-.5*np.sum(vel*vel,axis=1)*mdot
                Udot=Edot-Kdot
                mixing=.5*np.sum(abs(rate)*np.sum((vel[i]-vel[j])**2,axis=1))
                reference=-mixing if mode==0 else 0.
                assert abs(Kdot.sum()-reference)<1e-15
                if mode in (1,2):assert np.max(abs(Udot-e*mdot))<1e-15
                torque=float(np.sum(pos[:,0]*pdot[:,1]-pos[:,1]*pdot[:,0]))
                point_L=float(np.sum(rho*volume*(pos[ids,0]*vel[ids,1]-pos[ids,1]*vel[ids,0])))
                row['semidiscrete']={'global_mass_momentum_energy_rates':total.sum(axis=0).tolist(),
                    'kinetic_energy_rate':float(Kdot.sum()),'thermal_energy_rate':float(Udot.sum()),
                    'donor_mixing_reference':float(mixing),'uniform_specific_thermal_transport_max_residual':float(np.max(abs(Udot-e*mdot))),
                    'point_sample_angular_momentum_rate':torque,'point_sample_torque_over_patch_L':torque/point_L}
            reports[name][label]=row
    contacts=[]
    for ratio in (1.,100.):
        for sweep in (-.1,.1):
            for shear in (0.,.1):
                rhoi,rhoj=ratio,1.;ei=1./((gamma-1)*rhoi);ej=1./((gamma-1)*rhoj)
                vi=np.array([1.,.2]);vj=np.array([1.,.2+shear])
                flux=sweep*(rhoj if sweep>=0 else rhoi)
                data=np.array([[flux,*vi,*vj,ei,ej]])
                residual={}
                for mode in (0,1,2):
                    out=call(gpu.wr_cuda_batch,mode,data)[0]
                    ki=vi@out[1:3]-.5*(vi@vi)*out[0]
                    kj=-vj@out[1:3]+.5*(vj@vj)*out[0]
                    ui=out[3]-ki;uj=-out[3]-kj
                    prate=np.array([(gamma-1)*ui-sweep,(gamma-1)*uj+sweep])
                    residual[str(mode)]=float(np.max(abs(prate)))
                    if mode==2 or (mode==0 and shear==0):assert max(abs(prate))<1e-12
                contacts.append({'density_ratio':ratio,'weight_swept_volume_rate':sweep,'tangential_velocity_jump':shear,
                                 'max_pressure_rate_error_by_mode':residual})
    out={'status':'Standalone experimental face kernel only; not connected to production integrator.',
         'device':subprocess.check_output(['nvidia-smi','--query-gpu=name,compute_cap','--format=csv,noheader'],text=True).strip(),
         'reports':reports,'uniform_pressure_contact_checks':contacts,
         'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'weight_remap_flux.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so',Path(__file__))},
         'derivation':'dp_i=f*(vi+vj)/2, dm_i=f, opposite for j: sum Kdot=0. dE_i=f*(e_face+vi.vj/2) yields Udot_i=f*e_face. e_face=e_donor and f=rho_donor*swept_volume_rate preserve uniform pressure in pure weight remapping.',
         'primary_related_literature':'https://arxiv.org/abs/1209.4994',
         'scope_of_literature':'Related semidiscrete kinetic-energy preservation; this remap-only kernel is not the complete entropy-stable Euler flux from that paper.',
         'limitations':['No MPI owner exchange or full CUDA hydrodynamic integration is verified by this batch.',
                        'Central transport need not be positive, monotone or entropy stable.',
                        'Point-sampled angular momentum is not conserved by this flux; consistent cell moments/torque treatment remain necessary.',
                        'Finite time integration and gravity potential work are not included.',
                        'A zero-rate remap is inert, but the complete equal-weight/fixed-mass Voronoi limit is not implemented.']}
    (root/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
