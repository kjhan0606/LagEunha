"""Single-evaluation gravity from a moment Hamiltonian versus dense reference.
The reference quadrature is diagnostic only. Approximate potential is NOT
silently substituted for the saved physical potential energy.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();base=Path('/gpfs/kjhan/LagForce')
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py';geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py';readpath=base/'Codex_paired_response_20261017/run.py'
    mod=load('moments',modpath);geo=load('geo',geopath);reader=load('reader',readpath)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat';layout=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,t=reader.read(source,json.loads(layout.read_text()));a=a[sort];x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    fitpath=base/'Codex_checkpoint_weight_remap_20261018/lift_coefficients.npz';fit=np.load(fitpath);ids=fit['ids'];volume=fit['volume'];rho=a['mass'][ids]/volume
    coeff=fit['coefficient']/fit['feature_scale'];mean=fit['cell_feature_mean']@fit['coefficient'];phi0=mod.potential(x[ids]);polys=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    oldroot=base/'Codex_moment_gravity_20261018';gravity=ctypes.CDLL(str(oldroot/'cpu.so'))
    atomroot=base/'Codex_checkpoint_weight_remap_20261018';atoms=ctypes.CDLL(str(atomroot/'cpu.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.mq_batch,gpu.mq_cuda_batch,gravity.wg_batch,atoms.wm_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    gravity.wg_threads(1);atoms.wm_set_threads(1)
    def call(f,data,width):
        data=np.ascontiguousarray(data);out=np.empty((len(data),width));assert f(len(data),data,out)==0;assert np.isfinite(out).all();return out
    def evaluate(q):return call(cpu.mq_batch,np.column_stack((q,np.full(len(q),.5))),12)
    states=[];reference=[];atom_count=0
    for k,p in enumerate(polys):
        xy,dv=mod.quadrature(p,8);rel=xy-x[ids[k]];vel=v[ids[k]]
        feature=rel[:,0]*coeff[0]+rel[:,1]*coeff[1]+(rel[:,0]*vel[1]-rel[:,1]*vel[0])*coeff[2]+(mod.potential(xy)-phi0[k])*coeff[3]-mean[k]
        density=rho[k]*(1+rho[k]*feature);assert np.all(density>0);dm=dv*density
        thermal=a['ie'][ids[k]]/volume[k]/density
        ad=np.column_stack((dm,xy,np.tile(vel,(len(xy),1)),thermal,mod.potential(xy)))
        states.append(call(atoms.wm_batch,ad,17).sum(axis=0))
        gd=np.column_stack((dm,xy,np.tile(vel,(len(xy),1)),np.full(len(xy),.5)))
        reference.append(call(gravity.wg_batch,gd,20)[:,:17].sum(axis=0));atom_count+=len(xy)
    q=np.array(states);reference=np.array(reference)
    datasets=[('zero_weight_actual_lift',q)]
    paths=[Path(__file__),source,layout,fitpath,modpath,geopath,readpath,oldroot/'cpu.so',atomroot/'cpu.so',root/'moment_quadrupole_gravity.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so']
    for fraction in (.25,1.):
        path=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(path);datasets.append((f'weight_{fraction:g}',np.load(path)['state']))
    rows=[]
    for label,state in datasets:
        data=np.column_stack((state,np.full(len(state),.5)));assert cpu.mq_threads(1)==1;one=call(cpu.mq_batch,data,12);assert cpu.mq_threads(4)==4;four=call(cpu.mq_batch,data,12);device=call(gpu.mq_cuda_batch,data,12)
        assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
        torque=float(np.max(abs(one[:,3]-one[:,4])/state[:,0]));assert torque<1e-12
        differences=[]
        for eps in (1e-4,5e-5,2.5e-5):
            # Material tangent at fixed velocity: q'=p, Q'=D+D^T.
            direction=np.zeros_like(state);direction[:,1:3]=state[:,6:8];direction[:,3]=2*state[:,8];direction[:,4]=state[:,9]+state[:,10];direction[:,5]=2*state[:,11]
            fd=(evaluate(state+eps*direction)[:,7]-evaluate(state-eps*direction)[:,7])/(2*eps)
            work_error=float(np.max(abs(fd+one[:,6])/state[:,0]))
            # Translate each material cell, holding its central covariance fixed.
            force_errors=[]
            for axis in (0,1):
                plus=state.copy();minus=state.copy();m=state[:,0];centroid=state[:,1:3]/m[:,None]
                for sign,target in ((1,plus),(-1,minus)):
                    delta=np.zeros(2);delta[axis]=sign*eps;target[:,1:3]+=m[:,None]*delta
                    target[:,3]+=2*delta[0]*state[:,1]+m*delta[0]**2
                    target[:,4]+=delta[0]*state[:,2]+delta[1]*state[:,1]+m*delta[0]*delta[1]
                    target[:,5]+=2*delta[1]*state[:,2]+m*delta[1]**2
                forcefd=-(evaluate(plus)[:,7]-evaluate(minus)[:,7])/(2*eps)
                force_errors.append(float(np.max(abs(forcefd-one[:,axis])/m)))
            differences.append({'epsilon':eps,'work_identity_max_error_per_mass':work_error,'force_gradient_max_error_per_mass':max(force_errors)})
        assert differences[-1]['work_identity_max_error_per_mass']<1e-5 and differences[-1]['force_gradient_max_error_per_mass']<1e-5
        U=state[:,15]-.5*(state[:,12]+state[:,14]);row={'state':label,'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,'maximum_torque_per_mass':torque,'finite_difference':differences,
            'quadrupole_minus_saved_potential_over_total_U':float((one[:,7]-state[:,16]).sum()/U.sum()),
            'max_cell_potential_difference_over_cell_U':float(np.max(abs(one[:,7]-state[:,16])/U))}
        if label=='zero_weight_actual_lift':
            point=call(gravity.wg_batch,np.column_stack((a['mass'][ids],x[ids],v[ids],np.full(len(ids),.5))),20)
            mass=state[:,0];den=np.sqrt(np.average(np.sum((reference[:,6:8]/mass[:,None])**2,axis=1),weights=mass))
            def rms(F):return float(np.sqrt(np.average(np.sum(((F-reference[:,6:8])/mass[:,None])**2,axis=1),weights=mass))/den)
            row['force_error_relative_to_integrated_rms']={'old_generator_point':rms(point[:,6:8]),'moment_quadrupole':rms(one[:,:2])}
            row['gravity_power_difference_over_total_U']=float((one[:,6]-reference[:,15]).sum()/U.sum())
        rows.append(row);print(json.dumps(row),flush=True)
    result={'scope':'Single radial evaluation per cell for an approximate quadrupole Hamiltonian, actual t2 lifted/remapped states. Dense reference used only once for accuracy audit, not production.',
            'source_time':t,'cells':len(ids),'reference_atoms':atom_count,'rows':rows,
            'derivation':{'potential':'Psi=m phi(c)+0.5 H(c):S','force':'F=-m grad(phi)-0.5 grad(H):S','affine_acceleration':'A=-H','force_first_moment':'J=c F^T+S A^T','power':'W=F.vbar+trace(A C)','central_torque':'antisym(J)=0 for central radial phi; requires the quadrupole term in F'},
            'limits':['Hamiltonian truncated at second spatial moment; higher spatial moments not represented.', 'Approximate potential differs from the saved integrated potential. No silent checkpoint conversion or production adoption.', 'Zero spatial extent gives point gravity at the mass center, not automatically at the mesh generator; approved reduction still requires proof.', 'No coupled pressure/drift/remap integration or long-time Kepler stability claim.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
