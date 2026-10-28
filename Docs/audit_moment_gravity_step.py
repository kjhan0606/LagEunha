"""Finite coupled kick/drift of the cheap quadrupole gravity component.
This is material gravity evolution, NOT a Kepler fluid disk simulation.
"""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();base=Path('/gpfs/kjhan/LagForce')
    cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'));source=ctypes.CDLL(str(base/'Codex_quadrupole_gravity_20261018/cpu.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.ms_batch,gpu.ms_cuda_batch,source.mq_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    def run(f,q,dt,steps):
        q=q.copy();data=np.empty((len(q),19));data[:,17]=.5;data[:,18]=dt;out=np.empty_like(q)
        for _ in range(steps):
            data[:,:17]=q;assert f(len(q),data,out)==0;assert np.isfinite(out).all();q,out=out,q
        return q
    rows=[];paths=[Path(__file__),root/'moment_gravity_step.h',root/'moment_affine_kick.h',root/'moment_material_drift.h',root/'moment_quadrupole_gravity.h',root/'moment_affine_recovery.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so',base/'Codex_quadrupole_gravity_20261018/cpu.so']
    for fraction in (.25,1.):
        path=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(path);original=np.load(path)['state'];q=original.copy()
        g=np.empty((len(q),12));source.mq_threads(1);assert source.mq_batch(len(q),np.ascontiguousarray(np.column_stack((q,np.full(len(q),.5)))),g)==0
        q[:,16]=g[:,7] # Explicit diagnostic Hamiltonian initialization; original file preserved.
        U0=q[:,15]-.5*(q[:,12]+q[:,14]);E0=float((q[:,15]+q[:,16]).sum());L0=q[:,9]-q[:,10]
        errors=[]
        for dt in (.001,.0005,.00025):
            steps=round(.2/dt);assert cpu.ms_threads(1)==1;one=run(cpu.ms_batch,q,dt,steps);assert cpu.ms_threads(4)==4;four=run(cpu.ms_batch,q,dt,steps);device=run(gpu.ms_cuda_batch,q,dt,steps)
            assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
            U=one[:,15]-.5*(one[:,12]+one[:,14]);assert np.all(U>0)
            drift=float((one[:,15]+one[:,16]).sum()-E0);errors.append(abs(drift))
            reverse=run(cpu.ms_batch,one,-dt,steps);rev=float(np.max(abs(reverse-q)));assert rev<1e-11
            row={'fraction':fraction,'dt':dt,'horizon':.2,'steps':steps,'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,
                 'approximate_Hamiltonian_energy_change':drift,'energy_change_over_initial_U':drift/U0.sum(),
                 'max_cell_angular_change':float(np.max(abs(one[:,9]-one[:,10]-L0))),
                 'thermal_change_max_absolute':float(np.max(abs(U-U0))),'thermal_change_over_total_U':float((U-U0).sum()/U0.sum()),
                 'reverse_max_absolute_error':rev,'mass_unchanged':bool(np.array_equal(one[:,0],q[:,0]))}
            rows.append(row);print(json.dumps(row),flush=True)
        ratios=[errors[k]/errors[k+1] for k in range(2)];assert all(3.5<z<4.5 for z in ratios),ratios
    result={'scope':'Gravity-only material moment evolution for actual remapped t2 state, horizon .2. KDK of the truncated quadrupole Hamiltonian, no pressure, mesh evolution, weight evolution, remap or MPI.',
            'checkpoint_conversion':'Only in-memory diagnostic potential replaced by explicitly defined quadrupole Psi. Original checkpoints and their integrated potential unchanged. Not an approved production conversion.',
            'rows':rows,'limits':['Second-order component energy convergence does not prove full-fluid accuracy or long-time Kepler stability.', 'No zero-weight reduction or general positivity guarantee.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
