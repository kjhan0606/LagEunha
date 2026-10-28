"""Existing positive moment compression on selected actual remapped cells.

Compare preserving 17 state integrals alone with also preserving 9 independent
gravity-source entries. A compact atomic measure still is not a fluid closure.
"""
import ctypes
import importlib.util
import json
from pathlib import Path
import numpy as np


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def compress(F,weight):
    # Same positive nullspace elimination principle as the 2026-10-09
    # subcell_quadrature.py, generalized to extra source constraints.
    count=F.shape[0];active=np.arange(len(weight));w=weight.copy()
    while len(active)>count:
        idx=active[:count+1];f=F[:,idx];scale=np.max(abs(f),axis=1);scale[scale==0]=1
        _,_,vh=np.linalg.svd(f/scale[:,None],full_matrices=True)
        z=vh[-1].astype(np.longdouble);ww=w[idx].astype(np.longdouble)
        pos=np.flatnonzero(z>0)
        if not len(pos):z=-z;pos=np.flatnonzero(z>0)
        pivot=pos[np.argmin(ww[pos]/z[pos])]
        new=(ww*z[pivot]-ww[pivot]*z)/z[pivot];new[pivot]=0
        assert np.all(new>=0)
        w[idx]=new.astype(float);active=active[w[active]>0]
    return active,w[active]


def main():
    root=Path(__file__).parent;base=Path('/gpfs/kjhan/LagForce')
    mod=load('moments',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
    geo=load('geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
    reader=load('reader',base/'Codex_paired_response_20261017/run.py')
    kr=base/'Codex_checkpoint_weight_remap_20261018';gr=base/'Codex_moment_gravity_20261018'
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    _,a,sort,_=reader.read(source,json.loads((base/'Codex_paired_response_20261017/layout.json').read_text()));a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3));v=np.column_stack((a['vx'],a['vy']))
    fit=np.load(kr/'lift_coefficients.npz');ids=fit['ids'];vol=fit['volume'];rho=a['mass'][ids]/vol
    coeff=fit['coefficient']/fit['feature_scale'];mean=fit['cell_feature_mean']@fit['coefficient'];phi0=mod.potential(x[ids])
    saved=np.load(kr/'state_w1_q8.npz')['state']
    variance=.5*(saved[:,12]+saved[:,14])-.5*np.sum(saved[:,6:8]**2,axis=1)/saved[:,0]
    heat=saved[:,15]-.5*(saved[:,12]+saved[:,14])
    chosen=np.argsort(variance/heat)[-6:]
    old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip)
    w=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
    new=mod.polygons(x,w,ids[chosen],geo.clip)
    lo=np.array([p.min(axis=0) for p in old]);hi=np.array([p.max(axis=0) for p in old])
    ml=ctypes.CDLL(str(kr/'cpu.so'));mg=ctypes.CDLL(str(kr/'cuda.so'))
    gl=ctypes.CDLL(str(gr/'cpu.so'));gg=ctypes.CDLL(str(gr/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for fun in (ml.wm_batch,mg.wm_cuda_batch,gl.wg_batch,gg.wg_cuda_batch):
        fun.argtypes=[ctypes.c_size_t,arr,arr];fun.restype=ctypes.c_int
    ml.wm_set_threads(1);gl.wg_threads(1)
    def call(fun,data,width):
        data=np.ascontiguousarray(data);out=np.empty((len(data),width));assert fun(len(data),data,out)==0;return out
    rows=[]
    for local,receiver in zip(chosen,new):
        blocks=[];cover=0.
        candidates=np.flatnonzero(np.all(hi>=receiver.min(axis=0),axis=1)&np.all(lo<=receiver.max(axis=0),axis=1))
        for src in candidates:
            p=old[src].copy()
            for aa,bb in zip(receiver,np.roll(receiver,-1,axis=0)):
                edge=bb-aa;normal=np.array([edge[1],-edge[0]])
                p=geo.clip(p,normal,normal@aa)
                if len(p)<3:break
            if mod.area(p)<=1e-18:continue
            cover+=mod.area(p);xy,dv=mod.quadrature(p,8);r=xy-x[ids[src]];vel=v[ids[src]]
            feature=r[:,0]*coeff[0]+r[:,1]*coeff[1]+(r[:,0]*vel[1]-r[:,1]*vel[0])*coeff[2]+(mod.potential(xy)-phi0[src])*coeff[3]-mean[src]
            density=rho[src]*(1+rho[src]*feature);dm=density*dv
            assert np.all(dm>0)
            e=(a['ie'][ids[src]]/vol[src])/density
            blocks.append(np.column_stack((dm,xy,np.tile(vel,(len(xy),1)),e,mod.potential(xy))))
        assert abs(cover/mod.area(receiver)-1)<1e-8
        atoms=np.concatenate(blocks);weight=atoms[:,0].copy();unit=atoms.copy();unit[:,0]=1
        state_unit=call(ml.wm_batch,unit,17)
        gd=np.column_stack((unit[:,:5],np.full(len(unit),.5)));gravity_unit=call(gl.wg_batch,gd,20)[:,:17]
        state0=weight@state_unit;gravity0=weight@gravity_unit
        for label,F in [('state17',state_unit.T),('state_and_gravity26',np.column_stack((state_unit,gravity_unit[:,[6,7,8,9,10,11,12,13,14]])).T)]:
            select,weights=compress(F,weight);small=atoms[select].copy();small[:,0]=weights
            qcpu=call(ml.wm_batch,small,17);qgpu=call(mg.wm_cuda_batch,small,17)
            ggdata=np.column_stack((small[:,:5],np.full(len(small),.5)))
            gcpu=call(gl.wg_batch,ggdata,20);ggpu=call(gg.wg_cuda_batch,ggdata,20)
            assert np.array_equal(qcpu.view('u8'),qgpu.view('u8'))
            assert np.array_equal(gcpu.view('u8'),ggpu.view('u8'))
            state=qcpu.sum(axis=0);gravity=gcpu[:,:17].sum(axis=0)
            error=np.max(abs(state-state0)/(1+abs(state0)));assert error<1e-12
            source_error=float(np.max(abs(gravity-gravity0)/(1+abs(gravity0))))
            if label=='state_and_gravity26':assert source_error<1e-12
            rows.append({'id':int(ids[local]),'mode':label,'original_atoms':len(atoms),'retained_atoms':len(weights),
                         'state_max_scaled_error':float(error),'gravity_source_max_scaled_error':source_error,
                         'gravity_force_relative_error':float(np.linalg.norm(gravity[6:8]-gravity0[6:8])/np.linalg.norm(gravity0[6:8])),
                         'gravity_power_error_over_initial_cell_thermal':float((gravity[15]-gravity0[15])/heat[local]),
                         'potential_energy_error':float(state[16]-state0[16]),'minimum_retained_mass':float(weights.min()),
                         'CPU_CUDA_bitwise_equal':True})
            print(json.dumps(rows[-1]),flush=True)
    out={'scope':'Six actual remapped cells with largest unresolved-kinetic/thermal ratio. Existing positive compression principle; source constraints extended. Not whole-patch performance or physical-time closure.',
         'rows':rows,'history':'Docs/polybnd/subcell_quadrature.py (2026-10-09)',
         'interpretation':'Preserving potential energy at one time does not preserve its force/work derivatives. Additional source features can preserve these sampled derivatives during compression, but do not define future spatial reconstruction or fluid evolution.'}
    (root/'result.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
