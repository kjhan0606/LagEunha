"""Bounded-cost reconstruction of actual finite-remap states.

Spatial polynomial fit uses exact polynomial integration (degree <=4) on
polygons, not dense integration of gravity. It is a candidate, not a limiter:
negative polynomial densities are reported and never clipped or renormalized.
Velocity residual R is retained. No hydro time evolution is performed.
"""
import ctypes
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def basis(z):
    x,y=z.T
    return np.column_stack((np.ones(len(z)),x,y,x*x,x*y,y*y))

def minimum_quadratic(p,theta):
    # A quadratic minimum on a compact convex polygon lies on an edge or at
    # an interior stationary point. Degenerate flat minima reach an edge.
    h=np.array([[2*theta[3],theta[4]],[theta[4],2*theta[5]]])
    candidates=list(p)
    for aa,bb in zip(p,np.roll(p,-1,axis=0)):
        d=bb-aa;curvature=d@h@d
        if curvature>0:
            t=-d@(theta[1:3]+h@aa)/curvature
            if 0<t<1:candidates.append(aa+t*d)
    if abs(np.linalg.det(h))>1e-15*np.linalg.norm(h)**2:
        z=np.linalg.solve(h,-theta[1:3]);edges=np.roll(p,-1,axis=0)-p;rel=z-p
        if np.all(edges[:,0]*rel[:,1]-edges[:,1]*rel[:,0]>=-1e-13):candidates.append(z)
    return float(np.min(basis(np.array(candidates))@theta))

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists()
    base=Path('/gpfs/kjhan/LagForce');started=time.monotonic()
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py'
    geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py'
    readpath=base/'Codex_paired_response_20261017/run.py'
    mod=load('moments',modpath);geo=load('geo',geopath);reader=load('reader',readpath)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
    layout=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,t=reader.read(source,json.loads(layout.read_text()));a=a[sort]
    x=np.column_stack((a['x']-3,a['y']-3))
    targetpath=base/'Codex_evolving_weight_transport_20261017/weights_0.npz'
    target=np.load(targetpath)['weights'];cpu=ctypes.CDLL(str(root/'cpu.so'));gpu=ctypes.CDLL(str(root/'cuda.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
    for f in (cpu.ma_batch,gpu.ma_cuda_batch):f.argtypes=[ctypes.c_size_t,arr,arr];f.restype=ctypes.c_int
    def call(f,q):
        out=np.empty((len(q),21));assert f(len(q),q,out)==0;assert np.isfinite(out).all();return out
    rows=[];paths=[source,layout,readpath,modpath,geopath,targetpath,Path(__file__),root/'moment_affine_recovery.h',root/'probe.c',root/'probe.cu',root/'cpu.so',root/'cuda.so']
    for fraction in (.25,1.):
        path=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(path)
        saved=np.load(path);ids=saved['ids'];q=np.ascontiguousarray(saved['state'])
        cpu.ma_threads(1);one=call(cpu.ma_batch,q);cpu.ma_threads(4);four=call(cpu.ma_batch,q);device=call(gpu.ma_cuda_batch,q)
        assert np.array_equal(one.view('u8'),four.view('u8')) and np.array_equal(one.view('u8'),device.view('u8'))
        polygons=mod.polygons(x,fraction*target,ids,geo.clip)
        density_min=[];fit_error=[];condition=[];negative_ids=[];residual_eigs=[];angular_error=[];oracle_error=[]
        for k,(poly,z,o) in enumerate(zip(polygons,q,one)):
            m=z[0];S=np.array([[o[11],o[12]],[o[12],o[13]]]);C=o[14:18].reshape(2,2)
            T=np.array([[z[12],z[13]],[z[13],z[14]]])-np.outer(z[6:8],z[6:8])/m
            B=C.T@np.linalg.inv(S);R=T-B@C
            oracle_error.append(float(np.max(abs(B-o[4:8].reshape(2,2)))))
            residual_eigs.append(np.linalg.eigvalsh(np.array([[o[8],o[9]],[o[9],o[10]]]))/m)
            angular_error.append((S@o[4:8].reshape(2,2).T)[0,1]-(S@o[4:8].reshape(2,2).T)[1,0]-(C[0,1]-C[1,0]))
            # Center and scale before solving the 6x6 polynomial moment fit.
            origin=poly.mean(axis=0);V=mod.area(poly);h=np.sqrt(V)
            points,dv=mod.quadrature(poly,3);F=basis((points-origin)/h)
            G=F.T@(dv[:,None]*F)/V
            c=z[1:3]/m;Q=np.array([[z[3],z[4]],[z[4],z[5]]])/m
            centered=(Q-np.outer(origin,c)-np.outer(c,origin)+np.outer(origin,origin))/(h*h)
            b=np.r_[1.,(c-origin)/h,centered[0,0],centered[0,1],centered[1,1]]
            theta=np.linalg.solve(G,b);minimum=minimum_quadratic((poly-origin)/h,theta)
            density_min.append(minimum);condition.append(np.linalg.cond(G));fit_error.append(float(np.max(abs(G@theta-b))))
            if minimum < -1e-10:negative_ids.append(int(ids[k]))
        # The total kinetic content includes BOTH affine and residual motion.
        Kmean=.5*np.sum(q[:,6:8]**2,axis=1)/q[:,0];Ktotal=.5*(q[:,12]+q[:,14]);U=one[:,18]
        energy_residual=float(np.max(abs(Kmean+one[:,19]+one[:,20]-Ktotal)))
        assert energy_residual<1e-12 and max(oracle_error)<1e-8
        row={'weight_fraction':fraction,'cells':len(ids),'CPU_OMP1_OMP4_CUDA_bitwise_equal':True,
             'affine_B_numpy_max_absolute_error':max(oracle_error),'angular_recovery_max_absolute_error':float(np.max(np.abs(angular_error))),
             'kinetic_decomposition_max_absolute_error':energy_residual,'minimum_eigenvalue_R_per_mass':float(np.min(residual_eigs)),
             'affine_kinetic_over_total_thermal':float(one[:,19].sum()/U.sum()),'residual_kinetic_over_total_thermal':float(one[:,20].sum()/U.sum()),
             'residual_over_cell_thermal_quantiles':np.quantile(one[:,20]/U,[0,.5,.9,.99,1]).tolist(),
             'quadratic_density_minimum_factor':min(density_min),'negative_density_cells':len(negative_ids),'negative_density_ids':negative_ids,
             'quadratic_density_minimum_quantiles':np.quantile(density_min,[0,.01,.5,.99,1]).tolist(),
             'polynomial_moment_fit_max_error':max(fit_error),'max_polynomial_gram_condition':max(condition)}
        rows.append(row);print(json.dumps(row),flush=True)
        np.savez(root/f'reconstruction_w{fraction:g}.npz',ids=ids,affine=one,density_minimum_factor=density_min)
    bad=q[:1].copy();bad[:,0]=0;out=np.empty((1,21));assert cpu.ma_batch(1,bad,out)!=0 and gpu.ma_cuda_batch(1,bad,out)!=0
    result={'scope':'Actual t2 finite weight-remapped 17-moment states; compact static reconstruction only, no physical-time evolution, no force change, no MPI.',
            'source_time':t,'lineage':'2026-10-09 affine projection; new actual checkpoint and shared C/OpenMP/CUDA recovery. Quadratic density is a new bounded-cost candidate, not an adopted closure.',
            'density_integration':'Degree<=4 polynomial geometry moments, three-point tensor Gauss per centroid triangle. No subcell gravity evaluation.',
            'minimum_method':'Vertices, edge stationary minima, interior stationary point; floating arithmetic, not certified interval bounds.',
            'rows':rows,'invalid_mass_rejected_CPU_CUDA':True,'seconds':time.monotonic()-started,
            'limits':['Positive covariance alone does not prove spatial density feasibility.', 'R is unresolved kinetic covariance and cannot be silently discarded or heated. Its physical stress/energy transport remains open.', 'Quadratic density matches m,q,Q only; it does not preserve the nonlinear potential moment.', 'No long-term Kepler stability claim.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
