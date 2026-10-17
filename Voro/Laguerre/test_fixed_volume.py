"""3D geometry/capacity gate only; prescribed sites, no Euler evolution."""
import ctypes,json,time
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection,ConvexHull

lib=ctypes.CDLL(str(Path('libgate.so').resolve()))
array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib.cell.argtypes=[ctypes.c_int,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_double)];lib.cell.restype=ctypes.c_int
normals=np.vstack([np.eye(3),-np.eye(3)])


def cells(sites,weights,independent=False):
    result=[];vol=[];reference_error=0.
    for i,site in enumerate(sites):
        ids=np.flatnonzero(np.arange(len(sites))!=i);delta=sites[ids]-site
        distances=np.r_[1-site,site]
        ghosts=2*distances[:,None]*normals
        data=np.ascontiguousarray(np.vstack([np.column_stack([delta,weights[ids]]),np.column_stack([ghosts,np.full(6,weights[i])])]))
        out=np.zeros(1536);volume=ctypes.c_double();count=lib.cell(len(data),weights[i],data,out,ctypes.byref(volume))
        if count<4:raise ValueError('unsupported or empty 3D cell '+str((i,count)))
        vertices=out[:3*count].reshape(-1,3);vol.append(volume.value);result.append(vertices+site)
        if independent:
            a=data[:,:3];bound=.5*(np.sum(a*a,axis=1)+weights[i]-data[:,3])
            interior=vertices.mean(0)
            hull=HalfspaceIntersection(np.column_stack([a,-bound]),interior).intersections
            reference_error=max(reference_error,abs(ConvexHull(hull).volume-volume.value))
    return np.asarray(vol),result,reference_error


def jacobian(sites,weights,polys):
    N=len(sites);J=np.zeros((N,N));areas={}
    for i,poly in enumerate(polys):
        for j in range(N):
            if i==j:continue
            delta=sites[j]-sites[i];d=np.linalg.norm(delta);bound=.5*(delta@delta+weights[i]-weights[j]);r=poly-sites[i]
            face=r[abs(r@delta-bound)<2e-10]
            if len(face)<3:continue
            normal=delta/d;axis=np.eye(3)[np.argmin(abs(normal))];u=np.cross(normal,axis);u/=np.linalg.norm(u);v=np.cross(normal,u)
            area=ConvexHull(np.column_stack([face@u,face@v])).volume
            J[i,j]=-area/(2*d);J[i,i]+=area/(2*d);areas[i,j]=area
    reciprocal=max((abs(a-areas.get((j,i),0.)) for (i,j),a in areas.items()),default=0.)
    return J,reciprocal


def run(n):
    rng=np.random.default_rng(762+n);grid=np.meshgrid(*([(np.arange(n)+.5)/n]*3),indexing='ij');sites0=np.column_stack([a.ravel() for a in grid]);sites0+=rng.uniform(-.1,.1,sites0.shape)/n
    weights=np.zeros(len(sites0));target,_,ref0=cells(sites0,weights,True);assert abs(target.sum()-1)<1e-12
    velocity=np.zeros_like(sites0);velocity[:,0]=.05*np.sin(2*np.pi*sites0[:,1]);velocity[:,1]=.03*np.sin(2*np.pi*sites0[:,2])
    rows=[]
    for t in (.01,.05,.1):
        sites=sites0+t*velocity
        for iteration in range(12):
            volumes,polys,_=cells(sites,weights);residual=volumes-target
            if np.max(abs(residual)/target)<1e-10:break
            J,reciprocal=jacobian(sites,weights,polys);assert reciprocal<1e-10
            step=np.r_[np.linalg.solve(J[:-1,:-1],-residual[:-1]),0.]
            for power in range(15):
                trial=weights+step*2.**(-power)
                try:new,_,_=cells(sites,trial)
                except ValueError:continue
                if np.max(abs(new-target))<np.max(abs(residual)):weights=trial;break
            else:raise RuntimeError('capacity line search failed')
        else:raise RuntimeError('capacity iteration exhausted')
        volumes,polys,reference=cells(sites,weights,True);J,reciprocal=jacobian(sites,weights,polys)
        relative=float(np.max(abs(volumes-target)/target));assert relative<1e-10 and abs(volumes.sum()-1)<1e-12 and reference<1e-11
        # Independent directional derivative of the volume map.
        direction=rng.normal(size=len(weights));direction-=direction[-1];h=1e-7
        plus,_,_=cells(sites,weights+h*direction);minus,_,_=cells(sites,weights-h*direction)
        derivative_error=float(np.max(abs((plus-minus)/(2*h)-J@direction)))
        assert derivative_error<2e-8
        rows.append(dict(t=t,newton_iterations=iteration,max_relative_mass_error=relative,partition_error=float(abs(volumes.sum()-1)),independent_volume_error=reference,shared_face_area_mismatch=reciprocal,weight_derivative_error=derivative_error))
    return dict(n=n,cells=n**3,initial_partition_error=float(abs(target.sum()-1)),initial_reference_error=ref0,stages=rows)

if __name__=='__main__':
    start=time.monotonic();results=[run(n) for n in (3,4)]
    Path('fixed_volume_3d_report.json').write_text(json.dumps(dict(passed=True,results=results,seconds=time.monotonic()-start,scope='Actual C 3D geometry with Newton capacity solve at prescribed unsteered positions; not a fluid simulation.'),indent=2))
    print('PASS fixed-volume 3D geometry',flush=True)
