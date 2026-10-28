"""Stored-stage geometry audit of the actual Exam/exam_gpu.cu translation unit.

No RHS, massroot, evolution or native neighbor-search claim. Candidate CSR is
provided by the global regular triangulation and checked against ALL generators.
"""
import argparse,ctypes as ct,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull

base=Path('/gpfs/kjhan/LagForce');parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
build=base/'Codex_production_moment_geometry_build_20261018';library=build/'libproduction_geometry.so'
sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py');clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
CAP=128;D=ct.c_double;I=ct.c_int;L=ct.c_int64
class Polygon(ct.Structure):
    _fields_=[('epoch',ct.c_uint64),('count',I),('area',D),('centroid_relative',D*2),('second_relative',D*4),('x',D*CAP),('y',D*CAP),('neighbor',L*CAP),('neighbor_kind',I*CAP),('nx',D*CAP),('ny',D*CAP),('length',D*CAP),('wn_midpoint',D*CAP)]
lib=ct.CDLL(str(library));lib.gpu_moment_geometry_polygon_size.restype=ct.c_size_t
assert lib.gpu_moment_geometry_polygon_size()==ct.sizeof(Polygon)
lib.gpu_moment_geometry_export.argtypes=[I,I,I,I,ct.POINTER(L),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(I),ct.POINTER(L),ct.POINTER(L),ct.POINTER(L),ct.POINTER(D),ct.c_uint64,ct.POINTER(Polygon),ct.POINTER(I)]
lib.gpu_moment_geometry_export.restype=I

def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def call(data,backend,threads,expect=0):
    arrays=[np.ascontiguousarray(data[k],dtype=dtype) for k,dtype in [('ids',np.int64),('xy',float),('W',float),('v',float),('z',float),('kind',np.int32),('selected',np.int64),('offsets',np.int64),('candidates',np.int64),('bounds',float)]]
    out=(Polygon*len(arrays[6]))();ct.memset(out,0xA5,ct.sizeof(out));before=bytes(out);status=np.full(len(out),-999,dtype=np.int32)
    ptrs=[a.ctypes.data_as(ct.POINTER(t)) for a,t in zip(arrays,[L,D,D,D,D,I,L,L,L,D])]
    code=lib.gpu_moment_geometry_export(backend,threads,len(out),len(arrays[0]),*ptrs,int(data['epoch']),out,status.ctypes.data_as(ct.POINTER(I)))
    assert code==expect,(code,expect,status.tolist())
    if expect:assert bytes(out)==before;return dict(status=code,caller_unchanged=True)
    assert np.all(status==0)
    result=dict(count=np.array([a.count for a in out]),epoch=np.array([a.epoch for a in out],dtype=np.uint64),area=np.array([a.area for a in out]),centroid=np.array([a.centroid_relative for a in out]),second=np.array([a.second_relative for a in out]).reshape(-1,2,2))
    for key in ('x','y','neighbor','neighbor_kind','nx','ny','length','wn_midpoint'):result[key]=np.array([getattr(a,key) for a in out])
    # Unused entries have no geometry meaning; zero them before snapshots.
    for i,n in enumerate(result['count']):
        for key in ('x','y','neighbor','neighbor_kind','nx','ny','length','wn_midpoint'):result[key][i,n:]=0
    assert all(np.isfinite(a).all() for a in result.values()) and np.all(result['epoch']==data['epoch'])
    return result

def moments(p):
    # Independent origin-triangle integration of an oriented polygon.
    b=np.roll(p,-1,axis=0);a=.5*np.cross(p,b);area=a.sum()
    q=np.sum(a[:,None]*(p+b)/3,axis=0)
    Q=np.einsum('n,ni,nj->ij',a/12,p,p)+np.einsum('n,ni,nj->ij',a/12,b,b)+np.einsum('n,ni,nj->ij',a/12,p+b,p+b)
    return area,q/area,Q

def candidates(x,W,selected):
    hull=ConvexHull(np.column_stack([x,np.sum(x*x,axis=1)-W]));tri=hull.simplices[hull.equations[:,2]<-1e-12]
    pairs=np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]]]),axis=1),axis=0)
    nb=[[] for _ in x]
    for i,j in pairs:nb[i].append(j);nb[j].append(i)
    rows=[np.array(sorted(nb[i]),dtype=np.int64) for i in selected]
    assert all(len(row) for row in rows)
    return np.r_[0,np.cumsum([len(row) for row in rows])].astype(np.int64),np.concatenate(rows)

def audit(data,out,polygons):
    eps=np.finfo(float).eps;selected=data['selected'];xy=data['xy'];W=data['W'];n=len(selected)
    # Forward coordinate uncertainty from clipping planes and input-coordinate
    # subtraction, proportional to arithmetic scale and actual plane count.
    scale=np.max(abs(xy))+np.max(abs(data['bounds']))
    tolerance=64*eps*scale*(np.diff(data['offsets'])+1)
    stats=dict(maximum_coordinate_bound=float(tolerance.max()),maximum_vertex_difference=0.,maximum_area_difference=0.,maximum_centroid_difference=0.,maximum_second_difference=0.,maximum_power_plane_bound_fraction=0.,maximum_all_generator_bound_fraction=0.,maximum_normal_sum_bound_fraction=0.,maximum_shared_endpoint_difference=0.,maximum_wn_difference=0.,maximum_tangent_normal_bound_fraction=0.,external_faces=0,uncertain_short_faces=0,shared_faces=0)
    lookup={int(i):k for k,i in enumerate(selected)};face={};vertices=[];owners=[];vtols=[];hydro_bc=[];owner_ids=[];neighbor_ids=[];roles=[];geometric=[]
    for i,gid in enumerate(selected):
        nc=int(out['count'][i]);p=np.column_stack([out['x'][i,:nc],out['y'][i,:nc]]);reference=polygons[i]-xy[gid];tol=tolerance[i]
        diff=np.linalg.norm(p[:,None]-reference[None,:],axis=2);haus=max(diff.min(0).max(),diff.min(1).max());assert haus<=tol,('vertices',i,haus,tol)
        stats['maximum_vertex_difference']=max(stats['maximum_vertex_difference'],float(haus))
        area,c,Q=moments(reference);rmax=np.linalg.norm(reference,axis=1).max();perimeter=np.linalg.norm(np.roll(reference,-1,axis=0)-reference,axis=1).sum();da=2*perimeter*tol+nc*tol*tol
        dc=4*rmax*da/area+2*tol;dQ=4*rmax*rmax*da+4*area*rmax*tol
        assert abs(area-out['area'][i])<=da and np.linalg.norm(c-out['centroid'][i])<=dc and np.linalg.norm(Q-out['second'][i])<=dQ
        stats['maximum_area_difference']=max(stats['maximum_area_difference'],float(abs(area-out['area'][i])))
        stats['maximum_centroid_difference']=max(stats['maximum_centroid_difference'],float(np.linalg.norm(c-out['centroid'][i])))
        stats['maximum_second_difference']=max(stats['maximum_second_difference'],float(np.linalg.norm(Q-out['second'][i])))
        closure=np.linalg.norm(np.column_stack([out['nx'][i,:nc],out['ny'][i,:nc]]).T@out['length'][i,:nc]);bound=8*nc*tol;assert closure<=bound
        stats['maximum_normal_sum_bound_fraction']=max(stats['maximum_normal_sum_bound_fraction'],float(closure/bound))
        for e in range(nc):
            j=int(out['neighbor'][i,e]);assert 0<=j<len(xy) and out['neighbor_kind'][i,e]==data['kind'][j]
            a=p[e];b=p[(e+1)%nc];delta=xy[j]-xy[gid];d=np.linalg.norm(delta);normal=delta/d
            assert np.linalg.norm(normal-[out['nx'][i,e],out['ny'][i,e]])<=16*eps
            plane=abs(2*np.array([a,b])@delta-(delta@delta+W[gid]-W[j]))/(2*d)
            assert plane.max()<=tol,('plane',i,j,plane,tol)
            tangent_residual=abs((b-a)@normal);assert tangent_residual<=2*tol
            stats['maximum_tangent_normal_bound_fraction']=max(stats['maximum_tangent_normal_bound_fraction'],float(tangent_residual/(2*tol)))
            stats['maximum_power_plane_bound_fraction']=max(stats['maximum_power_plane_bound_fraction'],float(plane.max()/tol))
            point=(a+b)/2;wn=(delta@data['v'][j]-(data['v'][j]-data['v'][gid])@point)/d+(data['z'][gid]-data['z'][j])/(2*d)
            terms=(abs(delta*data['v'][j]).sum()+abs((data['v'][j]-data['v'][gid])*point).sum()+.5*(abs(data['z'][gid])+abs(data['z'][j])))/d
            assert abs(wn-out['wn_midpoint'][i,e])<=64*eps*terms
            stats['maximum_wn_difference']=max(stats['maximum_wn_difference'],float(abs(wn-out['wn_midpoint'][i,e])))
            external=j not in lookup
            if external:
                assert data['kind'][gid]==1 and data['kind'][j]==2
                assert np.all(data['v'][gid]==0) and np.all(data['v'][j]==0) and data['z'][gid]==0 and data['z'][j]==0
                assert wn==0 and out['wn_midpoint'][i,e]==0
                stats['external_faces']+=1
            hydro_bc.append(0. if external else wn);geometric.append(out['wn_midpoint'][i,e]);owner_ids.append(gid);neighbor_ids.append(j);roles.append((data['kind'][gid],data['kind'][j]))
            key=tuple(sorted((int(gid),j)))+(0,0)
            face.setdefault(key,[]).append((int(gid),p[e]+xy[gid],p[(e+1)%nc]+xy[gid],normal,float(out['length'][i,e]),tol))
        vertices.extend(p);owners.extend([gid]*nc);vtols.extend([tol]*nc)
    for key,rows in face.items():
        if key[0] not in lookup or key[1] not in lookup:assert len(rows)==1;continue
        tol=sum(r[5] for r in rows)
        if len(rows)==1:
            assert rows[0][4]<=2*tol,('missing finite shared face',key,rows)
            stats['uncertain_short_faces']+=1;continue
        assert len(rows)==2 and rows[0][0]!=rows[1][0],('duplicate face',key)
        a,b=rows;diff=max(np.linalg.norm(a[1]-b[2]),np.linalg.norm(a[2]-b[1]));assert diff<=tol,('shared',key,diff,tol)
        assert np.linalg.norm(a[3]+b[3])<=32*eps
        stats['maximum_shared_endpoint_difference']=max(stats['maximum_shared_endpoint_difference'],float(diff));stats['shared_faces']+=1
    # Independent tiled all-registry halfspace check. Never trusts the CSR.
    vertices=np.asarray(vertices);owners=np.asarray(owners);vtols=np.asarray(vtols)
    for begin in range(0,len(vertices),128):
        p=vertices[begin:begin+128];oi=owners[begin:begin+128];tol=vtols[begin:begin+128]
        for start in range(0,len(xy),2048):
            stop=min(start+2048,len(xy));delta=xy[start:stop][None,:,:]-xy[oi,None,:];d2=np.sum(delta*delta,axis=2);distance=np.sqrt(d2)
            value=2*np.einsum('mni,mi->mn',delta,p)-d2-W[oi,None]+W[None,start:stop]
            arithmetic=64*eps*(2*np.sum(abs(delta*p[:,None,:]),axis=2)+d2+abs(W[oi,None])+abs(W[None,start:stop]))
            bound=2*distance*tol[:,None]+arithmetic;value=np.maximum(value,0.)
            assert np.all(value<=bound),('missing candidate',begin,start,float(np.max(value-bound)))
            ratio=np.divide(value,bound,out=np.zeros_like(value),where=bound>0)
            stats['maximum_all_generator_bound_fraction']=max(stats['maximum_all_generator_bound_fraction'],float(ratio.max()))
    return stats,dict(owner_global_ids=np.array(owner_ids),neighbor_global_ids=np.array(neighbor_ids),roles=np.array(roles),geometry_wn=np.array(geometric),hydro_bc_wn=np.array(hydro_bc))

start=time.monotonic();ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
initialxy=np.column_stack([ic['x']-3,ic['y']-3]);reports={}
try:
    for stage,label in enumerate(('stage0','midpoint')):
        saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz');a=saved['state'][:-4].reshape(-1,20);selected=saved['ids'].astype(np.int64);n=len(selected);xy=initialxy.copy();xy[selected]=a[:,:2];W=np.zeros(len(xy));W[selected]=a[:,2];v=np.zeros_like(xy);v[selected]=saved['generator'];z=np.zeros(len(xy));z[selected]=saved['weight_rates'];kind=np.full(len(xy),2,dtype=np.int32);kind[selected]=saved['frozen'].astype(np.int32)
        assert np.all(v[selected[kind[selected]==1]]==0) and np.all(z[selected[kind[selected]==1]]==0)
        off,nb=candidates(xy,W,selected);data=dict(ids=np.arange(len(xy),dtype=np.int64),xy=xy,W=W,v=v,z=z,kind=kind,selected=selected,offsets=off,candidates=nb,bounds=np.array([-3.,3.,-3.,3.]),epoch=stage+1)
        np.savez(args.output/(label+'_inputs.npz'),**data)
        builder=EmbeddedLaguerreGeometry(initialxy,selected,clip.clip);polygons,_=builder(a[:,:2],a[:,2],float(saved['time']))
        results={};checks={}
        for name,backend,threads in [('cpu1',0,1),('omp4',0,4),('cuda',1,1)]:
            result=call(data,backend,threads);results[name]=result;check,faces=audit(data,result,polygons);checks[name]=check
            np.savez(args.output/(label+'_'+name+'.npz'),**result);np.savez(args.output/(label+'_'+name+'_faces.npz'),**faces)
        assert all(np.array_equal(results['cpu1'][k],results['omp4'][k]) for k in results['cpu1'])
        # Gauge fixture changes only the represented common W offset. Binary
        # rounding in the shifted inputs is included in the clipping bound.
        shifted=dict(data,W=W-.125)
        gauge=call(shifted,1,1);gaugecheck,_=audit(shifted,gauge,polygons)
        np.savez(args.output/(label+'_gauge.npz'),**gauge)
        reports[label]=dict(time=float(saved['time']),input_sha={key:sha(value) for key,value in data.items()},cells=n,candidates=len(nb),signed_weight_min=float(W.min()),signed_weight_max=float(W.max()),backends=checks,gauge=gaugecheck,cpu_omp_bitwise_equal=True,registry='Stable IDs are original IC row indices; registry outside patch is fixed with v=Wdot=0.',candidate_source='Global regular triangulation; independently checked against all registry power halfspaces; not native neighbor search.')
        print(json.dumps(dict(stage=label,status='passed',seconds=time.monotonic()-start,report=reports[label])),flush=True)
    # Capacity and coincident input transactional failures, same API/backends.
    angles=np.arange(200)*2*np.pi/200;xy=np.vstack([np.zeros(2),np.column_stack([np.cos(angles),np.sin(angles)])]);nr=len(xy)
    fixture=dict(ids=np.arange(nr,dtype=np.int64),xy=xy,W=np.zeros(nr),v=np.zeros_like(xy),z=np.zeros(nr),kind=np.zeros(nr,dtype=np.int32),selected=np.array([0],dtype=np.int64),offsets=np.array([0,nr-1],dtype=np.int64),candidates=np.arange(1,nr,dtype=np.int64),bounds=np.array([-2.,2.,-2.,2.]),epoch=9)
    fixtures={}
    for name,backend,threads in [('cpu1',0,1),('omp4',0,4),('cuda',1,1)]:
        fixtures[name]=dict(capacity=call(fixture,backend,threads,expect=2))
        coincident=dict(fixture,xy=xy.copy());coincident['xy'][1]=0
        fixtures[name]['coincident']=call(coincident,backend,threads,expect=6)
    report=dict(status='completed',scope=__doc__,seconds=time.monotonic()-start,library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in build.iterdir() if p.suffix in ('.h','.cu','.cuh')},stages=reports,fixtures=fixtures,limits='Stored geometry only. No native candidate search, massroot, RHS, MPI, evolution or long-term physical accuracy claim. Actual production translation unit compiled at O2 without fast-math for this strict API.')
    (args.output/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
except BaseException:
    (args.output/'partial_report.json').write_text(json.dumps(reports,indent=2));raise
