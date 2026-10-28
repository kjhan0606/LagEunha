"""Two stored1590-cell snapshots: HOST exact completion and GPU sign filters.

No GPU topology-builder, provider/RHS, new massroot, MPI or evolution claim.
"""
import argparse,hashlib,importlib.util,json,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull
base=Path('/gpfs/kjhan/LagForce');build=base/'Codex_production_geometry_v2_full_build_20261018'
sys.path[:0]=[str(build),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from moment_geometry_native_v2 import Backend
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py');clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
backend=Backend(build/'libgeometry_v2.so');started=time.monotonic()
def emit(row):
    print(json.dumps(row),flush=True)
    with (args.output/'events.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def enclosure(q):
    x=float(q);f=F(x)
    return (x,x) if f==q else ((x,np.nextafter(x,np.inf)) if f<q else (np.nextafter(x,-np.inf),x))
def candidate_csr(x,W,ids):
    h=ConvexHull(np.column_stack([x,np.sum(x*x,axis=1)-W]));t=h.simplices[h.equations[:,2]<-1e-12];pairs=np.unique(np.sort(np.concatenate([t[:,[0,1]],t[:,[0,2]],t[:,[1,2]]]),axis=1),axis=0);rows=[[] for _ in x]
    for i,j in pairs:rows[i].append(j);rows[j].append(i)
    rows=[np.array(sorted(rows[i]),dtype=np.int64) for i in ids]
    return np.r_[0,np.cumsum([len(r) for r in rows])].astype(np.int64),np.concatenate(rows)
def geometric_moments(p):
    b=np.roll(p,-1,axis=0);a=.5*np.cross(p,b);A=a.sum();c=np.sum(a[:,None]*(p+b)/3,axis=0)/A
    Q=np.einsum('n,ni,nj->ij',a/12,p,p)+np.einsum('n,ni,nj->ij',a/12,b,b)+np.einsum('n,ni,nj->ij',a/12,p+b,p+b)
    return A,c,Q

def clipping_filter(data,exact,folder):
    rows=np.array([q for c in exact['cells'] for q in c['clipping_predicate_queries']],dtype=np.int64);queries=np.ascontiguousarray(rows[:,:4]);oracle=rows[:,4];host_recorded=rows[:,5]
    cpu=backend.filter(data,queries,0,4);gpu=backend.filter(data,queries,1,1)
    assert np.array_equal(cpu,host_recorded)
    assert np.all((cpu==0)|(cpu==oracle)) and np.all((gpu==0)|(gpu==oracle))
    # The host exact oracle was computed at each actual clipping call, including
    # failed candidate passes before completion. Cross-check its separate ABI.
    independent=backend.exact(data,queries);assert np.array_equal(independent,oracle)
    np.savez(folder/'clipping_signs.npz',queries=queries,exact=oracle,cpu=cpu,gpu=gpu)
    return dict(count=len(rows),input_sha256=sha(queries),oracle_sha256=sha(oracle),cpu_sha256=sha(cpu),gpu_sha256=sha(gpu),cpu_unresolved=int((cpu==0).sum()),gpu_unresolved=int((gpu==0).sum()),wrong_nonzero_signs=0)

def feasibility_filter(data,exact,folder):
    idx={int(i):j for j,i in enumerate(data['ids'])};owners=[];supports=[];enc=[]
    for cell in exact['cells']:
        i=idx[cell['owner']];xi,yi=map(F,data['xy'][i])
        for edge in cell['edges']:
            pa=list(map(F,edge['A_global']));x=pa[0]-xi;y=pa[1]-yi;bounds=(*enclosure(x),*enclosure(y))
            assert F(bounds[0])<=x<=F(bounds[1]) and F(bounds[2])<=y<=F(bounds[3])
            a,b=edge['support'];supports.append([i,idx[a] if a>=0 else a,idx[b] if b>=0 else b]);owners.append(i);enc.append(bounds)
    owners=np.array(owners,dtype=np.int64);supports=np.array(supports,dtype=np.int64);enc=np.array(enc)
    np.savez(folder/'feasibility_cartesian_inputs.npz',owners=owners,supports=supports,enclosures=enc,registry_ids=data['ids'],xy=data['xy'],W=data['W'],bounds=data['bounds'])
    count=0;self_count=0;cpu_unknown=0;gpu_unknown=0;extra_exact=0;cpu_hash=hashlib.sha256();gpu_hash=hashlib.sha256();coverage=[];witness=[]
    for start in range(0,len(owners),128):
        stop=min(start+128,len(owners));sub=owners[start:stop];bounds=enc[start:stop]
        cpu=backend.filter(data,sub,0,4,bounds);gpu=backend.filter(data,sub,1,1,bounds)
        assert np.all(cpu<=0) and np.all(gpu<=0),'completed exact polygon cannot be outside'
        assert np.all((cpu==0)|(gpu==0)|(cpu==gpu))
        assert np.all(cpu[np.arange(len(sub)),sub]==0) and np.all(gpu[np.arange(len(sub)),sub]==0);self_count+=len(sub)
        a,b=np.where((cpu==0)&(gpu!=0))
        if len(a):
            q=np.column_stack([supports[start+a],b]);expected=backend.exact(data,q);assert np.all(gpu[a,b]==expected);extra_exact+=len(a)
            np.savez(folder/f'gpu_extra_resolved_{start:05d}.npz',queries=q,exact=expected,gpu=gpu[a,b])
        if len(witness)<64:
            a,b=np.where((cpu==0)|(gpu==0));mask=b!=sub[a]
            for ai,bi in zip(a[mask][:64-len(witness)],b[mask][:64-len(witness)]):witness.append(dict(vertex=int(start+ai),owner=int(sub[ai]),tested_registry_row=int(bi),cpu=int(cpu[ai,bi]),gpu=int(gpu[ai,bi])))
        count+=cpu.size;cpu_unknown+=int((cpu==0).sum());gpu_unknown+=int((gpu==0).sum());cpu_hash.update(cpu.tobytes());gpu_hash.update(gpu.tobytes());coverage.append([start,stop])
        if len(coverage)%16==0:emit(dict(stage=folder.name,feasibility_vertices=stop,total_vertices=len(owners),queries=count,seconds=time.monotonic()-started))
    assert count==len(owners)*len(data['ids']) and self_count==len(owners)
    report=dict(vertices=len(owners),registry_count=len(data['ids']),count=count,self_exact_zero_count=self_count,coverage=coverage,input_sha256={k:sha(a) for k,a in [('owners',owners),('supports',supports),('enclosures',enc),('registry_ids',data['ids']),('xy',data['xy']),('W',data['W']),('bounds',data['bounds'])]},cpu_sign_sha256=cpu_hash.hexdigest(),gpu_sign_sha256=gpu_hash.hexdigest(),cpu_unresolved=cpu_unknown,gpu_unresolved=gpu_unknown,gpu_extra_resolved_exact_checked=extra_exact,wrong_nonzero_signs=0,unresolved_witness=witness,certificate='Exact rational vertices were outward enclosed; all Cartesian CPU/GPU signs compared; CPU resolved signs are outward interval certificates, GPU extra resolved signs have native exact support oracle checks. Self plane is exact zero.')
    (folder/'feasibility_report.json').write_text(json.dumps(report,indent=2));return report

def export_report(data,out,exact,folder,polygons):
    selected=set(map(int,data['selected']));idrow={int(g):i for i,g in enumerate(data['ids'])};canonical=[];keys=[];added=[];lengths=[];reference=[];external=0
    for i,cell in enumerate(exact['cells']):
        gid=cell['owner'];owner=idrow[gid];n=int(out['count'][i]);assert n==len(cell['edges']) and out['area'][i]>0
        legacy=geometric_moments(polygons[i]-data['xy'][owner]);reference.append([abs(out['area'][i]-legacy[0]),np.linalg.norm(out['centroid'][i]-legacy[1]),np.linalg.norm(out['second'][i].reshape(2,2)-legacy[2])])
        added.extend(dict(owner=gid,**w) for w in cell['witnesses'])
        for e,edge in enumerate(cell['edges']):
            j=idrow[edge['neighbor']];q=F(edge['length_squared']);lo=out['length_lo'][i,e];hi=out['length_hi'][i,e];assert q>0 and F(lo)**2<=q<=F(hi)**2 and hi<=np.nextafter(lo,np.inf)
            A=list(map(F,edge['A_global']));B=list(map(F,edge['B_global']));ex=B[0]-A[0];ey=B[1]-A[1]
            assert out['edge_dx'][i,e]==float(ex) and out['edge_dy'][i,e]==float(ey)
            assert out['midpoint_x'][i,e]==float((A[0]+B[0])/2-F(data['xy'][owner,0])) and out['midpoint_y'][i,e]==float((A[1]+B[1])/2-F(data['xy'][owner,1]))
            lengths.append(out['length'][i,e]);is_external=j not in selected
            if is_external:
                assert data['kind'][owner]==1 and data['kind'][j]==2 and np.all(data['v'][[owner,j]]==0) and np.all(data['z'][[owner,j]]==0) and out['geometry_wn'][i,e]==0
                external+=1
            if not is_external and gid>edge['neighbor']:continue
            key=(min(gid,edge['neighbor']),max(gid,edge['neighbor']),0,0);assert key not in keys;keys.append(key)
            canonical.append([gid,edge['neighbor'],float((A[0]+B[0])/2),float((A[1]+B[1])/2),float(ex),float(ey),out['length'][i,e],out['nx'][i,e],out['ny'][i,e],out['geometry_wn'][i,e],0. if is_external else out['geometry_wn'][i,e],data['kind'][owner],data['kind'][j]])
    np.savez(folder/'canonical_owner_faces.npz',keys=np.array(keys,dtype=np.int64),data=np.array(canonical),columns=np.array(['owner_global_id','neighbor_global_id','midpoint_global_x','midpoint_global_y','edge_dx','edge_dy','length','nx','ny','geometry_wn','hydro_bc_wn','owner_kind','neighbor_kind']))
    stats=dict(positive_oriented_faces=len(lengths),canonical_faces=len(keys),external_fixed_wall_faces=external,length_quantiles=np.quantile(lengths,[0,.01,.5,.99,1]).tolist(),faces_below_1e15=int((np.array(lengths)<1e-15).sum()),added_candidates=added,endpoint_rounding_difference=int(out['endpoint_rounding_difference'].sum()),endpoint_collapse=int(out['endpoint_collapse'].sum()),endpoint_reversal=int(out['endpoint_reversal'].sum()),maximum_relative_endpoint_vector_error=float(out['max_relative_vector_error'].max()),legacy_numeric_reference_maximum_difference=np.max(reference,axis=0).tolist(),legacy_reference_scope='Area/relative centroid/second moment numeric differences only. Exact micro-topology is not forced to match Qhull.',shared_global_endpoints_exact_reversed_equal=exact['shared_global_endpoints_exact_reversed_equal'],native_filtered_checks=int(out['filtered_checks'].sum()),native_exact_resolution_checks=int(out['exact_checks'].sum()))
    return stats

reports={}
try:
    initial=dict(np.load(base/'Codex_production_moment_geometry_20261018/stage0_inputs.npz'));initialxy=initial['xy'].copy()
    for epoch,label in enumerate(('stage0','midpoint'),1):
        folder=args.output/label;folder.mkdir();saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz');a=saved['state'][:-4].reshape(-1,20);ids=saved['ids'];xy=initialxy.copy();xy[ids]=a[:,:2];W=np.zeros(len(xy));W[ids]=a[:,2];v=np.zeros_like(xy);v[ids]=saved['generator'];z=np.zeros(len(xy));z[ids]=saved['weight_rates'];kind=np.full(len(xy),2,dtype=np.int32);kind[ids]=saved['frozen'].astype(np.int32)
        off,nb=(initial['offsets'],initial['candidates']) if label=='stage0' else candidate_csr(xy,W,ids)
        data=dict(ids=np.arange(len(xy),dtype=np.int64),xy=xy,W=W,v=v,z=z,kind=kind,selected=ids,offsets=off,candidates=nb,bounds=initial['bounds'],epoch=epoch);np.savez(folder/'inputs.npz',**data)
        outputs={};times={};exact=None
        for name,threads in [('cpu1',1),('omp4',4)]:
            emit(dict(stage=label,operation='explicit_host_exact_resolve',backend=name,status='starting'))
            outputs[name],audit,times[name]=backend.resolve(data,folder/name,threads)
            if exact is None:exact=audit
            else:assert audit==exact
            emit(dict(stage=label,operation='explicit_host_exact_resolve',backend=name,status='completed',seconds=times[name]))
        assert all(np.array_equal(outputs['cpu1'][k],outputs['omp4'][k]) for k in outputs['cpu1'])
        polys,_=EmbeddedLaguerreGeometry(initialxy,ids,clip.clip)(a[:,:2],a[:,2],float(saved['time']))
        exports=export_report(data,outputs['cpu1'],exact,folder,polys)
        clipping=clipping_filter(data,exact,folder);emit(dict(stage=label,clipping_filter=clipping))
        feasibility=feasibility_filter(data,exact,folder)
        reports[label]=dict(time=float(saved['time']),host_seconds=times,cpu_omp_bitwise_equal=True,export=exports,clipping_filter=clipping,feasibility_filter=feasibility)
        (folder/'report.json').write_text(json.dumps(reports[label],indent=2));emit(dict(stage=label,status='completed',export=exports))
    # Filter singular, almost parallel, subnormal and overflowing intermediates.
    fixtures={};queries=np.array([[0,1,2,3],[0,1,1,3],[0,1,3,1]],dtype=np.int64)
    positions=np.array([[0.,0.],[1.,0.],[1.,1e-12],[0.,1.]])
    for name,scale in [('nearly_parallel',1.),('subnormal',1e-310),('overflow_intermediates',1e200)]:
        data=dict(xy=positions*scale,W=np.zeros(4),bounds=np.array([-3.,3.,-3.,3.]));cpu=backend.filter(data,queries,0,1);gpu=backend.filter(data,queries,1,1);q=backend.exact(data,queries)
        assert np.all((cpu==0)|(cpu==q)) and np.all((gpu==0)|(gpu==q));assert np.all(cpu[q==2]==0) and np.all(gpu[q==2]==0)
        if name in ('subnormal','overflow_intermediates'):assert np.all(cpu==0) and np.all(gpu==0)
        fixtures[name]=dict(cpu=cpu.tolist(),gpu=gpu.tolist(),exact=q.tolist())
    report=dict(status='completed',scope=__doc__,seconds=time.monotonic()-started,library_sha256=hashlib.sha256((build/'libgeometry_v2.so').read_bytes()).hexdigest(),source_sha256={name:hashlib.sha256((build/name).read_bytes()).hexdigest() for name in ('moment_geometry_v2_exact.cpp','moment_geometry_v2_filter.cu','moment_geometry_interval.h','moment_geometry_v2.h','moment_geometry_native_v2.py','check_full_native_geometry_v2.py')},stages=reports,filter_fixtures=fixtures,limits='Host exact completion and GPU interval SIGN filter only. GPU topology construction, native neighbor search, density/face provider, massroot, MPI, physical evolution and long accuracy remain incomplete. Endpoint-only provider remains unapproved.')
    (args.output/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status='completed',seconds=report['seconds'])),flush=True)
except BaseException:
    (args.output/'partial_report.json').write_text(json.dumps(reports,indent=2));raise
