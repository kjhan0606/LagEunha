"""Saved-predicate replay only: diagnose CUDA nextafter overload and filter fix.
No geometry resolution, RHS, massroot or evolution is re-executed.
"""
import argparse,ctypes as ct,hashlib,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');build=base/'Codex_geometry_v2_filter_replay_build_20261018';source=base/'Codex_full_native_geometry_v2_20261018';sys.path.insert(0,str(build))
from moment_geometry_native_v2 import Backend,D,L,I,pointer
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
backend=Backend(build/'libgeometry_v2_fixed.so');started=time.monotonic();stages={}
def emit(row):
    print(json.dumps(row),flush=True)
    with (args.output/'events.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def probe(data,queries,library,kind):
    lib=ct.CDLL(str(library));lib.moment_geometry_interval_probe.argtypes=[I,I,ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),I,ct.POINTER(L),ct.POINTER(D)];arrays=[np.ascontiguousarray(data[k],dtype=float) for k in ('xy','W','bounds')];queries=np.ascontiguousarray(queries,dtype=np.int64);assert len(queries)<=128;out=np.empty((len(queries),37))
    code=lib.moment_geometry_interval_probe(kind,len(data['xy']),*[pointer(a,D) for a in arrays],len(queries),pointer(queries,L),pointer(out,D));assert code==0;return out

data=dict(np.load(source/'stage0/inputs.npz'));record=np.load(source/'stage0/clipping_signs.npz');queries=record['queries'];take=np.flatnonzero(record['cpu']!=0)[:8];q=queries[take];probes={}
for variant in ('original','fixed'):
    for name,kind in [('cpu',0),('gpu',1)]:probes[variant+'_'+name]=probe(data,q,build/('libprobe_'+variant+'.so'),kind)
np.savez(args.output/'intermediate_probe.npz',queries=q,fields=np.array(['legacy_next_up','legacy_next_down','explicit_next_up','explicit_next_down','add_lo','add_hi','mul_lo','mul_hi']+[f'plane_{i}_{k}' for i in range(3) for k in ('nx_lo','nx_hi','ny_lo','ny_hi','b_lo','b_hi')]+['D_lo','D_hi','X_lo','X_hi','Y_lo','Y_hi','N_lo','N_hi','sign_D','sign_N','sign_result']),**probes)
expected=np.array([np.nextafter(1.,np.inf),np.nextafter(1.,-np.inf)]).view(np.uint64)
for name in ('fixed_cpu','fixed_gpu'):
    out=probes[name];assert np.all(out[:,2:4].view(np.uint64)==expected)
    assert np.all(out[:,4]<2.) and np.all(out[:,5]>2.) and np.all(out[:,6]<1.) and np.all(out[:,7]>1.)
    assert np.array_equal(out[:,36].astype(int),record['exact'][take])
probe_report={name:dict(first_nextafter_bits=value[0,:4].view(np.uint64).tolist(),first_simple_interval_hex=[x.hex() for x in value[0,4:8]],first_planes_hex=[x.hex() for x in value[0,8:26]],first_D_N_hex=[value[0,k].hex() for k in (26,27,32,33)],resolved_count=int((value[:,36]!=0).sum())) for name,value in probes.items()}
(args.output/'probe_report.json').write_text(json.dumps(probe_report,indent=2));emit(dict(probe=probe_report))
try:
    original_summary=json.loads((source/'summary.json').read_text())
    for label in ('stage0','midpoint'):
        folder=args.output/label;folder.mkdir();data=dict(np.load(source/label/'inputs.npz'));old=np.load(source/label/'clipping_signs.npz');queries=old['queries'];expected=old['exact'];cpu=backend.filter(data,queries,0,4);gpu=backend.filter(data,queries,1,1)
        assert np.array_equal(cpu,old['cpu']) and np.all((cpu==0)|(cpu==expected)) and np.all((gpu==0)|(gpu==expected))
        assert (gpu!=0).any();np.savez(folder/'clipping_signs.npz',queries=queries,exact=expected,cpu=cpu,gpu=gpu)
        clipping=dict(count=len(queries),input_sha256=sha(queries),exact_sha256=sha(expected),cpu_sha256=sha(cpu),gpu_sha256=sha(gpu),cpu_resolved=int((cpu!=0).sum()),gpu_resolved=int((gpu!=0).sum()),old_gpu_resolved=int((old['gpu']!=0).sum()),cpu_resolved_gpu_unresolved=int(((cpu!=0)&(gpu==0)).sum()),gpu_resolved_cpu_unresolved=int(((gpu!=0)&(cpu==0)).sum()),wrong_nonzero_signs=0)
        emit(dict(stage=label,clipping_filter=clipping))
        saved=np.load(source/label/'feasibility_cartesian_inputs.npz');owners=saved['owners'];supports=saved['supports'];enc=saved['enclosures'];cpu_hash=hashlib.sha256();gpu_hash=hashlib.sha256();count=0;cpu_unknown=0;gpu_unknown=0;extra=0;coverage=[];witness=[]
        assert all(np.array_equal(saved[k],data[k]) for k in ('xy','W','bounds')) and np.array_equal(saved['registry_ids'],data['ids'])
        for start in range(0,len(owners),128):
            stop=min(start+128,len(owners));sub=owners[start:stop];bounds=enc[start:stop];cpu=backend.filter(data,sub,0,4,bounds);gpu=backend.filter(data,sub,1,1,bounds)
            assert np.all(cpu<=0) and np.all(gpu<=0) and np.all((cpu==0)|(gpu==0)|(cpu==gpu))
            assert np.all(cpu[np.arange(len(sub)),sub]==0) and np.all(gpu[np.arange(len(sub)),sub]==0)
            a,b=np.where((cpu==0)&(gpu!=0))
            if len(a):
                q=np.column_stack([supports[start+a],b]);exact=backend.exact(data,q);assert np.array_equal(exact,gpu[a,b]);extra+=len(a);np.savez(folder/f'extra_resolved_{start:05d}.npz',queries=q,exact=exact,gpu=gpu[a,b])
            if len(witness)<64:
                a,b=np.where(cpu!=gpu)
                for ai,bi in zip(a[:64-len(witness)],b[:64-len(witness)]):witness.append(dict(vertex=int(start+ai),registry_row=int(bi),cpu=int(cpu[ai,bi]),gpu=int(gpu[ai,bi])))
            cpu_hash.update(cpu.tobytes());gpu_hash.update(gpu.tobytes());count+=cpu.size;cpu_unknown+=int((cpu==0).sum());gpu_unknown+=int((gpu==0).sum());coverage.append([start,stop])
            if len(coverage)%16==0:emit(dict(stage=label,vertices=stop,total_vertices=len(owners),queries=count,seconds=time.monotonic()-started))
        previous=original_summary['stages'][label]['feasibility_filter'];assert count==previous['count'] and coverage==previous['coverage'] and cpu_hash.hexdigest()==previous['cpu_sign_sha256'];assert count>gpu_unknown
        feasibility=dict(count=count,vertices=len(owners),registry_count=len(data['ids']),coverage=coverage,input_sha256=previous['input_sha256'],cpu_sign_sha256=cpu_hash.hexdigest(),gpu_sign_sha256=gpu_hash.hexdigest(),cpu_unresolved=cpu_unknown,gpu_unresolved=gpu_unknown,gpu_resolved=count-gpu_unknown,old_gpu_resolved=count-previous['gpu_unresolved'],gpu_extra_resolved_exact_checked=extra,self_exact_zero_count=len(owners),wrong_nonzero_signs=0,unresolved_difference_witness=witness)
        stages[label]=dict(clipping_filter=clipping,feasibility_filter=feasibility);(folder/'report.json').write_text(json.dumps(stages[label],indent=2));emit(dict(stage=label,status='completed',gpu_resolved=feasibility['gpu_resolved']))
    fixtures={};q=np.array([[0,1,2,3],[0,1,1,3],[0,1,3,1]],dtype=np.int64);xy=np.array([[0.,0.],[1.,0.],[1.,1e-12],[0.,1.]])
    for label,scale in [('nearly_parallel',1.),('subnormal',1e-310),('overflow_intermediates',1e200)]:
        data=dict(xy=xy*scale,W=np.zeros(4),bounds=np.array([-3.,3.,-3.,3.]));cpu=backend.filter(data,q,0,1);gpu=backend.filter(data,q,1,1);exact=backend.exact(data,q)
        assert np.all((cpu==0)|(cpu==exact)) and np.all((gpu==0)|(gpu==exact)) and np.all(gpu[exact==2]==0) and np.all(cpu[exact==2]==0)
        if label!='nearly_parallel':assert np.all(cpu==0) and np.all(gpu==0)
        fixtures[label]=dict(cpu=cpu.tolist(),gpu=gpu.tolist(),exact=exact.tolist())
    summary=dict(status='completed',scope=__doc__,seconds=time.monotonic()-started,library_sha256=hashlib.sha256((build/'libgeometry_v2_fixed.so').read_bytes()).hexdigest(),original_library_sha256=original_summary['library_sha256'],probe=probe_report,stages=stages,fixtures=fixtures,limits='Saved sign-filter replay; no repeat of exact geometry, no GPU topology-builder or provider/MPI/evolution claim.')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2));emit(dict(status='completed',seconds=summary['seconds']))
except BaseException:
    (args.output/'partial_report.json').write_text(json.dumps(stages,indent=2));raise
