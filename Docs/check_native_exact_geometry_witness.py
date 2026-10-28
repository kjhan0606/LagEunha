"""Explicit native HOST exact-resolution of saved cells 3/5; no new evolution."""
import argparse,ctypes as ct,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');build=base/'Codex_production_geometry_v2_build_20261018';parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
D=ct.c_double;L=ct.c_int64;U=ct.c_uint64;I=ct.c_int;CAP=128
class Polygon(ct.Structure):
    _fields_=[('epoch',U),('count',I),('endpoint_rounding_difference',I),('endpoint_collapse',I),('endpoint_reversal',I),('max_relative_vector_error',D),('filtered_checks',U),('exact_checks',U),('added_candidates',U),('area',D),('centroid',D*2),('second',D*4)]+[(name,typ*CAP) for name,typ in [('x',D),('y',D),('support_a',L),('support_b',L),('neighbor',L),('neighbor_kind',I),('midpoint_x',D),('midpoint_y',D),('edge_dx',D),('edge_dy',D),('length',D),('length_lo',D),('length_hi',D),('nx',D),('ny',D),('geometry_wn',D)]]
lib=ct.CDLL(str(build/'libgeometry_exact_v2.so'));lib.moment_geometry_v2_size.restype=ct.c_size_t;assert ct.sizeof(Polygon)==lib.moment_geometry_v2_size()
lib.moment_geometry_v2_resolve.argtypes=[I,I,I,ct.POINTER(L),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(D),ct.POINTER(I),ct.POINTER(L),ct.POINTER(L),ct.POINTER(L),ct.POINTER(D),U,ct.POINTER(Polygon),ct.POINTER(I),ct.c_char_p];lib.moment_geometry_v2_resolve.restype=I
original=np.load(base/'Codex_production_moment_geometry_20261018/stage0_inputs.npz');data={k:original[k].copy() for k in original.files};lookup={int(g):i for i,g in enumerate(data['selected'])};rows=[data['candidates'][data['offsets'][lookup[g]]:data['offsets'][lookup[g]+1]] for g in (3,5)];data.update(selected=np.array([3,5],dtype=np.int64),offsets=np.r_[0,np.cumsum([len(r) for r in rows])].astype(np.int64),candidates=np.concatenate(rows))
np.savez(args.output/'inputs.npz',**data)

def call(data,name,threads,expected=0):
    typed=[np.ascontiguousarray(data[k],dtype=d) for k,d in [('ids',np.int64),('xy',float),('W',float),('v',float),('z',float),('kind',np.int32),('selected',np.int64),('offsets',np.int64),('candidates',np.int64),('bounds',float)]]
    out=(Polygon*len(typed[6]))();ct.memset(out,0xA5,ct.sizeof(out));before=bytes(out);status=np.full(len(out),-999,dtype=np.int32);path=args.output/(name+'_exact.json')
    pointers=[a.ctypes.data_as(ct.POINTER(t)) for a,t in zip(typed,[L,D,D,D,D,I,L,L,L,D])];start=time.monotonic()
    code=lib.moment_geometry_v2_resolve(threads,len(out),len(typed[0]),*pointers,int(data['epoch']),out,status.ctypes.data_as(ct.POINTER(I)),str(path).encode());elapsed=time.monotonic()-start
    assert code==expected,(name,code,status.tolist())
    if expected:assert bytes(out)==before;return dict(code=code,caller_unchanged=True,seconds=elapsed)
    assert np.all(status==0)
    arrays={name:np.array([getattr(p,name) for p in out]) for name,_ in Polygon._fields_}
    for i,n in enumerate(arrays['count']):
        for field,typ in Polygon._fields_:
            if getattr(typ,'_length_',0)==CAP:arrays[field][i,n:]=0
    np.savez(args.output/(name+'_export.npz'),**arrays)
    return dict(seconds=elapsed,arrays=arrays,exact=json.loads(path.read_text()))

oracle=json.loads((base/'Codex_production_moment_geometry_20261018/exact_contact_all_registry_diagnostic.json').read_text())
def rational(pair):return F(int(pair[0]),int(pair[1]))
def check(result):
    arrays=result['arrays'];native=result['exact'];rows=[]
    assert native['shared_global_endpoints_exact_reversed_equal']
    for k,(ref,actual) in enumerate(zip(oracle['cells'],native['cells'])):
        gid=ref['owner'];assert actual['owner']==gid;expected={e['neighbor']:e for e in ref['edges']};assert set(expected)=={e['neighbor'] for e in actual['edges']}
        xi,yi=map(F,data['xy'][gid]);points=[]
        for e in actual['edges']:
            old=expected[e['neighbor']];A=list(map(F,e['A_global']));B=list(map(F,e['B_global']));assert A==list(map(rational,old['A_global_exact'])) and B==list(map(rational,old['B_global_exact']))
            assert F(e['length_squared'])==rational(old['length_squared_exact'])>0;points.append((A[0]-xi,A[1]-yi))
        twice=F(0);qx=F(0);qy=F(0);xx=F(0);yy=F(0);xy=F(0)
        for e,(ax,ay) in enumerate(points):
            bx,by=points[(e+1)%len(points)];cross=ax*by-bx*ay;twice+=cross;qx+=(ax+bx)*cross;qy+=(ay+by)*cross;xx+=(ax*ax+ax*bx+bx*bx)*cross;yy+=(ay*ay+ay*by+by*by)*cross;xy+=(2*ax*ay+ax*by+bx*ay+2*bx*by)*cross
        assert arrays['area'][k]==float(twice/2)
        assert np.array_equal(arrays['centroid'][k],[float(qx/(3*twice)),float(qy/(3*twice))])
        assert np.array_equal(arrays['second'][k],[float(xx/12),float(xy/24),float(xy/24),float(yy/12)])
        for e,edge in enumerate(actual['edges']):
            A=list(map(F,edge['A_global']));B=list(map(F,edge['B_global']));ex=B[0]-A[0];ey=B[1]-A[1];q=F(edge['length_squared']);lo=arrays['length_lo'][k,e];hi=arrays['length_hi'][k,e]
            assert F(lo)**2<=q<=F(hi)**2 and hi<=np.nextafter(lo,np.inf)
            assert arrays['edge_dx'][k,e]==float(ex) and arrays['edge_dy'][k,e]==float(ey)
            assert arrays['midpoint_x'][k,e]==float((A[0]+B[0])/2-xi) and arrays['midpoint_y'][k,e]==float((A[1]+B[1])/2-yi)
            assert arrays['x'][k,e]==float(A[0]-xi) and arrays['y'][k,e]==float(A[1]-yi)
        added=[w['added_id'] for w in actual['witnesses']];assert added==([] if gid==3 else [2])
        rows.append(dict(owner=gid,edges=len(actual['edges']),added_global_ids=added,filtered_checks=actual['filtered_checks'],exact_checks=actual['exact_checks'],endpoint_rounding_difference=int(arrays['endpoint_rounding_difference'][k]),endpoint_collapse=int(arrays['endpoint_collapse'][k]),endpoint_reversal=int(arrays['endpoint_reversal'][k]),max_relative_vector_error=float(arrays['max_relative_vector_error'][k])))
    return rows

results={name:call(data,name,threads) for name,threads in [('cpu1',1),('omp4',4)]};checks={name:check(r) for name,r in results.items()}
assert all(np.array_equal(results['cpu1']['arrays'][k],results['omp4']['arrays'][k]) for k in results['cpu1']['arrays'])
# Exact gauge: original W all zero, shift -.125 exactly preserves every Wi-Wj.
assert np.all(data['W']==0);gauge=dict(data,W=np.full_like(data['W'],-.125));gr=call(gauge,'negative_common_gauge',1);check(gr)
assert all(np.array_equal(results['cpu1']['arrays'][k],gr['arrays'][k]) for k in results['cpu1']['arrays'])
# Deterministic registry permutation with candidate order by stable IDs.
perm=np.arange(len(data['ids']))[::-1];inverse=np.argsort(perm);reordered=dict(data)
for key in ('ids','xy','W','v','z','kind'):reordered[key]=data[key][perm]
reordered['selected']=inverse[data['selected']];reordered['candidates']=inverse[data['candidates']]
rr=call(reordered,'reordered_registry',1);check(rr)
assert all(np.array_equal(results['cpu1']['arrays'][k],rr['arrays'][k]) for k in results['cpu1']['arrays'])
# Coincident generator: transactional explicit failure, no replacement geometry.
bad=dict(data,xy=data['xy'].copy());bad['xy'][2]=bad['xy'][3]
fixture=call(bad,'coincident',1,6)
report=dict(status='completed',scope=__doc__,library_sha256=hashlib.sha256((build/'libgeometry_exact_v2.so').read_bytes()).hexdigest(),checks=checks,seconds={k:v['seconds'] for k,v in results.items()},cpu_omp_bitwise_equal=True,negative_exact_gauge_bitwise_equal=True,registry_permutation_bitwise_equal=True,coincident=fixture,shared_edge_3_5_length=next(float(rational(e['length_squared_exact']))**.5 for e in oracle['cells'][0]['edges'] if e['neighbor']==5),limits='Two cells in explicit native HOST rational resolution. GPU fast path, whole1590-cell export, production provider, massroot/MPI/evolution and long accuracy remain incomplete. Endpoint-only downstream cannot consume exports without a separate representation gate.')
(args.output/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
