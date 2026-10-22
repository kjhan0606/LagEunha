from pathlib import Path
import argparse,json,struct
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--version',default='v1');args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[1];dest=ROOT/f'rk4_diagnostics/enriched_{args.version}'
manifest=json.loads((dest/'manifest.json').read_text());cases={};result=dict(scope=manifest['scope'],cases=[])
def make_dtype(nm):
 return np.dtype(dict(names=['id','x','y','geometry','generator','scale','L','mean','diagnostics','rank','basis_status','coeff','integral','metrics','flags'],
    formats=['<u8','<f8','<f8',('<f8',6),('<f8',2),('<f8',4),('<f8',(10,10)),('<f8',10),('<f8',4),'<i4','<i4',('<f8',(4,10)),('<f8',4),('<f8',nm),('<i4',4)],
    offsets=[0,8,16,24,72,88,120,920,1000,1032,1036,1040,1360,1392,1392+8*nm],itemsize=1408+8*nm))
def load(row):
    pieces=[]
    for r in range(row['ranks']):
        path=dest/row['name']/f'enriched_project_r{r}.bin'
        with path.open('rb') as f:h=struct.unpack('<8Q',f.read(64))
        assert h[0]==0x474653454e523031 and (h[1],h[3],h[4],h[5]) in ((1,10,14,1520),(2,10,15,1528))
        dtype=make_dtype(h[4]);assert path.stat().st_size==64+h[2]*dtype.itemsize
        pieces.append(np.fromfile(path,dtype=dtype,offset=64,count=h[2]))
    a=np.concatenate(pieces);a=a[np.argsort(a['id'])];assert len(np.unique(a['id']))==len(a);return a
for row in manifest['cases']:
    a=load(row);cases[row['name']]=a;m=a['metrics'];d=a['diagnostics'];ranks,counts=np.unique(a['rank'],return_counts=True)
    worst=int(np.argmax(m[:,2]));orth=int(np.argmax(d[:,2]));quad=int(np.argmax(m[:,13]))
    v=dict(name=row['name'],cells=len(a),all_build_valid=bool(np.all(a['flags'][:,0]==1)),
        rank_counts={str(k):int(n) for k,n in zip(ranks,counts)},
        nonpositive_rho_cells=int(np.count_nonzero(a['flags'][:,1])),nonpositive_P_cells=int(np.count_nonzero(a['flags'][:,2])),
        min_rho=float(m[:,0].min()),min_P=float(m[:,1].min()),max_relative_P_error=float(m[:,2].max()),
        worst_P_id=int(a[worst]['id']),worst_P_r=float(np.hypot(a[worst]['x']-12,a[worst]['y']-12)),
        volume_L1_relative_P_error=float(np.sum(a['geometry'][:,0]*m[:,3])/np.sum(a['geometry'][:,0])),
        max_relative_field_error=float(m[:,4].max()),max_boost_relative_P_error=float(m[:,5].max()),
        max_independent_conserved_integral_error=float(m[:,6:10].max()),
        max_quadrature_conserved_error=float(m[:,10:13].max()),max_quadrature_E_over_thermal=float(m[:,13].max()),
        quadrature_E_worst_id=int(a[quad]['id']),maximum_orthogonality_error=float(d[:,2].max()),
        orthogonality_worst_id=int(a[orth]['id']),min_retained_residual_norm=float(d[:,0].min()),
        max_excluded_residual_norm=float(d[:,1].max()),max_area_error=float(d[:,3].max()))
    if m.shape[1]>14:v['stable_evaluation_max_relative_rho_change']=float(m[:,14].max())
    if row['name'].startswith('fine'):
        v['focus']=[]
        for pid in (67943,68000,68527,68584,67935,68001,68526,68592):
            q=a[a['id']==pid][0];v['focus'].append(dict(id=pid,rank=int(q['rank']),relative_P_error=float(q['metrics'][2]),
                orthogonality_error=float(q['diagnostics'][2]),quad_E_over_thermal=float(q['metrics'][13])))
    result['cases'].append(v)
if 'fine_fresh' in cases and 'fine_mpi4' in cases:
    a,b=cases['fine_fresh'],cases['fine_mpi4'];assert a.shape==b.shape
    result['MPI1_4_all_exported_fields_bitmatch']=bool(a.tobytes()==b.tobytes())
    result['MPI1_4_field_matches']={k:bool(a[k].tobytes()==b[k].tobytes()) for k in a.dtype.names}
result['classification']='these projection records alone do not certify Euler transport, runtime state, shocks or long-time behavior; see separate coupled-RHS artifacts when present'
(dest/'analysis.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
