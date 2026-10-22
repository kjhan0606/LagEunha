"""Independent read/assessment of native Euler projection packets."""
from pathlib import Path
import argparse,json,struct
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=2);parser.add_argument('--version',default='v1');args=parser.parse_args()
dest=ROOT/f'rk4_diagnostics/dg_projection_{args.version}_p{args.order}'
manifest=json.loads((dest/'manifest.json').read_text())
def load_case(d):
    arrays=[]
    for path in sorted(d.glob('dg_project_r*.bin')):
        with path.open('rb') as f:magic,version,n,order,nb,nm,size,rule=struct.unpack('<8Q',f.read(64))
        assert magic==0x4746534447503031 and version==1 and order==args.order and nm==17 and size==288+40*nb
        dt=np.dtype(dict(names=['id','x','y','basis','mean','c','integral','coarse_integral','metrics','flags'],
            formats=['<u8','<f8','<f8',('<f8',6),('<f8',nb),('<f8',(4,nb)),('<f8',4),('<f8',4),('<f8',17),('<i4',4)],
            offsets=[0,8,16,24,72,72+8*nb,72+40*nb,104+40*nb,136+40*nb,272+40*nb],itemsize=size))
        assert path.stat().st_size==64+n*size
        arrays.append(np.fromfile(path,dtype=dt,count=n,offset=64))
    a=np.concatenate(arrays);a.sort(order='id');assert np.unique(a['id']).size==len(a)
    return a
out=dict(order=args.order,scope='initial conservative Euler projection, not time evolution',cases=[])
loaded={}
for case in manifest['cases']:
    a=load_case(dest/case['name']);loaded[case['name']]=a;m=a['metrics'];status=a['flags'][:,0]
    codes,counts=np.unique(status,return_counts=True)
    row=dict(name=case['name'],n=len(a),statuses={str(k):int(v) for k,v in zip(codes,counts)},
        min_density=float(m[:,0].min()),min_pressure=float(m[:,1].min()),max_pressure=float(m[:,2].max()),
        nonpositive_density_cells=int((a['flags'][:,1]>0).sum()),nonpositive_pressure_cells=int((a['flags'][:,2]>0).sum()),
        max_relative_pressure_error=float(m[:,3].max()),volume_L1_relative_pressure_error=float(np.average(m[:,4],weights=a['basis'][:,0])),
        max_pressure_integration_control_relative=float(m[:,6].max()),
        max_conserved_integral_relative_error=[float(x) for x in m[:,13:17].max(axis=0)],
        max_quad_mass_relative_change=float(m[:,9].max()),max_quad_E_change_over_physical_U=float(m[:,12].max()),
        conservative_totals=[float(x) for x in a['integral'].sum(axis=0)],radial=[])
    if case['name'].startswith('fine'):
        r=np.hypot(a['x']-12,a['y']-12)
        for lo,hi in [(0,.047),(.047,.1),(.1,.2),(.2,.5),(.5,1),(1,2),(2,4),(4,8),(8,20)]:
            sel=(r>=lo)&(r<hi)
            if not sel.any():continue
            row['radial'].append(dict(rlo=lo,rhi=hi,n=int(sel.sum()),nonpositive_pressure_cells=int((a['flags'][sel,2]>0).sum()),
                max_P_error=float(m[sel,3].max()),median_cell_max_P_error=float(np.median(m[sel,3])),
                max_quad_pressure_change=float(m[sel,6].max()),max_conservation_error=float(m[sel,13:17].max()),
                volume_L1_P_error=float(np.average(m[sel,4],weights=a['basis'][sel,0]))))
        worst=np.argsort(m[:,3])[-12:][::-1]
        row['worst_cells']=[dict(id=int(a['id'][i]),r=float(r[i]),min_P=float(m[i,1]),max_P=float(m[i,2]),
            max_P_error=float(m[i,3]),quad_P_change=float(m[i,6]),worst_x=float(a['x'][i]+m[i,7]),
            worst_y=float(a['y'][i]+m[i,8]),worst_kind=int(a['flags'][i,3])) for i in worst]
    out['cases'].append(row)
if 'fine_fresh' in loaded and 'fine_mpi4' in loaded:
    a,b=loaded['fine_fresh'],loaded['fine_mpi4'];assert np.array_equal(a['id'],b['id'])
    out['MPI1_4_bitwise_fields']={key:bool(np.array_equal(a[key],b[key])) for key in a.dtype.names}
    assert all(out['MPI1_4_bitwise_fields'].values())
out['classification']='spatial candidate fails positivity' if any(c['nonpositive_pressure_cells'] or c['nonpositive_density_cells'] for c in out['cases']) else 'positivity observed; pressure accuracy and reference integration still require evaluation'
(dest/'analysis.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
