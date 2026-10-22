"""Independent native moment audit; execute on lageunha with one CPU thread."""
from pathlib import Path
import argparse
import json
import struct
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--version',default='v1');args=parser.parse_args()
dest=ROOT/('rk4_diagnostics/moment_initial_'+args.version)
manifest=json.loads((dest/'manifest.json').read_text())

def load_case(d):
    arrays=[];rk_offset=None
    for path in sorted(d.glob('moment_initial_r*.bin')):
        magic,version,n,audit_size,rk_offset,payload_size=struct.unpack('<6Q',path.read_bytes()[:48])
        assert magic==0x4746534d4f4d3031 and version==1 and payload_size==192 and audit_size==360
        dtype=np.dtype(dict(names=['id','x','y','payload','coarse','thermal','thermal_coarse','nonaffine','min_density','moment_residual','identity_error','valid','iterations','restart'],
            formats=['<u8','<f8','<f8',('<f8',24),('<f8',13)]+['<f8']*6+['<i4']*3,
            offsets=[0,8,16,24,216,320,328,336,344,352,360,368,372,376],itemsize=384))
        assert path.stat().st_size==48+n*dtype.itemsize
        arrays.append(np.fromfile(path,dtype=dtype,offset=48,count=n))
    if not arrays:return None,rk_offset
    a=np.concatenate(arrays);a.sort(order='id');assert np.unique(a['id']).size==len(a)
    return a,rk_offset

def primitive(p):
    a=p.astype(np.longdouble);m=a[:,0];cx=a[:,1]/m;cy=a[:,2]/m
    vx=a[:,6]/m;vy=a[:,7]/m
    sxx=a[:,3]-a[:,1]*cx;sxy=a[:,4]-a[:,1]*cy;syy=a[:,5]-a[:,2]*cy
    cxx=a[:,8]-a[:,1]*vx;cxy=a[:,9]-a[:,1]*vy
    cyx=a[:,10]-a[:,2]*vx;cyy=a[:,11]-a[:,2]*vy
    det=sxx*syy-sxy*sxy
    bxx=(cxx*syy-cyx*sxy)/det;bxy=(-cxx*sxy+cyx*sxx)/det
    byx=(cxy*syy-cyy*sxy)/det;byy=(-cxy*sxy+cyy*sxx)/det
    km=.5*(a[:,6]*vx+a[:,7]*vy)
    ka=.5*(bxx*cxx+bxy*cyx+byx*cxy+byy*cyy)
    return dict(U=a[:,12]-km-ka,Kmean=km,Kaffine=ka,B=np.array([bxx,bxy,byx,byy]).T,
        det=det,cx=cx,cy=cy,vx=vx,vy=vy)

result=dict(version=args.version,cases=[],scope='initial integration, realization, algebra, and native restart; no evolution or stability result')
loaded={}
for row in manifest['cases']:
    d=dest/row['name'];a,rk_offset=load_case(d)
    if a is None:continue
    loaded[row['name']]=a
    codes,counts=np.unique(a['valid'],return_counts=True)
    out=dict(name=row['name'],n=len(a),exit_code=row['exit_code'],
        statuses={str(k):int(v) for k,v in zip(codes,counts)},
        invalid_ids=a['id'][a['valid']!=1].tolist(),
        max_density_iterations=int(a['iterations'].max()),max_moment_residual=float(a['moment_residual'].max()))
    usable=(a['payload'][:,0]>0)
    s=a['payload'][usable,:13];r=primitive(s);aa=a[usable]
    if len(aa):
        out.update(mass=float(s[:,0].sum()),min_S_determinant=float(r['det'].min()),
            min_model_U=float(r['U'].min()),min_U_in_roundoff_units=float(np.min(r['U']/(np.finfo(float).eps*(abs(s[:,12])+r['Kmean']+r['Kaffine'])))),
            all_mass_weight_relations_exact=bool(np.array_equal(s[:,0]+aa['payload'][:,14]*aa['payload'][:,15],aa['payload'][:,13])),
            tag_exact=bool(np.all(aa['payload'][:,18]==2026102001)))
        fresh=aa['restart']==0
        if fresh.any():
            f=aa[fresh];fr={k:v[fresh] for k,v in r.items()};fs=s[fresh]
            contamination=np.asarray((fr['U']-f['thermal'])/f['thermal'],dtype=float)
            radius=np.hypot(f['x']-12,f['y']-12)
            identity=np.abs(f['identity_error'])
            scale=np.finfo(float).eps*(abs(fs[:,12])+fr['Kmean']+fr['Kaffine'])
            out.update(thermal_physical_sum=float(f['thermal'].sum()),
                nonaffine_kinetic_sum_direct=float(f['nonaffine'].sum()),
                model_U_sum=float(fr['U'].sum()),
                max_energy_identity_in_roundoff_units=float(np.max(identity/scale)),
                max_E_quadrature_change_over_physical_U=float(np.max(abs(fs[:,12]-f['coarse'][:,12])/f['thermal'])),
                max_mass_quadrature_relative_change=float(np.max(abs(fs[:,0]-f['coarse'][:,0])/fs[:,0])),radial=[])
            if row['name'] in ('rigid','shear','shear_boost'):
                exact=np.array([0,-.7,.7,0] if row['name']=='rigid' else [0,0,-.7,0])
                out['max_B_error']=float(np.max(abs(fr['B']-exact)))
                out['max_thermal_relative_error']=float(np.max(abs(contamination)))
            else:
                for lo,hi in [(0,.047),(.047,.1),(.1,.2),(.2,.5),(.5,1),(1,2),(2,4),(4,8),(8,20)]:
                    sel=(radius>=lo)&(radius<hi)
                    if not sel.any():continue
                    cc=contamination[sel]
                    out['radial'].append(dict(rlo=lo,rhi=hi,n=int(sel.sum()),
                        nonaffine_over_thermal_median=float(np.median(cc)),
                        nonaffine_over_thermal_max=float(cc.max()),
                        physical_U_sum=float(f['thermal'][sel].sum()),model_U_sum=float(fr['U'][sel].sum()),
                        quadrature_E_change_over_U_max=float(np.max(abs(fs[sel,12]-f['coarse'][sel,12])/f['thermal'][sel]))))
                worst=np.argsort(contamination)[-8:][::-1]
                out['worst_initial_thermal_contamination']=[dict(id=int(f['id'][i]),r=float(radius[i]),ratio=float(contamination[i]),
                    U_model=float(fr['U'][i]),U_physical=float(f['thermal'][i]),K_nonaffine_direct=float(f['nonaffine'][i])) for i in worst]
    # Verify the native KH full-stride writer preserved every tagged double.
    snapshots=sorted(path for path in d.glob('khout.*.dat') if not path.is_symlink())
    if snapshots:
        path=snapshots[-1]
        with path.open('rb') as fh:magic,n,stride,stress=struct.unpack('<4i',fh.read(16))
        assert magic==-1 and n==len(a) and stride==576 and stress==1
        cp_dtype=np.dtype(dict(names=['id','payload'],formats=['<u8',('<f8',24)],offsets=[8,rk_offset],itemsize=stride))
        cp=np.fromfile(path,dtype=cp_dtype,offset=16,count=n);cp.sort(order='id')
        assert np.array_equal(cp['id'],a['id'])
        out['native_snapshot_payload_bitwise_equal']=bool(np.array_equal(cp['payload'].view('<u8'),a['payload'].copy().view('<u8')))
        assert out['native_snapshot_payload_bitwise_equal']
    result['cases'].append(out)

if 'fine_fresh' in loaded and 'fine_restart_mpi4' in loaded:
    first=loaded['fine_fresh'];restart=loaded['fine_restart_mpi4']
    assert np.array_equal(first['id'],restart['id'])
    result['restart_MPI4_state_first19_doubles_bitwise_equal']=bool(np.array_equal(
        first['payload'][:,:19].copy().view('<u8'),restart['payload'][:,:19].copy().view('<u8')))
    result['restart_MPI4_full_payload_bitwise_equal']=bool(np.array_equal(
        first['payload'].copy().view('<u8'),restart['payload'].copy().view('<u8')))
    assert result['restart_MPI4_state_first19_doubles_bitwise_equal']

(dest/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
