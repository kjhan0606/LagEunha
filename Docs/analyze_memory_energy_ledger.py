"""Read completed memory steps and sidecars without altering a live simulation.

Energy is hydro K+U+Phi plus auxiliary energy. Boundary corrections here cover
only the explicitly reported auxiliary and HG terms, not all hydro boundaries.
"""
import argparse,json,math,re,struct,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('folder',type=Path)
args=parser.parse_args()
def fields(line):
    return dict(re.findall(r'(\w+)=([^\s]+)',line))

pending=None;paired=[];aux_work=[];hg_work=[];hg_heat=[];initial=None
with (args.folder/'log').open() as stream:
    for line in stream:
        if not line.endswith('\n'):break # writer may still be completing a line
        if line.startswith('[ACCEPTED_STATE]'):
            assert pending is None,'Unpaired accepted step'
            pending=fields(line)
        elif line.startswith('[RMEM_ACCEPTED]'):
            q=fields(line);assert pending is not None
            assert int(pending['finite'])==1 and float(pending['dM0'])==0
            assert float(q['dt'])==float(pending['dt'])
            assert int(q['retries'])==int(pending['retries'])
            eh=float(pending['E']);ez=float(q['Eaux'])
            assert all(math.isfinite(float(q[k])) for k in ('Eaux','dEaux','dEhydro','boundary_work','augmented_residual'))
            assert ez>=0
            if initial is None:
                initial=eh-float(q['dEhydro'])+ez-float(q['dEaux'])
            aux_work.append(float(q['boundary_work']))
            paired.append({'step':int(pending['step']),'elapsed':float(pending['elapsed']),
                'hydro_energy':eh,'auxiliary_energy':ez,'augmented_energy':eh+ez,
                'augmented_energy_change':eh+ez-initial,
                'auxiliary_boundary_work_cumulative':math.fsum(aux_work),
                'change_minus_auxiliary_boundary_work':eh+ez-initial-math.fsum(aux_work),
                'retries':int(q['retries'])})
            pending=None
        elif line.startswith('[ACCEPTED_HG_BOUNDARY]'):
            q=fields(line);hg_work.append(float(q['work']));hg_heat.append(float(q['heat_out']))

ic,tic=read(Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat'))
assert tic==0
r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'inner_background':(r0<.5)&(ic['den']<.5),
       'inner_disk_edge':(r0<.6)&(ic['den']>.5),
       'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5),
       'inner_shear_band':(r0>=.3)&(r0<.5)}
assert {key:int(mask.sum()) for key,mask in masks.items()}=={
    'inner_background':804,'inner_disk_edge':390,'disk_core':10304,'inner_shear_band':521}
sidecars=[];skipped=[]
for path in sorted(args.folder.glob('khout.*.dat')):
    side=Path(str(path)+'.rmem')
    if not side.exists():skipped.append(path.name);continue
    data=side.read_bytes();assert len(data)>=56
    magic,version,endian,n,hydro_hash,checksum=struct.unpack('<6Q',data[:48])
    assert (magic,version,endian)==(0x524d454d53544154,1,0x0102030405060708)
    assert len(data)==56+24*n
    h=14695981039346656037
    for byte in data[56:]:h=((h^byte)*1099511628211)&0xffffffffffffffff
    assert h==checksum,'Sidecar record checksum mismatch'
    a,t=read(path);assert t==struct.unpack('<d',data[48:56])[0] and len(a)==n
    z=np.frombuffer(data,dtype=[('id','<u8'),('x','<f8'),('y','<f8')],offset=56)
    assert np.array_equal(z['id'],a['id']&np.uint64(72057594037927935))
    assert np.array_equal(z['id'],ic['id']&np.uint64(72057594037927935))
    assert np.isfinite(z['x']).all() and np.isfinite(z['y']).all()
    aux=.5*(z['x']**2+z['y']**2);total=float(aux.sum())
    xx=a['x']-3;yy=a['y']-3
    omega=np.sqrt(gravity_coefficient(np.linalg.norm(np.column_stack((xx,yy)),axis=1)))
    departure=.5*a['mass']*((a['vx']+omega*yy)**2+(a['vy']-omega*xx)**2)
    regions={}
    for name,mask in masks.items():
        ez=float(aux[mask].sum());ed=float(departure[mask].sum());thermal=float(a['ie'][mask].sum())
        regions[name]={'cells':int(mask.sum()),'auxiliary_energy':ez,
                       'departure_kinetic_energy':ed,'thermal_energy':thermal,
                       'auxiliary_over_departure':ez/ed if ed>0 else None,
                       'auxiliary_over_thermal':ez/thermal if thermal>0 else None,
                       'fraction_of_global_auxiliary':ez/total if total>0 else None}
    sidecars.append({'snapshot':str(path),'time':t,
        'auxiliary_energy':total,'material_regions':regions,
        'record_checksum_verified':True,'time_and_ids_verified':True,
        'stored_hydro_fingerprint':hex(hydro_hash),'hydro_fingerprint_independently_verified':False})
out={'completed_memory_steps':len(paired),'latest':paired[-1] if paired else None,
     'total_retries':sum(q['retries'] for q in paired),'unpaired_tail':pending is not None,
     'hg_boundary_records':len(hg_work),'hg_boundary_work_cumulative':math.fsum(hg_work),
     'hg_heat_out_cumulative':math.fsum(hg_heat),'sidecars':sidecars,
     'snapshots_without_sidecar':skipped,
     'scope':'Partial boundary accounting, not a closed-system conservation or physical stability pass. Regions use original IC material IDs and overlap; auxiliary energy is assigned to its owner, not a local energy-flux balance. Missing sidecar at original t0 is expected; later files may be in publication. Hydro fingerprint validation remains performed by simulator restart, not this reader.'}
print(json.dumps(out,indent=2))
