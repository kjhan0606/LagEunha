#!/usr/bin/env python3
"""Saved full-state/data audit only. No physical quadrature or solver replay."""
import ctypes as C
import hashlib
import json
import struct
import sys
from pathlib import Path
import numpy as np

run, base, out = map(Path, sys.argv[1:]); out.mkdir()
r = json.loads((run / 'run/full_initial_result.json').read_text())
assert r['initial_operator_algebra_pass'] and r['accepted_steps'] == 0 and not r['long_Kepler_validated']
assert max(r['true_equation_residual']) < 1e-16
assert r['global_energy_rate_error_over_U'] < 1e-8 and r['global_angular_rate_error_over_2K'] < 1e-11
assert r['mass_metric_vs_point_K_error_over_U'] < 1e-8
assert r['CPU_CUDA_mass_shadow_relative'] < 1e-11 and r['CPU_CUDA_force_shadow_relative'] < 1e-11
assert np.dtype(np.longdouble).itemsize == 16 and np.finfo(np.longdouble).nmant == 63

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()
for name in ['source.sha256', 'retained_input.sha256', 'binary.sha256', 'retained_stage_builds.sha256']:
    for line in (run / name).read_text().splitlines():
        expected, p = line.split('  ', 1); assert sha(Path(p)) == expected, p

class Edge(C.Structure):
    _fields_ = [('A', C.c_double*2), ('D', C.c_double*2), ('qnode', (C.c_double*2)*4),
                ('eps', C.c_double), ('encode', (C.c_double*2)*2), ('decode', (C.c_double*2)*2),
                ('rank', C.c_int), ('id', C.c_int*2), ('poly', C.c_int*4)]
class Tri(C.Structure):
    _fields_ = [('vertex', (C.c_double*2)*3), ('grad', (C.c_double*2)*3), ('qnode', (C.c_double*2)*10),
                ('eps', C.c_double), ('bubble_shift', C.c_double*2), ('bubble_encode', (C.c_double*2)*2),
                ('bubble_decode', (C.c_double*2)*2), ('edge', Edge*3), ('endpoint', (C.c_int*2)*3),
                ('id', C.c_int*18), ('offset', C.c_int*3), ('bubble_offset', C.c_int), ('n', C.c_int), ('cell', C.c_int)]
assert C.sizeof(Edge) == 200 and C.sizeof(Tri) == 1064
tri_type = np.dtype({'names':['vertex','id','n','cell'], 'formats':[('<f8',(3,2)),('<i4',18),'<i4','<i4'],
                     'offsets':[Tri.vertex.offset,Tri.id.offset,Tri.n.offset,Tri.cell.offset], 'itemsize':1064})
dd = np.dtype([('hi','<f8'),('lo','<f8')])
raw_type = np.dtype([('mass',dd,324),('force',dd,(2,18)),('metric',dd,16),('minJ','<f8'),
                     ('points','<u8'),('bad_point','<i4'),('bad_kind','<i4')])
assert raw_type.itemsize == 6040

reference = base/'rk4_diagnostics/material_full_reference_v2_h200_mpi4/run/material_canonical_reference.bin'
with reference.open('rb') as f: refh = struct.unpack('<8Q', f.read(64))
vertices = np.memmap(reference,mode='r',offset=64,dtype='<f8',shape=(refh[3],2))
ref_tri = np.memmap(reference,mode='r',offset=64+refh[3]*16,dtype='<u8',shape=(refh[4],4))
nc, nt, nd = refh[2], refh[4], r['global_scalar_dofs']
seen = np.full(nd,-1,dtype=np.int8); triangle_seen = np.zeros(nt,dtype=bool)
full_coeff = np.zeros((nd,4),dtype=np.float64)
mass = np.zeros(nc,dtype=np.longdouble); quad_K = mass.copy(); quad_W = mass.copy()
positive = 0; minimum_J = 1e99; files=[]; views=[]
UUID=set()

for rank in range(4):
    runtime=json.loads((run/f'run/full_initial_runtime_r{rank}.json').read_text());UUID.add(runtime['uuid'])
    p=run/f'run/full_conditioned_basis_r{rank}.bin';files.append(p)
    with p.open('rb') as f:h=struct.unpack('<12Q',f.read(96))
    assert h[:2] == (0x4d464f574e454431,1) and h[2:6] == (refh[3],refh[5],nt,nd) and h[9:] == (1064,32,288)
    ncell,ntri,nlocal=h[6:9];at=96
    gid=np.memmap(p,mode='r',offset=at,dtype='<u8',shape=(nlocal,));at+=nlocal*8
    owner=np.memmap(p,mode='r',offset=at,dtype='<i4',shape=(nlocal,));at+=nlocal*4
    cell_id=np.memmap(p,mode='r',offset=at,dtype='<u8',shape=(ncell,));at+=ncell*8
    source=np.memmap(p,mode='r',offset=at,dtype='<u8',shape=(ntri,));at+=ntri*8
    source_cell=np.memmap(p,mode='r',offset=at,dtype='<u8',shape=(ntri,));at+=ntri*8
    at+=ntri*3*8
    triangle=np.memmap(p,mode='r',offset=at,dtype=tri_type,shape=(ntri,));at+=ntri*1064
    coeff=np.memmap(p,mode='r',offset=at,dtype='<f8',shape=(nlocal,4));at+=nlocal*32
    at+=ntri*288+refh[5]
    assert at==p.stat().st_size and np.isfinite(coeff).all()
    assert np.all(np.diff(gid)>0) and np.all(gid<nd) and np.all((owner>=0)&(owner<4))
    own=owner==rank;assert np.all(seen[gid[own]]==-1);seen[gid[own]]=rank;full_coeff[gid[own]]=coeff[own]
    assert not triangle_seen[source].any();triangle_seen[source]=True
    assert np.array_equal(source_cell,ref_tri[source,0])
    assert np.array_equal(cell_id[triangle['cell']],source_cell)
    assert np.all((triangle['n']>=15)&(triangle['n']<=18))
    for start in range(0,ntri,4096):
        end=min(start+4096,ntri);assert np.array_equal(triangle['vertex'][start:end],vertices[ref_tri[source[start:end],1:]])
    views.append((gid,owner,coeff))
    rawp=run/f'run/full_raw_initial_blocks_r{rank}.bin';files.append(rawp)
    with rawp.open('rb') as f:rh=struct.unpack('<6Q',f.read(48))
    assert rh == (0x4d46494e49545231,1,ntri,6040,rank,nd) and rawp.stat().st_size==48+ntri*6040
    raw=np.memmap(rawp,mode='r',offset=48,dtype=raw_type,shape=(ntri,))
    assert np.all(raw['bad_point']==-1) and np.all(raw['minJ']>0)
    positive+=int(np.sum(raw['points'],dtype=np.uint64));minimum_J=min(minimum_J,float(np.min(raw['minJ'])))
    for start in range(0,ntri,4096):
        end=min(start+4096,ntri);metric=raw['metric'][start:end]['hi'].astype(np.longdouble)+raw['metric'][start:end]['lo'].astype(np.longdouble)
        assert np.isfinite(metric).all() and np.all(metric[:,0]>0)
        ids=source_cell[start:end];np.add.at(mass,ids,metric[:,0]);np.add.at(quad_W,ids,metric[:,1]);np.add.at(quad_K,ids,metric[:,5])
    mp=run/f'run/full_fixed_mass_r{rank}.bin';files.append(mp)
    with mp.open('rb') as f:mh=struct.unpack('<4Q',f.read(32))
    assert mh==(0x4d464d4153533031,1,ntri*324,16) and mp.stat().st_size==32+ntri*324*16
    matrix=np.memmap(mp,mode='r',offset=32,dtype=np.longdouble,shape=(ntri,18,18))
    for start in range(0,ntri,4096):
        block=matrix[start:start+4096];assert np.isfinite(block).all() and np.array_equal(block,block.transpose(0,2,1))
        for j,nval in enumerate(triangle['n'][start:start+4096]):assert np.all(np.diag(block[j])[:nval]>0)
    files.extend([run/f'run/full_unconditioned_basis_r{rank}.bin',run/f'run/full_generating_thermal_r{rank}.bin',run/f'run/full_initial_acceleration_r{rank}.bin'])

assert np.all(seen>=0) and triangle_seen.all() and positive==2040256512==r['actual_positive_points']
for gid,owner,coeff in views:assert np.array_equal(coeff,full_coeff[gid])
assert (mass>0).all()
legacy_mass=np.zeros(nc,dtype=np.float64)
native=base/'rk4_diagnostics/material_full_reference_v2_h200_mpi4/run'
for rank in range(4):
    with (native/f'material_native_geometry_r{rank}.bin').open('rb') as f:
        h=struct.unpack('<8Q',f.read(64))
        for _ in range(h[2]):
            row=struct.unpack('<2Q8d',f.read(80));legacy_mass[row[0]]=row[4];f.seek(row[1]*56,1)
relative=np.max(abs(mass-legacy_mass)/legacy_mass)
report={'scope':'independent all full-state owned/ghost/topology/raw-positive-query/fixed-mass saved-output audit; no field quadrature or mass solve replay',
        'global_dofs':nd,'source_triangles':nt,'cells':nc,'positive_points':positive,
        'all_global_DOFs_have_unique_owner':True,'all_ghost_coefficients_match_owners_bitwise':True,
        'all_owned_reference_triangles_match_saved_canonical_reference_bitwise':True,
        'all_fixed_mass_blocks_finite_symmetric_positive_diagonal':True,
        'all_raw_GPU_point_failures_absent':True,'min_J':minimum_J,
        'continuous_initial_mass':str(np.sum(mass)),'legacy_integrated_mass':str(np.sum(legacy_mass.astype(np.longdouble))),
        'max_cell_mass_relative_difference_from_retained_integrated_IC':str(relative),
        'existing_initial_mass_relative_bound':'1e-8','all_cell_initial_mass_accuracy_pass':bool(relative<1e-8),
        'distinct_GPU_UUIDs':len(UUID),'audit_pass':bool(relative<1e-8),
        'initial_curvature_evaluated_by_this_audit':False,'accepted_steps':0,'long_Kepler_validated':False}
(out/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
with (out/'verified_output.sha256').open('x') as f:
    for p in files+[run/'run/full_initial_result.json']:f.write(f'{sha(p)}  {p}\n')
print(json.dumps(report,indent=2));sys.exit(0 if report['audit_pass'] else 2)
