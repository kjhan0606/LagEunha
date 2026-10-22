#!/usr/bin/env python3
"""Read-only saved full native geometry/plan audit; no hydro/field evaluation."""
import csv
import hashlib
import json
import os
import struct
import sys
from pathlib import Path
import numpy as np

run, base, output = map(Path, sys.argv[1:])
output.mkdir()
raw = run / 'run'
manifest = json.loads((raw / 'material_full_reference.json').read_text())
assert manifest['cells'] == 71080 and manifest['accepted_steps'] == 0
assert manifest['reference_preparation_complete'] and not manifest['full_field_initialized']
assert not manifest['distributed_mass_operator_implemented'] and not manifest['long_Kepler_validated']
def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()
for m in ['source.sha256', 'retained_input.sha256', 'binary.sha256', 'retained_build_objects.sha256']:
    for line in (run / m).read_text().splitlines():
        expected, name = line.split('  ', 1)
        assert sha(Path(name)) == expected, name
geom = np.load(base / 'core_refinement/fine_integrated/geometry.npz')
xy, old_volume = geom['xy'], geom['volume']
n = len(xy)
assert n == 71080
seen = np.zeros(n, dtype=bool)
owned = np.empty(n, dtype=np.int32)
native_volume = np.zeros(n)
native_mass = np.zeros(n)
native_U = np.zeros(n)
raw_face_ends = 0
external_ends = 0
UUID = set()
for rank in range(4):
    runtime = json.loads((raw / f'material_runtime_r{rank}.json').read_text())
    assert runtime['rank'] == rank and runtime['MPI_ranks'] == 4
    assert runtime['geometry_from_original_CUDA']
    UUID.add(runtime['device_UUID_hex'])
    with (raw / f'material_native_geometry_r{rank}.bin').open('rb') as f:
        h = struct.unpack('<8Q', f.read(64))
        assert h[0] == 0x4d46524157474531 and h[1] == 1 and h[3:6] == (80, 56, rank)
        assert h[2] == runtime['owned_cells']
        for _ in range(h[2]):
            ident, count, x, y, mass, volume, vx, vy, U, P = struct.unpack('<2Q8d', f.read(80))
            assert ident < n and not seen[ident] and count >= 3
            assert (x, y) == tuple(xy[ident])
            assert mass > 0 and U > 0 and volume > 0
            seen[ident] = True
            owned[ident] = rank
            native_volume[ident], native_mass[ident], native_U[ident] = volume, mass, U
            for _ in range(count):
                neighbor, nx, ny, ax, ay, bx, by = struct.unpack('<Q6d', f.read(56))
                assert neighbor < n and all(np.isfinite([nx, ny, ax, ay, bx, by]))
                assert nx * 512 == round(nx * 512) and ny * 512 == round(ny * 512)
                shift = (np.array([nx, ny]) - xy[neighbor]) / 24
                assert np.array_equal(shift, np.rint(shift)) and np.max(abs(shift)) <= 1
                if np.any(shift):
                    external_ends += 1
                raw_face_ends += 1
        assert f.read() == b''
assert seen.all() and len(UUID) == 1  # four MPI ranks sharing one card
assert raw_face_ends == manifest['raw_face_ends']
with (raw / 'material_canonical_reference.bin').open('rb') as f:
    h = struct.unpack('<8Q', f.read(64))
    assert h[0] == 0x4d4643414e4f4e31 and h[1] == 1 and h[2] == n
    assert h[6:] == (16, 32)
    vertex = np.fromfile(f, '<f8', h[3] * 2).reshape(-1, 2)
    tri = np.fromfile(f, '<u8', h[4] * 4).reshape(-1, 4)
    assert f.read() == b''
assert np.isfinite(vertex).all() and np.max(tri[:, 0]) < n and np.max(tri[:, 1:]) < len(vertex)
a, b, c = (vertex[tri[:, j]] for j in [1, 2, 3])
area = ((b[:, 0].astype(np.longdouble) - a[:, 0]) * (c[:, 1].astype(np.longdouble) - a[:, 1]) -
        (b[:, 1].astype(np.longdouble) - a[:, 1]) * (c[:, 0].astype(np.longdouble) - a[:, 0])) / 2
assert np.min(area) > 0
cell_area = np.zeros(n, dtype=np.longdouble)
np.add.at(cell_area, tri[:, 0], area)
assert abs(np.sum(cell_area) - 576) < 1e-11
relative = np.max(abs(cell_area - old_volume) / old_volume)
assert relative < 1e-12, relative
assert np.max(abs(native_volume - old_volume) / old_volume) < 1e-12
# Whole refined topology, using integer vertex IDs; no coordinate matching.
edges = np.concatenate([tri[:, [1, 2]], tri[:, [2, 3]], tri[:, [3, 1]]])
edges.sort(axis=1)
edge, multiplicity = np.unique(edges, axis=0, return_counts=True)
assert np.all((multiplicity == 1) | (multiplicity == 2))
assert len(edge) == h[5] == manifest['reference_edges']
assert len(vertex) - len(edge) + len(tri) == 1
boundary = edge[multiplicity == 1]
left, right = vertex[boundary[:, 0]], vertex[boundary[:, 1]]
on_box = ((left[:, 0] == -12) & (right[:, 0] == -12) |
          (left[:, 0] == 12) & (right[:, 0] == 12) |
          (left[:, 1] == -12) & (right[:, 1] == -12) |
          (left[:, 1] == 12) & (right[:, 1] == 12))
assert on_box.all() and len(boundary) == 2 * manifest['external_physical_faces']
del a, b, c, edges, edge
# Read the compact plans without expanding a single physical quadrature point.
planfile = raw / 'material_reference_quadrature_plan.bin'
with planfile.open('rb') as f:
    qh = struct.unpack('<6Q', f.read(48))
assert qh == (0x4d46515541443131, 1, len(tri), 64, 28, 288)
assert planfile.stat().st_size == 48 + len(tri) * 288
dtype = np.dtype([('anchor', '<f8', 2), ('cut', '<f8', 33), ('ns', '<i4'), ('anchor_edge', '<i4')])
plan = np.memmap(planfile, mode='r', offset=48, dtype=dtype, shape=(len(tri),))
assert np.isfinite(plan['anchor']).all() and np.all((plan['ns'] >= 1) & (plan['ns'] <= 32))
slots = int(np.sum(plan['ns'], dtype=np.uint64)) * 3 * 64 * 28
assert slots == manifest['fixed_quadrature_slots']
assert np.all(plan['cut'][:, 0] == 0)
assert np.all(plan['cut'][np.arange(len(plan)), plan['ns']] == 1)
files = [raw / f'material_native_geometry_r{r}.bin' for r in range(4)]
files += [raw / f'material_runtime_r{r}.json' for r in range(4)]
files += [raw / n for n in ['material_canonical_reference.bin', 'material_reference_quadrature_plan.bin', 'material_reference_cells.csv', 'material_full_reference.json']]
files += [run / 'native.log', run / 'native_exit_code.txt']
with (output / 'verified_output.sha256').open('x') as f:
    for p in files:
        f.write(f'{sha(p)}  {p}\n')
# Save the original MPI owner map for subsequent distributed state construction.
np.save(output / 'cell_native_owner.npy', owned)
report = {
    'scope': 'independent saved full native geometry/topology/compact-plan audit; no field quadrature or evolution',
    'owned_cells': n, 'MPI_ranks': 4, 'distinct_GPU_UUIDs': len(UUID),
    'raw_face_ends': raw_face_ends,
    'canonical_refined_triangles': len(tri),
    'integer_ID_manifold_and_Euler_pass': True,
    'refined_external_edges': len(boundary),
    'all_external_edges_on_original_square': True,
    'all71080_cell_areas_match_saved_geometry_relative_max': str(relative),
    'all71080_native_centers_match_saved_IC_bitwise': True,
    'actual_positive_triangle_area_min': str(np.min(area)),
    'compact_plan_slots': slots,
    'physical_field_or_RHS_queries_replayed': 0,
    'audit_pass': True,
    'full_material_field_initialized': False,
    'accepted_steps': 0,
    'long_Kepler_validated': False,
}
(output / 'contract.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
