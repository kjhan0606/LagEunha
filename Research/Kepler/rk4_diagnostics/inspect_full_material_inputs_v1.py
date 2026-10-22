#!/usr/bin/env python3
"""Bounded inspection of retained full native geometry, not a simulation.

Exact equality only. A singly incident interior edge is a missing canonical
topology connection, not permission to merge coordinates with a tolerance.
"""
from pathlib import Path
import collections
import hashlib
import json
import struct
import sys
import numpy as np

root, dest = map(Path, sys.argv[1:]); dest.mkdir()
saved = root / 'dg_projection_v2_p2/fine_fresh'
inputs = [saved / 'dg_project_r0.bin', saved / 'dg_polygons_r0.bin',
          root.parent / 'core_refinement/fine_integrated/geometry.npz',
          root.parent / 'core_refinement/fine_integrated/khout.000001.dat']
hashes = {}
for p in inputs:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    hashes[str(p)] = h.hexdigest()
assert hashes[str(inputs[3])] == '46253652acfdd282904f77b1ac230708aa29a315aca7dcc96c89cb07342d2437'
xy = np.load(inputs[2])['xy']; assert xy.shape == (71080, 2)
centres = {}
with inputs[0].open('rb') as f:
    h = struct.unpack('<8Q', f.read(64)); assert h[0] == 0x4746534447503031 and h[2] == 71080
    for _ in range(h[2]):
        row = f.read(h[6]); assert len(row) == h[6]
        ident, x, y = struct.unpack_from('<Qdd', row)
        assert ident not in centres and np.array_equal(xy[ident], [x,y])
        centres[ident] = (x,y)
    assert not f.read(1)
polygons = []; edge_count = 0; areas = []; zero_edges = 0
with inputs[1].open('rb') as f:
    h = struct.unpack('<4Q', f.read(32)); assert h[0] == 0x4746534447504f31 and h[2] == 71080
    for _ in range(h[2]):
        ident, count, valid = struct.unpack('<Qii', f.read(16)); assert valid == 1 and 3 <= count <= 48
        local = np.frombuffer(f.read(count*16), dtype='<f8').reshape(count,2)
        p = local + np.array(centres[ident]) - 12
        assert np.isfinite(p).all()
        keep = np.any(p != np.roll(p,-1,axis=0),axis=1)
        zero_edges += int((~keep).sum()); p = p[keep]; assert len(p) >= 3
        area = .5 * np.sum(local[:,0]*np.roll(local[:,1],-1)-local[:,1]*np.roll(local[:,0],-1))
        assert area > 0; areas.append(float(area));polygons.append((ident,p));edge_count += len(p)
    assert not f.read(1)
corners = np.concatenate([p for _,p in polygons]); vertices, inverse = np.unique(corners, axis=0, return_inverse=True)
edges = collections.defaultdict(list); vertex_cells = collections.defaultdict(set)
offset = 0
for ident,p in polygons:
    ids = inverse[offset:offset+len(p)]; offset += len(p)
    for j,a in enumerate(ids):
        b = ids[(j+1)%len(ids)]; key = tuple(sorted((int(a),int(b))))
        edges[key].append(int(ident));vertex_cells[int(a)].add(int(ident))
inside_single = []; exterior = []; nonmanifold = []
for (a,b),cells in edges.items():
    if len(cells) == 1:
        p,q = vertices[[a,b]]
        on_box = any(p[k] == q[k] == side for k in (0,1) for side in (-12.,12.))
        if on_box: exterior.append((a,b,cells[0]))
        else: inside_single.append({'vertices':[p.tolist(),q.tolist()],'cell':cells[0]})
    elif len(cells) != 2: nonmanifold.append({'edge':[a,b],'cells':cells})
# Four subtriangles per native generator fan triangle is the retained patch
# reference refinement rule. Counts below are geometry-derived estimates;
# actual streamed quadrature counts require its native fixed-rule constructor.
ntri = edge_count * 4
raw_blocks = ntri * 18*18*16
report = {
    'scope':'all retained71080 native polygons; exact coordinate incidence and full-state memory sizing only; no new tessellation, field queries or evolution',
    'input_sha256':hashes,
    'cells':len(polygons),'native_nonzero_directed_edges':edge_count,
    'native_zero_length_edges_removed_exactly':zero_edges,
    'polygon_area_sum':float(np.sum(areas)),
    'exact_coordinate_vertices':len(vertices),'exact_coordinate_edges':len(edges),
    'paired_exact_edges':sum(len(c)==2 for c in edges.values()),
    'box_boundary_edges':len(exterior),'unpaired_interior_edges':len(inside_single),
    'nonmanifold_exact_edges':len(nonmanifold),
    'exact_coordinate_topology_ready':not inside_single and not nonmanifold,
    'interior_edge_witnesses':inside_single[:16],'nonmanifold_witnesses':nonmanifold[:16],
    'level1_triangle_count_if_same_fan_refinement':ntri,
    'one_global_raw_binary80_ST18_mass_block_bytes':raw_blocks,
    'minimum_regular64x28_points_assuming_one_segment_per_triangle':ntri*64*28,
    'materialized_STPoint_plus_two_GPOutput_bytes_per_point_estimate':1600,
    'must_stream_not_materialize_full_point_output_arrays':True,
    'canonical_neighbor_or_support_ID_registry_available_in_polygon_file':False,
    'proximity_merging_performed':False,
    'native_material_mode_implemented':False,
    'full_disk_evolution_performed':False,
}
(dest/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if 'witness' not in k},indent=2))
