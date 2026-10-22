#!/usr/bin/env python3
"""Saved-output analysis only; no field/force evaluations or evolution."""
from pathlib import Path
import csv
import hashlib
import json
import mmap
import struct
import numpy as np

ROOT=Path(__file__).resolve().parent
RUN=ROOT/'material_compact_initial_v1_ada_mpi4/run'
DEST=ROOT/'material_compact_initial_failure_analysis_v1'
DEST.mkdir(exist_ok=False)
reference=ROOT/'material_full_reference_v2_h200_mpi4/run/material_canonical_reference.bin'
with reference.open('rb') as f:
    header=struct.unpack('<8Q',f.read(64))
    vertex=np.fromfile(f,dtype='<f8',count=2*header[3]).reshape(-1,2)
    triangle=np.fromfile(f,dtype='<u8',count=4*header[4]).reshape(-1,4)
parent=np.column_stack((triangle[::4,0],triangle[::4,1],triangle[1::4,2],triangle[2::4,3]))
core=np.full(header[2],-1,dtype=np.int8)
for rank in range(4):
    p=ROOT/f'material_compact_mass_v1_mpi4/run/compact_mass_elements_r{rank}.bin'
    with p.open('rb') as f,mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as data:
        h=struct.unpack_from('<8Q',data,0);at=64
        for _ in range(h[4]):
            first,count,np_=struct.unpack_from('<3Q',data,at);at+=24
            cell=int(triangle[first,0]);flag=int(count==1)
            assert core[cell] in (-1,flag)
            core[cell]=flag;n=np_+2;at+=4*np_+16*(n*n+2*n)
        assert at+8==len(data)
assert np.count_nonzero(core==1)==5524 and np.all(core>=0)
edges=np.concatenate([parent[:,[1,2]],parent[:,[2,3]],parent[:,[3,1]]])
parents=np.tile(np.arange(parent.shape[0]),3)
edges.sort(axis=1)
order=np.lexsort((edges[:,1],edges[:,0]));edges=edges[order];parents=parents[order]
adjacent=np.all(edges[1:]==edges[:-1],axis=1)
left=parents[:-1][adjacent];right=parents[1:][adjacent]
different=core[parent[left,0]]!=core[parent[right,0]]
interface_cells=set(map(int,np.concatenate((parent[left[different],0],parent[right[different],0]))))
assert np.count_nonzero(different)==240
centroid=np.zeros((header[2],2))
centroid[parent[:,0]]=vertex[parent[:,1]]
radius=np.hypot(centroid[:,0],centroid[:,1])
rows=[]
for rank in range(4):
    with (RUN/f'full_initial_curvature_cells_r{rank}.csv').open() as f:
        for record in csv.DictReader(f):
            c=int(record['cell']);rows.append({'cell':c,'rank':rank,'x':float(centroid[c,0]),'y':float(centroid[c,1]),'radius':float(radius[c]),'protected_core':bool(core[c]),'adjacent_to_core_interface':c in interface_cells,'max_chi':float(record['max_chi']),'mean_chi':float(record['mean_chi'])})
assert len(rows)==header[2] and len({r['cell'] for r in rows})==len(rows)
failed=[r for r in rows if r['max_chi']>.002]
categories={}
for key,predicate in [('protected',lambda r:r['protected_core']),('coarse',lambda r:not r['protected_core']),('interface_adjacent',lambda r:r['adjacent_to_core_interface']),('not_interface_adjacent',lambda r:not r['adjacent_to_core_interface'])]:
    selected=[r for r in rows if predicate(r)]
    categories[key]={'cells':len(selected),'point_max_failures':sum(r['max_chi']>.002 for r in selected),'max_chi':max(r['max_chi'] for r in selected),'maximum_cell_mean_chi':max(r['mean_chi'] for r in selected)}
bins=[]
for lo,hi in [(0,.45),(.45,.49),(.49,.5),(.5,.51),(.51,.55),(.55,2),(2,8),(8,float('inf'))]:
    selected=[r for r in rows if lo<=r['radius']<hi]
    bins.append({'r_min':lo,'r_max':hi if np.isfinite(hi) else None,'cells':len(selected),'point_max_failures':sum(r['max_chi']>.002 for r in selected),'max_chi':max((r['max_chi'] for r in selected),default=None),'maximum_cell_mean_chi':max((r['mean_chi'] for r in selected),default=None)})
old=json.loads((ROOT/'material_full_curvature_v1_h200_mpi4/run/full_initial_curvature_result.json').read_text())
result={'scope':'all saved compact cell summaries + retained topology/mass group metadata, no new operator query','cells_checked':len(rows),'point_max_failures':len(failed),'cell_mean_failures':sum(r['mean_chi']>.001 for r in rows),'maximum_over_bound':max(r['max_chi'] for r in rows)/.002,'old_full_space_max_chi':old['max_chi'],'categories':categories,'radial_bins_by_cell_centroid':bins,'worst_cells':sorted(rows,key=lambda r:r['max_chi'],reverse=True)[:12],'failure_is_not_a_GPU_CPU_disagreement':True,'causality_limit':'Interface localization does not prove that hanging constraints alone cause every error; the polynomial and residual trial spaces also changed. No core/order/enrichment sweep is run.','fixed_candidate_rejected':True,'new_accepted_step':False,'physical_long_term_heating_resolved':False}
(DEST/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
with (DEST/'failed_cells.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(sorted(failed,key=lambda r:r['cell']))
with (DEST/'evidence.sha256').open('x') as f:
    for p in [Path(__file__),DEST/'analysis.json',DEST/'failed_cells.csv',RUN/'compact_initial_solve_result.json',RUN/'full_initial_curvature_result.json']:
        f.write(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p)+'\n')
print(json.dumps({'failed_cells':len(failed),'categories':categories,'radial_bins':bins},indent=2))
