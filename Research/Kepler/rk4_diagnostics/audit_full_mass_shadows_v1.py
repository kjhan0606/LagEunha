#!/usr/bin/env python3
"""Recompute normalized matrix comparisons from saved witnesses only.

No quadrature, density, basis-point evaluation or simulation is repeated.
"""
import hashlib
import json
import struct
import sys
from pathlib import Path
import numpy as np

run, dest = map(Path, sys.argv[1:]);dest.mkdir()
assert np.dtype(np.longdouble).itemsize == 16
LD = np.longdouble
summary = json.loads((run/'contract.json').read_text())
assert summary['actual_positive_mass_queries'] == 2040256512
assert summary['global_cells'] == 71080 and summary['global_dofs'] == 10816441
assert summary['all_local_mass_Cholesky_SPD'] and not summary['full_disk_evolution_performed']
ranks = [json.loads((run/f'rank{r}.json').read_text()) for r in range(4)]
assert sum(r['actual_positive_mass_queries'] for r in ranks) == 2040256512
assert sum(r['local_triangles'] for r in ranks) == summary['global_triangles']
assert sum(r['owned_dofs'] for r in ranks) == summary['global_dofs']
devices = [json.loads((run/f'device_rank{r}.json').read_text()) for r in range(4)]
assert all(d['visible_devices'] == 1 for d in devices)
assert len({d['UUID_hex'] for d in devices}) == 1
rows = []
# Exactly one known run directory, no recursive file search.
for p in sorted(run.iterdir()):
    if not (p.name.startswith('mass_shadow_rank') and p.suffix == '.bin'):
        continue
    data = p.read_bytes();h = struct.unpack_from('<6Q',data)
    assert h[0:2] == (0x4746534d53484431,1) and h[3:] == (1064,16,18)
    t = data[48:48+1064];bubble,n,cell = struct.unpack_from('<3i',t,1052)
    assert 15 <= n <= 18 and bubble == n-2
    off=48+1064
    gpu=np.frombuffer(data,dtype=np.longdouble,count=324,offset=off).copy().reshape(18,18)
    off+=324*16;cpu=np.frombuffer(data,dtype=np.longdouble,count=324,offset=off).copy().reshape(18,18)
    off+=324*16;mass=np.frombuffer(data,dtype=np.longdouble,count=1,offset=off)[0]
    assert len(data)==off+16 and np.isfinite(gpu).all() and np.isfinite(cpu).all() and mass>0
    shift=np.asarray(gpu[9,bubble:bubble+2]/gpu[9,9],dtype=np.float64)
    H=np.zeros((2,2),dtype=np.longdouble)
    for i in range(2):
        for j in range(2):
            H[i,j]=gpu[bubble+i,bubble+j]-LD(shift[i])*gpu[9,bubble+j]-LD(shift[j])*gpu[bubble+i,9]+LD(shift[i])*LD(shift[j])*gpu[9,9]
    pivot=int(H[0,0]<H[1,1]);other=1-pivot
    assert H[pivot,pivot]>0
    s0=float(np.sqrt(H[pivot,pivot]/mass));beta=float(H[pivot,other]/H[pivot,pivot])
    variance=H[other,other]-2*LD(beta)*H[pivot,other]+LD(beta)*LD(beta)*H[pivot,pivot]
    assert variance>0;s1=float(np.sqrt(variance/mass))
    encode=np.zeros((2,2),dtype=np.float64);encode[0,pivot]=1/s0;encode[1,pivot]=-beta/s1;encode[1,other]=1/s1
    T=np.eye(18,dtype=np.longdouble)
    for a in range(2):
        row=bubble+a;T[row,:]=0
        for b in range(2):
            T[row,bubble+b]=LD(encode[a,b]);T[row,9]-=LD(encode[a,b])*LD(shift[b])
    G=T@gpu@T.T;C=T@cpu@T.T
    ratio=np.max(np.abs(C-G))/np.max(np.abs(C))
    assert ratio < LD('1e-11'),(p.name,str(ratio))
    rows.append({'source':h[2],'witness':str(p),'sha256':hashlib.sha256(data).hexdigest(),
                 'normalized_mass_CPU_CUDA_relative_difference':str(ratio)})
assert len(rows)==14 and len({r['source'] for r in rows})==14
report={'scope':'saved full-state mass shadow matrices and rank/device accounting; no physical query replay',
        'sources_audited':len(rows),'actual_positive_mass_queries_accounted':2040256512,
        'one_real_H200_shared_by_four_MPI_ranks':True,
        'unique_owned_dofs':summary['global_dofs'],
        'max_per_source_normalized_mass_CPU_CUDA_difference':str(max(LD(r['normalized_mass_CPU_CUDA_relative_difference']) for r in rows)),
        'rows':rows,'audit_pass':True,'full_disk_evolution_verified':False}
(dest/'contract.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
