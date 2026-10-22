#!/usr/bin/env python3
"""Saved full-RHS/cell/shadow audit and legacy-IC comparison only."""
from pathlib import Path
from decimal import Decimal as D, getcontext
import csv
import hashlib
import json
import struct
import sys
import numpy as np
getcontext().prec=70
run, root, dest = map(Path,sys.argv[1:]);dest.mkdir()
report=json.loads((run/'contract.json').read_text(),parse_float=D)
assert report['native_full_RHS_pass'] and report['actual_physical_RHS_evaluated']==1
assert not report['full_disk_RK4_connected'] and not report['long_Kepler_validated']
assert report['true_mass_equation_residual_max']<D('1e-16')
assert report['initial_chi_max']<=D('.002') and report['initial_cell_mean_chi_max']<=D('.001')
seen=set();U=K=global_E_rate=global_L_rate=D(0);max_E=max_L=D(0);shadow=[];total_point_mass=np.longdouble(0)
for rank in range(4):
    with (run/f'cell_physics_rank{rank}.csv').open() as f:
        for row in csv.DictReader(f):
            ident=int(row['cell']);assert ident not in seen;seen.add(ident)
            v={k:D(x) for k,x in row.items() if k!='cell'}
            e=v['Kdot']+v['Udot']-v['gwork']-v['external_work']-v['hybrid_work']
            l=v['Ldot']-v['gtorque']-v['external_torque']-v['hybrid_torque']
            assert v['U']>0 and v['K_metric']>0 and v['inertia']>0
            assert abs(e)/v['U']<D('1e-8') and abs(l)/(2*v['K_metric'])<D('1e-11')
            # Printed fields carry21significant digits; allow binary80
            # operation rounding when independently recombining them.
            e_scale=sum(abs(v[k]) for k in ['Kdot','Udot','gwork','external_work','hybrid_work'])
            l_scale=sum(abs(v[k]) for k in ['Ldot','gtorque','external_torque','hybrid_torque'])
            assert abs(e-v['energy_residual'])<=D('1e-18')*e_scale+D('1e-80')
            assert abs(l-v['angular_residual'])<=D('1e-18')*l_scale+D('1e-80')
            rotation=v['intrinsic_L']*v['intrinsic_L']/(2*v['inertia'])
            assert abs(rotation-v['rotation_energy'])<=D('1e-18')*max(rotation,abs(v['rotation_energy']))+D('1e-80')
            U+=v['U'];K+=v['K_metric'];global_E_rate+=v['Kdot']+v['Udot']-v['gwork']-v['external_work'];global_L_rate+=v['Ldot']-v['gtorque']-v['external_torque']
            max_E=max(max_E,abs(e)/v['U']);max_L=max(max_L,abs(l)/(2*v['K_metric']))
    p=run/f'source_moments_rank{rank}.bin'
    with p.open('rb') as f:
        h=struct.unpack('<6Q',f.read(48));assert h[0:3]==(0x474653504d415331,1,6) and h[4:]==(16,18)
        values=np.frombuffer(f.read(),dtype=np.longdouble);assert len(values)==h[3]
    moments=values.reshape(-1,16);total_point_mass+=np.sum(moments[:,0],dtype=np.longdouble)
    curvature=json.loads((run/f'curvature_rank{rank}.json').read_text(),parse_float=D)
    assert curvature['chi_max']<=D('.002') and curvature['cell_mean_chi_max']<=D('.001')
assert len(seen)==71080 and abs(global_E_rate)/U<D('1e-8') and abs(global_L_rate)/(2*K)<D('1e-11')
assert abs(U-report['initial_U'])/U<D('1e-16') and abs(K-report['initial_K_metric'])/K<D('1e-16')
for p in sorted(run.iterdir()):
    if not(p.name.startswith('force_shadow_rank') and p.suffix=='.bin'):continue
    data=p.read_bytes();h=struct.unpack_from('<5Q',data)
    assert h[0:2]==(0x4746534653484431,1) and h[3:]==(16,36) and len(data)==40+72*16
    actual=np.frombuffer(data,dtype=np.longdouble,count=36,offset=40)
    cpu=np.frombuffer(data,dtype=np.longdouble,count=36,offset=40+36*16)
    error=np.max(np.abs(actual-cpu))/np.max(np.abs(cpu));assert error<np.longdouble('1e-11')
    shadow.append({'source':h[2],'relative_force_difference':str(error),'sha256':hashlib.sha256(data).hexdigest()})
assert len(shadow)==14
# Read already saved native initial data. No Voronoi reconstruction or IC
# regeneration; contrast generator kinetic energy with continuum-field K.
sys.path.insert(0,str(root.parent))
from audit_snapshot import read_snapshot
legacy_path=root.parent/'core_refinement/fine_integrated/khout.000001.dat'
p,meta=read_snapshot(legacy_path);assert len(p)==71080 and meta['t']==0
ld=np.longdouble;m=p['mass'].astype(ld);x=p['x'].astype(ld)-12;y=p['y'].astype(ld)-12;vx=p['vx'].astype(ld);vy=p['vy'].astype(ld)
legacy_K=np.sum(ld('.5')*m*(vx*vx+vy*vy),dtype=ld);legacy_U=np.sum(p['ie'].astype(ld),dtype=ld);legacy_L=np.sum(m*(x*vy-y*vx),dtype=ld);legacy_potential=-np.sum(m/np.sqrt(x*x+y*y+ld('.047')**2),dtype=ld)
out={
 'scope':'independent saved full-cell ledgers, force shadows and legacy initial-data comparison; no physical replay',
 'source_cells_audited':len(seen),'force_shadows_audited':len(shadow),
 'recomputed_cell_energy_rate_error_over_U_max':str(max_E),
 'recomputed_cell_angular_rate_error_over_2K_max':str(max_L),
 'recomputed_global_energy_rate_error_over_U':str(abs(global_E_rate)/U),
 'recomputed_global_angular_rate_error_over_2K':str(abs(global_L_rate)/(2*K)),
 'max_CPU_CUDA_force_shadow_difference':str(max(ld(s['relative_force_difference']) for s in shadow)),
 'actual_point_mass':str(total_point_mass),
 'legacy_generator_IC':{'mass':str(np.sum(m,dtype=ld)),'U':str(legacy_U),'K':str(legacy_K),'L':str(legacy_L),'potential':str(legacy_potential),'total_E':str(legacy_K+legacy_U+legacy_potential)},
 'new_continuum_material_IC':{k:str(report[v]) for k,v in [('U','initial_U'),('K_metric','initial_K_metric'),('L','initial_L'),('total_E','initial_total_E')]},
 'new_minus_legacy_K':str(ld(str(K))-legacy_K),
 'legacy_is_generator_quadrature_not_the_same_subcell_velocity_metric':True,
 'heat_adjusted_to_reconcile_K_difference':False,
 'audit_pass':True,'full_disk_time_evolution_verified':False,'shadow_witnesses':shadow}
(dest/'contract.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='shadow_witnesses'},indent=2))
