"""Independently check selected native Q1 pressure force/work at exact t=10."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from scipy.spatial import cKDTree
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from audit_kepler_wall_native_band_probe_20261029 import polygon, reconstruct


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(root):
    if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
        raise RuntimeError('Run independent geometry on a compute host')
    cfg=json.loads((root/'probe.json').read_text());initial,_=ordered_global(cfg['endpoint'])
    ids=initial['u4if']&ID_MASK;probe=root/'khdiag.step020790.wall_rhs_probe.dat'
    p,meta=ordered_global(probe,ids)
    if (meta['t']!=10. or sha(cfg['endpoint'])!=cfg['endpoint_sha256']
            or (root/'exit_status.txt').read_text().strip()!='0'
            or not all(np.array_equal(initial[k],p[k]) for k in ['x','y','vx','vy','mass','ie','w2'])):
        raise ValueError('RHS probe did not retain the exact accepted state')
    with (root/'rk4_stage_r0.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    expected={row['id']:row for row in cfg['selected']}
    if (len(rows)!=8 or {int(row['id']) for row in rows}!=set(expected)
            or any(row['stage']!='k1_rhs' or row['rhs_valid']!='1' or int(row['step'])!=cfg['restart_step']
                   or float(row['t_start'])!=10. for row in rows)):
        raise ValueError('Unexpected native selected-ID or stage output')
    xy=np.column_stack((p['x'],p['y']));tree=cKDTree(xy);cases=[]
    for row in rows:
        uid=int(row['id']);i=int(np.searchsorted(ids,uid));poly,area,support=polygon(i,xy,tree)
        vertices=np.array([v[:2] for v in poly]);length=np.linalg.norm(np.roll(vertices,-1,axis=0)-vertices,axis=1)
        radius=float(np.max(np.linalg.norm(vertices,axis=1)));resolution=64*np.finfo(float).eps*max(1.,radius)
        omitted=float(np.sum(length[length<=resolution]));near=float(tree.query(xy[i],k=2)[0][1])
        pmax=float(np.max(p['pressure']));vmax=float(np.max(np.hypot(p['vx'],p['vy'])))
        acc_bound=2*pmax*omitted*(.5+2*radius/near)/float(p['mass'][i])
        work_bound=4*pmax*vmax*omitted*(1+2*radius/near)
        accel,edot,udot,closure,_=reconstruct(i,xy,p,poly,edge_resolution=resolution)
        native=np.array([float(row['ax_hydro']),float(row['ay_hydro'])]);native_edot=float(row['die_stored'])
        error=float(np.linalg.norm(accel-native));work_error=abs(edot-native_edot)
        x=float(row['x'])-3.;y=float(row['y'])-3.;r=float(np.hypot(x,y));eps=.01175
        omega=(r*r+eps*eps)**(-.75);gravity=r*(r*r+eps*eps)**(-1.5)
        pressure=float(np.linalg.norm(native));u=float(row['U_raw']);ud=float(row['Udot'])
        checks={'area_relative_below_1e_minus6':abs(area/float(p['volume'][i])-1)<1e-6,
                'independent_pressure_acceleration':error+acc_bound<1e-9*max(1.,pressure),
                'independent_total_energy_work':work_error+work_bound<1e-11,
                'all_site_support_and_closure':support<1e-12 and closure<1e-10}
        cases.append({'id':uid,'purpose':expected[uid]['purpose'],'initial_radius':expected[uid]['initial_radius'],
                      'current_radius':r,'native_pressure_acceleration':pressure,'central_gravity_acceleration':gravity,
                      'pressure_acceleration_over_gravity':pressure/gravity if gravity>1e-12 else None,
                      'absolute_thermal_rate_over_internal_energy':abs(ud)/u,
                      'thermal_fraction_per_local_dynamical_time':abs(ud)/(u*omega),
                      'independent_acceleration_absolute_error':error,'independent_work_absolute_error':work_error,
                      'omitted_roundoff_edge_acceleration_bound':acc_bound,'omitted_roundoff_edge_work_bound':work_bound,
                      'checks':checks,'passed':all(checks.values())})
    return {'passed':all(v['passed'] for v in cases),'time':meta['t'],'cases':cases,
            'source_hashes':{str(p):sha(p) for p in [Path(cfg['endpoint']),probe,root/'rk4_stage_r0.tsv',root/'log',Path(__file__)]},
            'scope':'Eight purposefully selected frozen-state cells, including extremes. Force/work are independently reconstructed from bounded Voronoi polygons. Thermal rates are instantaneous compressional work, not entropy production; no population inference, pressure-history attribution, or proof of slow/fast splitting accuracy.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:raise SystemExit(1)
