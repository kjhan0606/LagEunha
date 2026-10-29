"""All-core/disk frozen W0 pressure force and work, cross-checked with native K1."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import socket
import time
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_native_radial_20261026 import mean


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def describe(value,mass):
    if not np.isfinite(value).all():raise ValueError('Nonfinite physical sample')
    return {'mass_mean':mean(value,mass),'cell_median':float(np.median(value)),
            'cell_p90':float(np.percentile(value,90)),'cell_p99':float(np.percentile(value,99)),
            'cell_max':float(np.max(value))}


def audit(root):
    if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):
        raise RuntimeError('Run Voronoi reconstruction on a compute host')
    started=time.monotonic();cfg=json.loads((root/'probe.json').read_text())
    selected=json.loads((root/'audit.json').read_text())
    if not selected['passed'] or sha(cfg['initial'])!=cfg['initial_sha256'] or sha(cfg['endpoint'])!=cfg['endpoint_sha256']:
        raise ValueError('Accepted source or independent native selected-force audit failed')
    p0,_=ordered_global(cfg['initial']);ids=p0['u4if']&ID_MASK
    probe=root/'khdiag.step020790.wall_rhs_probe.dat';p,meta=ordered_global(probe,ids)
    if meta['t']!=10. or not np.array_equal(p0['mass'],p['mass']) or np.any(p['w2']!=0):
        raise ValueError('Requires exact fixed-mass W0 native t10 RHS state')
    x0=np.column_stack((p0['x']-3.,p0['y']-3.));r0=np.hypot(x0[:,0],x0[:,1])
    xy=np.column_stack((p['x']-3.,p['y']-3.));n=len(p)
    active=np.flatnonzero(r0<2.);flag=np.zeros(n,dtype=bool);flag[active]=True
    lookup=np.full(n,-1,dtype=int);lookup[active]=np.arange(len(active))
    area,geom,triangles,centers=power_areas(xy,np.zeros(n),active,lookup,return_geometry=True)
    if geom['hidden_active'] or np.any(area<=0):raise ValueError('Invalid interior Voronoi geometry')
    area_error=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if area_error>=1e-5:raise ValueError('Independent interior cell area does not match native')
    i,j,center,length=finite_faces_touching_active(triangles,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1)
    if np.any(dist<=0):raise ValueError('Coincident generators')
    normal=delta/dist[:,None];offset=center-.5*(xy[i]+xy[j])
    bi=length[:,None]*(.5*normal+offset/dist[:,None]);bj=length[:,None]*(.5*normal-offset/dist[:,None])
    pressure=p['pressure'].astype(float);dp=pressure[i]-pressure[j]
    force=np.zeros((n,2));volume_rate=np.zeros(n)
    ii=flag[i];jj=flag[j]
    np.add.at(force,i[ii],dp[ii,None]*bi[ii]);np.add.at(force,j[jj],dp[jj,None]*bj[jj])
    velocity=np.column_stack((p['vx'],p['vy']))
    normal_speed=.5*np.einsum('ij,ij->i',velocity[i]+velocity[j],normal)-np.einsum('ij,ij->i',velocity[j]-velocity[i],offset)/dist
    flux=length*normal_speed
    np.add.at(volume_rate,i[ii],flux[ii]);np.add.at(volume_rate,j[jj],-flux[jj])
    accel=force/p['mass'].astype(float)[:,None];udot=-pressure*volume_rate
    with (root/'rk4_stage_r0.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    expected={v['id'] for v in cfg['selected']}
    if len(rows)!=8 or {int(v['id']) for v in rows}!=expected:raise ValueError('Native sample incomplete')
    errors=[]
    for v in rows:
        uid=int(v['id']);idx=int(np.searchsorted(ids,uid));native=np.array([float(v['ax_hydro']),float(v['ay_hydro'])])
        errors.append({'id':uid,'acceleration_absolute':float(np.linalg.norm(accel[idx]-native)),
                       'thermal_work_absolute':abs(float(udot[idx])-float(v['Udot']))})
    if any(v['acceleration_absolute']>=1e-9 or v['thermal_work_absolute']>=1e-11 for v in errors):
        raise ValueError('All-cell reconstruction did not reproduce selected native force/work')
    radius=np.hypot(xy[:,0],xy[:,1]);omega=(radius*radius+.01175**2)**(-.75)
    gravity=radius*(radius*radius+.01175**2)**(-1.5);amag=np.linalg.norm(accel,axis=1)
    ratio=np.divide(amag,gravity,out=np.full(n,np.nan),where=gravity>1e-12)
    thermal=np.abs(udot)/(p['ie'].astype(float)*omega);mass=p['mass'].astype(float)
    regions=[]
    for lo,hi in [(0.,.1),(.1,.5),(.5,.75),(.75,1.),(1.,2.),(.5,2.)]:
        mask=(r0>=lo)&(r0<hi);valid=mask&np.isfinite(ratio)
        if np.count_nonzero(valid)!=np.count_nonzero(mask):raise ValueError('Zero gravity makes ratio undefined')
        regions.append({'rlo':lo,'rhi':hi,'cells':int(mask.sum()),
                        'pressure_acceleration':describe(amag[mask],mass[mask]),
                        'pressure_acceleration_over_gravity':describe(ratio[mask],mass[mask]),
                        'thermal_fraction_per_dynamical_time':describe(thermal[mask],mass[mask])})
    return {'passed':True,'time':meta['t'],'active_cells':len(active),'faces_touching_active':len(i),
            'native_volume_max_relative_error':area_error,'native_selected_force_work_errors':errors,'regions':regions,
            'analysis_seconds':time.monotonic()-started,
            'source_hashes':{str(q):sha(q) for q in [Path(cfg['initial']),Path(cfg['endpoint']),probe,root/'rk4_stage_r0.tsv',root/'audit.json',Path(__file__),Path(__file__).with_name('laguerre_nonlinear_volume_probe_t6_20261026.py'),Path(__file__).with_name('audit_kepler_bounded_pressure_force_t10_20261027.py')]},
            'scope':'One frozen, already core-damaged t10 native state; all original r0<2 IDs, reconstructed W0 geometry. Ratios measure instantaneous pressure acceleration and |Udot|/(U Omega). Compression work is not entropy production. A pointwise snapshot does not prove the temporal variation is slow, establish an integrator error order, or represent a core-accurate trajectory.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=audit(args.root)
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2))
