"""Estimate weight-induced conservative mass remap on native one-step faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N=55351


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('case',choices=('aligned','jitter'))
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    docs=Path(__file__).resolve().parent
    folder=docs/'KEPLER_Q2_NATIVE_WEIGHTED_T6501_2026-10-27'
    weighted=json.loads((folder/f'{args.case}_weighted_audit.json').read_text())
    w0=json.loads((folder/f'{args.case}_W0_audit.json').read_text())
    initial_name=('KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json' if args.case=='aligned'
                  else 'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json')
    initial_receipt=json.loads((docs/initial_name).read_text())
    if not all(weighted['checks'].values()) or not all(w0['checks'].values()):
        raise ValueError('Expected accepted native one-step receipts')
    if weighted['restart_sha256']!=w0['restart_sha256']:
        raise ValueError('W0 and weighted native starts do not match')
    for receipt in (weighted,w0):
        for key in ('restart','endpoint'):
            if digest(Path(receipt[key]))!=receipt[key+'_sha256']:
                raise ValueError(f'Changed {key} snapshot')
    if digest(Path(initial_receipt['initial']))!=initial_receipt['initial_sha256']:
        raise ValueError('Changed t=0 initial state')
    dt=float(weighted['endpoint_time'])-float(weighted['start_time'])
    if abs(dt-.001)>1e-10 or abs(float(w0['endpoint_time'])-float(weighted['endpoint_time']))>1e-10:
        raise ValueError('Expected matching t=6.5→6.501 interval')
    initial,_=ordered(Path(initial_receipt['initial']),N)
    start,_=ordered(Path(weighted['restart']),N)
    end,_=ordered(Path(weighted['endpoint']),N)
    base,_=ordered(Path(w0['endpoint']),N)
    if not (np.array_equal(start['mass'],end['mass']) and np.array_equal(start['mass'],base['mass'])):
        raise ValueError('Expected exact fixed native masses')
    r0=np.hypot(initial['x'].astype(float)-3.,initial['y'].astype(float)-3.)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    flag=np.zeros(N,dtype=bool);flag[active]=True
    lookup=np.full(N,-1,dtype=int);lookup[active]=np.arange(len(active))
    xy=.5*np.column_stack((start['x'].astype(float)+end['x'].astype(float)-6.,
                           start['y'].astype(float)+end['y'].astype(float)-6.))
    weight=.5*(start['w2'].astype(float)+end['w2'].astype(float))
    rate=(end['w2'].astype(float)-start['w2'].astype(float))/dt
    area,meta,triangles,centers=power_areas(xy,weight,active,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):
        raise ValueError('Invalid midpoint power geometry')
    i,j,c,length=finite_faces_touching_active(triangles,centers,flag)
    distance=np.linalg.norm(xy[j]-xy[i],axis=1)
    if np.any(distance<=0):
        raise ValueError('Duplicate generator pair')
    density=.5*(start['den'].astype(float)+end['den'].astype(float))
    face_mass=dt*.5*(density[i]+density[j])*length*(rate[i]-rate[j])/(2*distance)
    dm=np.zeros(N)
    np.add.at(dm,i,face_mass)
    np.add.at(dm,j,-face_mass)
    mass=start['mass'].astype(float)
    vstart=start['volume'].astype(float)
    vend=end['volume'].astype(float)
    vbase=base['volume'].astype(float)
    if np.any(mass+dm<=0):
        raise ValueError('Hypothetical remap gives nonpositive mass')
    rho0=mass/vstart
    rho_weighted=mass/vend
    rho_remap=(mass+dm)/vend
    rho_w0=mass/vbase
    w=mass[active]
    def mean_abs(a):
        return float(np.average(np.abs(a[active]),weights=w))
    drift_weighted=rho_weighted/rho0-1.
    drift_remap=rho_remap/rho0-1.
    drift_w0=rho_w0/rho0-1.
    gross=float(np.sum(np.abs(face_mass),dtype=np.longdouble))
    imbalance=float(np.sum(dm,dtype=np.longdouble))
    midpoint_volume=.5*(vstart[active]+vend[active])
    area_error=float(np.max(np.abs(area/midpoint_volume-1.)))
    checks={'pairwise_mass_conservation':abs(imbalance)<1e-12*max(gross,1e-30),
            'positive_remapped_mass':bool(np.all(mass+dm>0)),
            'midpoint_power_area':area_error<1e-4}
    result={'case':args.case,'start_sha256':weighted['restart_sha256'],
            'weighted_endpoint_sha256':weighted['endpoint_sha256'],
            'w0_endpoint_sha256':w0['endpoint_sha256'],'dt':dt,'disk_sites':len(active),
            'faces_touching_disk':len(i),'faces_disk_disk':int(np.count_nonzero(flag[i]&flag[j])),
            'midpoint_area_error_max':area_error,'gross_pair_mass_transfer':gross,
            'sum_all_cell_mass_change':imbalance,
            'disk_cell_mass_fraction_p99':float(np.quantile(np.abs(dm[active])/mass[active],.99)),
            'disk_mean_abs_density_change_w0':mean_abs(drift_w0),
            'disk_mean_abs_density_change_weighted_fixed_mass':mean_abs(drift_weighted),
            'disk_mean_abs_density_change_weighted_with_weight_remap':mean_abs(drift_remap),
            'disk_mean_abs_remap_minus_w0_density_change':mean_abs(drift_remap-drift_w0),
            'checks':{k:bool(v) for k,v in checks.items()},
            'scope':'Midpoint, arithmetic-density estimate of mass swept solely by weight-induced '
                    'face motion. This is not a native RK4 flux, a full ALE Riemann solve, '
                    'or an energy/momentum/spin-conservative evolution.'}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:result[k] for k in ('case','disk_mean_abs_density_change_w0',
          'disk_mean_abs_density_change_weighted_fixed_mass',
          'disk_mean_abs_density_change_weighted_with_weight_remap','checks')},indent=2))
    if not all(checks.values()):
        raise SystemExit('Weight-only remap audit failed')


if __name__=='__main__':
    main()
