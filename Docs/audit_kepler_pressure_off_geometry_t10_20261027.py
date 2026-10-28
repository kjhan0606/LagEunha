"""Check whether pressure-off t10 density residual is inherited site geometry."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS=Path(__file__).resolve().parent
ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
AUDIT=ROOT/'pressure_force_off_jitter_to10_a10_v1/audit_t10_native.json'
ORBIT=DOCS/'KEPLER_Q2_JITTER_INHERITED_ORBIT_T10_2026-10-27.json'


def weighted_abs(value,weight):
    return float(np.sum(np.abs(value)*weight,dtype=np.longdouble)
                 /np.sum(weight,dtype=np.longdouble))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    audit=json.loads(AUDIT.read_text())
    orbit=json.loads(ORBIT.read_text())
    p0_path=Path(audit['initial']);p_path=Path(audit['endpoint'])
    if (audit['time']!=10. or digest(p0_path)!=audit['initial_sha256']
            or digest(p_path)!=audit['endpoint_sha256']
            or orbit['pressure_off_sha256']!=audit['endpoint_sha256']):
        raise ValueError('Source audits or snapshots changed')
    p0,_=ordered(p0_path,55351)
    p,_=ordered(p_path,55351)
    r0=np.hypot(p0['x'].astype(float)-3.,p0['y'].astype(float)-3.)
    selected=np.flatnonzero((r0>=.5)&(r0<2.))
    if len(selected)!=49038:
        raise ValueError('Unexpected disk sites')
    xy=np.column_stack((p['x'].astype(float),p['y'].astype(float)))
    vor=Voronoi(xy,qhull_options='Qbb Qc Qz')
    area=np.empty(len(selected))
    for k,i in enumerate(selected):
        region=vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f'Unbounded pressure-off disk cell {i}')
        v=vor.vertices[region]
        area[k]=.5*abs(np.dot(v[:,0],np.roll(v[:,1],-1))
                     -np.dot(v[:,1],np.roll(v[:,0],-1)))
    _,exact_area,safety=disk_geometry(p0,10.,selected)
    mass=p0['mass'][selected].astype(float)
    native_volume=p['volume'][selected].astype(float)
    geometric_vs_native=area/native_volume-1.
    inherited_residual=exact_area/area-1.
    native_residual=exact_area/native_volume-1.
    geometric_error=weighted_abs(geometric_vs_native,mass)
    predicted=weighted_abs(inherited_residual,mass)
    recorded=weighted_abs(native_residual,mass)
    if abs(recorded-audit['disk_native_minus_ballistic_abs_mass_weighted'])>1e-10:
        raise ValueError('Native residual differs from accepted auditor')
    result={'time':10.,'disk_sites':len(selected),
            'pressure_off_endpoint_sha256':audit['endpoint_sha256'],
            'inherited_orbit_receipt_sha256':digest(ORBIT),
            'independent_voronoi_vs_native_volume_abs_mass_weighted':geometric_error,
            'independent_voronoi_vs_native_volume_max_relative':float(np.max(np.abs(geometric_vs_native))),
            'inherited_orbit_predicted_density_residual_abs_mass_weighted':predicted,
            'native_pressure_off_density_residual_abs_mass_weighted':recorded,
            'difference_of_pointwise_residuals_abs_mass_weighted':weighted_abs(
                inherited_residual-native_residual,mass),
            'minimum_exact_orbit_image_safety_ratio':safety,
            'scope':'Independent ordinary Voronoi construction at pressure-off native site positions. With the point-force orbit receipt, this tests whether the late residual is inherited state/geometry, not an accepted hydro solution.'}
    with args.output.open('x') as out:
        json.dump(result,out,indent=2)
        out.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
