"""Measure midpoint local pressure-weight work in a fresh weighted RK4 step."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weight_work_t65_20261027 import faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
RUN=ROOT/'gpu_energy_neutral_capacity64_v5_fresh_localwork_t0826_ada_v1'
RECEIPT=RUN/'audit_t0826.json'
OUTPUT=Path(__file__).with_name('KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T0826_2026-10-27.json')
N=55351


def rms(values):
    return float(np.sqrt(np.mean(values*values)))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--receipt',type=Path,default=RECEIPT)
    ap.add_argument('--output',type=Path,default=OUTPUT)
    ap.add_argument('--start',type=float,default=.825)
    ap.add_argument('--stop',type=float,default=.826)
    args=ap.parse_args()
    audit=json.loads(args.receipt.read_text())
    if (audit['time']!=args.stop or audit['parent_time']!=args.start or
            len(audit['checks'])!=21 or not all(audit['checks'].values())):
        raise ValueError('Weighted one-step native audit did not pass')
    start_path=Path(audit['parent'])
    end_path=Path(audit['endpoint'])
    if (digest(start_path)!=audit['parent_sha256'] or
            digest(end_path)!=audit['endpoint_sha256']):
        raise ValueError('Weighted one-step snapshots changed')
    start,meta0=ordered(start_path,N)
    end,meta1=ordered(end_path,N)
    dt=float(meta1['t'])-float(meta0['t'])
    if abs(dt-(args.stop-args.start))>1e-10 or not np.array_equal(start['mass'],end['mass']):
        raise ValueError('Unexpected interval or mass change')
    xy0=np.column_stack((start['x'].astype(float)-3.,start['y'].astype(float)-3.))
    xy1=np.column_stack((end['x'].astype(float)-3.,end['y'].astype(float)-3.))
    xy=.5*(xy0+xy1)
    w0=start['w2'].astype(float)
    w1=end['w2'].astype(float)
    weight=.5*(w0+w1)
    rate=(w1-w0)/dt
    pressure=.5*(start['pressure'].astype(float)+end['pressure'].astype(float))
    internal=.5*(start['ie'].astype(float)+end['ie'].astype(float))
    volume0=start['volume'].astype(float)
    volume1=end['volume'].astype(float)
    r0=np.linalg.norm(xy0,axis=1)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    active_flag=np.zeros(N,dtype=bool)
    active_flag[active]=True
    lookup=np.full(N,-1,dtype=int)
    lookup[active]=np.arange(len(active))
    area,meta,triangles,centers=power_areas(
        xy,weight,active,lookup,return_geometry=True)
    if meta['hidden_active'] or not np.all(area>0):
        raise ValueError('Invalid midpoint power diagram')
    midpoint_volume=.5*(volume0[active]+volume1[active])
    midpoint_area_error=float(np.max(np.abs(midpoint_volume/area-1)))
    if midpoint_area_error>1e-5:
        raise ValueError('Midpoint reconstructed area differs from native volumes')
    i,j,length=faces_touching_active(triangles,centers,active_flag)
    separation=np.linalg.norm(xy[j]-xy[i],axis=1)
    if np.any(separation<=0):
        raise ValueError('Coincident midpoint sites')
    conductance=length/(2*separation)
    face_weight_volume_rate=conductance*(rate[i]-rate[j])
    cell_weight_volume_rate=np.zeros(N)
    np.add.at(cell_weight_volume_rate,i,face_weight_volume_rate)
    np.add.at(cell_weight_volume_rate,j,-face_weight_volume_rate)
    local_power=pressure*cell_weight_volume_rate
    face_power=(pressure[i]-pressure[j])*face_weight_volume_rate
    local_sum=float(np.sum(local_power,dtype=np.longdouble))
    face_sum=float(np.sum(face_power,dtype=np.longdouble))
    face_abs=float(np.sum(np.abs(face_power),dtype=np.longdouble))
    identity_error=abs(local_sum-face_sum)
    if identity_error>1e-12*max(face_abs,1e-30):
        raise ValueError('Face/cell local pressure-work identity failed')
    fraction=dt*np.abs(local_power[active])/internal[active]
    modeled=dt*cell_weight_volume_rate[active]/volume0[active]
    observed=(volume1[active]-volume0[active])/volume0[active]
    result={
        'receipt':str(args.receipt),'receipt_sha256':digest(args.receipt),
        'start_sha256':audit['parent_sha256'],'end_sha256':audit['endpoint_sha256'],
        'dt':dt,'active_cells':len(active),'disk_touching_faces':len(i),
        'midpoint_area_max_relative_difference':midpoint_area_error,
        'face_cell_power_identity_abs_error':identity_error,
        'signed_weight_pressure_power_rate':face_sum,
        'absolute_face_weight_pressure_power_rate':face_abs,
        'signed_one_step_work_over_domain_internal_energy':dt*face_sum/float(np.sum(internal,dtype=np.longdouble)),
        'absolute_face_one_step_work_over_domain_internal_energy':dt*face_abs/float(np.sum(internal,dtype=np.longdouble)),
        'active_local_abs_weight_work_over_cell_internal_energy_p50':float(np.quantile(fraction,.5)),
        'active_local_abs_weight_work_over_cell_internal_energy_p99':float(np.quantile(fraction,.99)),
        'active_local_abs_weight_work_over_cell_internal_energy_max':float(np.max(fraction)),
        'weight_component_relative_volume_change_rms':rms(modeled),
        'actual_total_relative_volume_change_rms':rms(observed),
        'native_max_abs_global_projected_stage_work':audit['gpu_weight_max_abs_stage_work'],
        'scope':'Midpoint reconstruction over one native weighted RK4 step. It measures local pressure-weight work hidden by the stage-global projection, but not exact four-stage local energy closure or material-face flux.'}
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
