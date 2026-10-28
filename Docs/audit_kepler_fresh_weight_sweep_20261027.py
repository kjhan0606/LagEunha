"""Estimate weight-only swept mass on an accepted fresh weighted step."""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_kepler_bounded_weight_work_t65_20261027 import faces_touching_active
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('receipt',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    audit=json.loads(args.receipt.read_text())
    if len(audit['checks'])!=21 or not all(audit['checks'].values()):
        raise ValueError('Native weighted step did not pass')
    for key in ('initial','parent','endpoint'):
        if digest(Path(audit[key]))!=audit[key+'_sha256']:raise ValueError(f'Changed {key}')
    initial,_=ordered(Path(audit['initial']),N)
    start,m0=ordered(Path(audit['parent']),N)
    end,m1=ordered(Path(audit['endpoint']),N)
    dt=float(m1['t'])-float(m0['t'])
    if abs(dt-.001)>1e-10 or not np.array_equal(start['mass'],end['mass']):
        raise ValueError('Expected one fixed-mass native step')
    r0=np.hypot(initial['x'].astype(float)-3.,initial['y'].astype(float)-3.)
    active=np.flatnonzero((r0>=.5)&(r0<2.))
    flag=np.zeros(N,dtype=bool);flag[active]=True
    lookup=np.full(N,-1,dtype=int);lookup[active]=np.arange(len(active))
    xy=.5*np.column_stack((start['x'].astype(float)+end['x'].astype(float)-6.,
                           start['y'].astype(float)+end['y'].astype(float)-6.))
    w=.5*(start['w2'].astype(float)+end['w2'].astype(float))
    rate=(end['w2'].astype(float)-start['w2'].astype(float))/dt
    area,meta,triangles,centers=power_areas(xy,w,active,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid power cells')
    i,j,length=faces_touching_active(triangles,centers,flag)
    d=np.linalg.norm(xy[j]-xy[i],axis=1)
    if np.any(d<=0):raise ValueError('Duplicate sites')
    sweep=dt*length*(rate[i]-rate[j])/(2*d)
    density=.5*(start['den'].astype(float)+end['den'].astype(float))
    mass=start['mass'].astype(float);v0=start['volume'].astype(float);v1=end['volume'].astype(float)
    weights=mass[active]
    def meanabs(x):return float(np.average(np.abs(x[active]),weights=weights))
    def remap(face_rho):
        face_mass=sweep*face_rho
        dm=np.zeros(N);np.add.at(dm,i,face_mass);np.add.at(dm,j,-face_mass)
        gross=float(np.sum(abs(face_mass),dtype=np.longdouble))
        imbalance=float(np.sum(dm,dtype=np.longdouble))
        if np.any(mass+dm<=0) or abs(imbalance)>1e-12*max(gross,1e-30):
            raise ValueError('Mass positivity or pair conservation failure')
        drift=((mass+dm)/v1)/(mass/v0)-1.
        return {'disk_mean_abs_density_change':meanabs(drift),
                'disk_mass_fraction_change_p99':float(np.quantile(abs(dm[active])/mass[active],.99)),
                'sum_cell_mass_change':imbalance,'gross_face_mass_transfer':gross},dm
    avg,dm_avg=remap(.5*(density[i]+density[j]))
    upwind,dm_upwind=remap(np.where(sweep>0.,density[j],density[i]))
    fixed=(mass/v1)/(mass/v0)-1.
    area_error=float(np.max(abs(area/(.5*(v0[active]+v1[active]))-1.)))
    if area_error>=1e-5:raise ValueError('Midpoint geometry mismatch')
    result={'initial_sha256':audit['initial_sha256'],'start_sha256':audit['parent_sha256'],
            'end_sha256':audit['endpoint_sha256'],'start_time':float(m0['t']),
            'end_time':float(m1['t']),'dt':dt,'active_cells':len(active),
            'faces_touching_active':len(i),'midpoint_area_relative_error_max':area_error,
            'fixed_mass_disk_mean_abs_density_change':meanabs(fixed),
            'arithmetic_density_remap':avg,'weight_sweep_upwind_remap':upwind,
            'upwind_minus_arithmetic_cell_mass_p99_fraction':float(np.quantile(abs(dm_upwind[active]-dm_avg[active])/mass[active],.99)),
            'checks':{'native_step':True,'midpoint_area':True,'mass_conservation_and_positivity':True},
            'scope':'First-order frozen-state midpoint estimate of weight-only swept mass on an already evolved fixed-mass trajectory. Not an evolved RK4 face flux, full ALE flux, or proof of a physically conservative weighted method.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()
