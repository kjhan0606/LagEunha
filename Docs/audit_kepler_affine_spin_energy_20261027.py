"""Measure reconstructed affine spin/strain using native velocity gradients."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_polygon_moment_native_faces_20261027 import oriented_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest

N=55351


def tensor_from_faces(a,b,owner,offset,counts,mass):
    starts=offset[:-1]
    anchor=np.column_stack((np.add.reduceat(a[:,0],starts)/counts,
                            np.add.reduceat(a[:,1],starts)/counts))
    u=anchor[owner]
    area=.5*((a[:,0]-u[:,0])*(b[:,1]-u[:,1])-(a[:,1]-u[:,1])*(b[:,0]-u[:,0]))
    if np.min(area)<-1e-14:raise ValueError('Inconsistent triangle orientation')
    total=np.add.reduceat(area,starts)
    if np.any(total<=0):raise ValueError('Nonpositive polygon area')
    centroid=np.column_stack((np.add.reduceat(area*(u[:,0]+a[:,0]+b[:,0])/3.,starts)/total,
                              np.add.reduceat(area*(u[:,1]+a[:,1]+b[:,1])/3.,starts)/total))
    u=u-centroid[owner];v=a-centroid[owner];w=b-centroid[owner]
    scale=mass/total
    qxx=scale*np.add.reduceat(area*(u[:,0]**2+v[:,0]**2+w[:,0]**2+
                                u[:,0]*v[:,0]+v[:,0]*w[:,0]+w[:,0]*u[:,0])/6.,starts)
    qyy=scale*np.add.reduceat(area*(u[:,1]**2+v[:,1]**2+w[:,1]**2+
                                u[:,1]*v[:,1]+v[:,1]*w[:,1]+w[:,1]*u[:,1])/6.,starts)
    qxy=scale*np.add.reduceat(area*((u[:,0]+v[:,0]+w[:,0])*(u[:,1]+v[:,1]+w[:,1])+
                              u[:,0]*u[:,1]+v[:,0]*v[:,1]+w[:,0]*w[:,1])/12.,starts)
    return total,centroid,qxx,qxy,qyy


def energies(axx,axy,ayx,ayy,qxx,qxy,qyy):
    inertia=qxx+qyy
    spin=ayx*qxx-axy*qyy+(ayy-axx)*qxy
    total=.5*((axx*axx+ayx*ayx)*qxx+2*(axx*axy+ayx*ayy)*qxy+
              (axy*axy+ayy*ayy)*qyy)
    rotation=spin*spin/(2*inertia)
    return spin,rotation,total,total-rotation


def main():
    ap=argparse.ArgumentParser();ap.add_argument('receipt',type=Path)
    ap.add_argument('--thermal-risk',action='store_true',
                    help='Quantify hypothetical energy reallocation positivity risks')
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    receipt=json.loads(args.receipt.read_text())
    if not all(receipt['checks'].values()):raise ValueError('Parent physics gate failed')
    for name in ('initial','endpoint'):
        if digest(Path(receipt[name]))!=receipt[name+'_sha256']:
            raise ValueError(f'Changed {name} checkpoint')
    p0,m0=ordered(Path(receipt['initial']),N);p,m=ordered(Path(receipt['endpoint']),N)
    if abs(m0['t'])>1e-12 or abs(m['t']-receipt['time'])>1e-10:
        raise ValueError('Checkpoint time mismatch')
    xy0=np.column_stack((p0['x'].astype(float)-3.,p0['y'].astype(float)-3.))
    xy=np.column_stack((p['x'].astype(float)-3.,p['y'].astype(float)-3.))
    r0=np.linalg.norm(xy0,axis=1)
    ids=np.flatnonzero((r0>=.5)&(r0<2.));lookup=np.full(N,-1,dtype=int);lookup[ids]=np.arange(len(ids))
    weights=p['w2'].astype(float)
    area,meta,triangles,centers=power_areas(xy,weights,ids,lookup,return_geometry=True)
    if meta['hidden_active'] or np.any(area<=0):raise ValueError('Invalid power diagram')
    owner,a,b,offset,counts=oriented_faces(triangles,centers,xy,ids,lookup)
    mass=p['mass'][ids].astype(float)
    area_tri,centroid,qxx,qxy,qyy=tensor_from_faces(a,b,owner,offset,counts,mass)
    inertia=qxx+qyy
    grad=[p['stress.'+name][ids].astype(float) for name in ('gUxx','gUxy','gUyx','gUyy')]
    if not all(np.all(np.isfinite(x)) for x in grad):raise ValueError('Nonfinite stored gradients')
    spin,rot,total,residual=energies(*grad,qxx,qxy,qyy)
    x=xy[ids,0];y=xy[ids,1];rad2=x*x+y*y;omega=(rad2+.01175**2)**(-.75)
    deriv=-1.5*(rad2+.01175**2)**(-1.75)
    exact=[-y*deriv*x,-omega-y*deriv*y,omega+x*deriv*x,x*deriv*y]
    graderr=np.sqrt(sum((g-e)**2 for g,e in zip(grad,exact)))
    gradscale=np.sqrt(sum(e*e for e in exact))
    bands=[]
    for lo,hi in ((.5,.75),(.75,1.),(1.,1.5),(1.5,2.)):
        mask=(r0[ids]>=lo)&(r0[ids]<hi)
        bands.append({'r_min':lo,'r_max':hi,'cells':int(mask.sum()),
                      'rotational_to_affine_energy':float(np.sum(rot[mask])/np.sum(total[mask]))})
    checks={'finite_positive_tensor':bool(np.all(np.isfinite(inertia)) and np.all(qxx>0) and
                                          np.all(qyy>0) and np.all(qxx*qyy-qxy*qxy>0)),
            'independent_area':bool(np.max(np.abs(area_tri/area-1.))<1e-9),
            'native_volume':bool(np.max(np.abs(area_tri/p['volume'][ids].astype(float)-1.))<1e-5),
            'positive_affine_energy':bool(np.all(total>0)),
            'spin_energy_bounded':bool(np.min(residual)>=-1e-12*np.max(total)),
            'finite_gradients':True}
    result={'receipt':str(args.receipt.resolve()),'receipt_sha256':digest(args.receipt),
            'initial_sha256':receipt['initial_sha256'],'endpoint_sha256':receipt['endpoint_sha256'],
            'time':receipt['time'],'weighted':bool(np.any(weights!=0)),
            'disk_sites':len(ids),'directed_faces':len(owner),
            'stored_gradient_relative_rms_vs_analytic':float(np.sqrt(np.sum(graderr**2)/np.sum(gradscale**2))),
            'total_affine_energy':float(np.sum(total,dtype=np.longdouble)),
            'rigid_spin_energy':float(np.sum(rot,dtype=np.longdouble)),
            'residual_affine_energy':float(np.sum(residual,dtype=np.longdouble)),
            'rigid_spin_energy_fraction':float(np.sum(rot)/np.sum(total)),
            'radial_bands':bands,'checks':checks,
            'scope':'Postprocessing of saved native velocity gradients and independently '
                    'reconstructed uniform-density polygon moments. Rigid spin is a '
                    'projection of the affine velocity field, not an evolved native '
                    'spin variable or a measured irreversible heating rate. Gradient '
                    'stage timing and nonuniform subcell density remain limitations.'}
    if args.thermal_risk:
        thermal=p['ie'][ids].astype(float)
        if not np.all(np.isfinite(thermal)&(thermal>0)):
            raise ValueError('Nonpositive native thermal energy')
        mass_total=np.sum(mass,dtype=np.longdouble)
        ratio=rot/thermal
        bad=rot>=thermal
        full_bad=total>=thermal
        result['hypothetical_thermal_reallocation']={
            'native_thermal_energy_sum':float(np.sum(thermal,dtype=np.longdouble)),
            'rigid_rotation_over_thermal_sum':float(np.sum(rot,dtype=np.longdouble)/np.sum(thermal,dtype=np.longdouble)),
            'full_affine_over_thermal_sum':float(np.sum(total,dtype=np.longdouble)/np.sum(thermal,dtype=np.longdouble)),
            'rigid_rotation_over_thermal_median':float(np.median(ratio)),
            'rigid_rotation_over_thermal_p95':float(np.quantile(ratio,.95)),
            'cells_nonpositive_if_spin_subtracted':int(np.count_nonzero(bad)),
            'mass_fraction_nonpositive_if_spin_subtracted':float(np.sum(mass[bad],dtype=np.longdouble)/mass_total),
            'cells_nonpositive_if_affine_subtracted':int(np.count_nonzero(full_bad)),
            'mass_fraction_nonpositive_if_affine_subtracted':float(np.sum(mass[full_bad],dtype=np.longdouble)/mass_total),
            'scope':'Counterfactual subtraction from current native thermal U at fixed total '
                    'energy. The reconstructed affine energy is not a stored conserved '
                    'variable, so this is a positivity risk analysis, not a required '
                    'or physically valid initialization prescription.'}
    with args.output.open('x') as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
    print(json.dumps({'time':result['time'],'weighted':result['weighted'],
                      'fraction':result['rigid_spin_energy_fraction'],'checks':checks},indent=2))
    if not all(checks.values()):raise SystemExit('Affine spin diagnostic failed')


if __name__=='__main__':main()
