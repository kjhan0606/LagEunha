"""Static spatial consistency of the actual smooth owner operator in 2D.

Repeated irregular motif preserves relative disorder under refinement. Tests
smooth fixed physical wavelengths and a grid-scale alternating mode separately.
This is not Euler/Kepler evolution or a bound for arbitrary mesh distortion.
"""
import json
import numpy as np
from scipy.spatial import Voronoi,cKDTree
from scipy.sparse import csr_matrix
from smooth_owner_memory import build_row

rng=np.random.default_rng(20261018)
motif=rng.uniform(-.2,.2,(4,4,2))
# Compute exact periodic motif Voronoi areas, reused under geometric scaling.
base=np.array([(i+.5+motif[i,j,0],j+.5+motif[i,j,1])
               for i in range(4) for j in range(4)])
tiles=np.concatenate([base+np.array([i,j])*4
                      for i,j in [(0,0)]+[(i,j) for i in (-1,0,1)
                      for j in (-1,0,1) if (i,j)!=(0,0)]])
vor=Voronoi(tiles);areas=[]
for i in range(16):
    region=vor.regions[vor.point_region[i]];assert -1 not in region
    p=vor.vertices[region]
    areas.append(.5*abs(np.dot(p[:,0],np.roll(p[:,1],-1))-
                        np.dot(p[:,1],np.roll(p[:,0],-1))))
areas=np.array(areas).reshape(4,4)
assert np.all(areas>0) and abs(areas.sum()-16)<1e-12
results=[]
for ns in (16,32,64):
    ij=np.array([(i,j) for i in range(ns) for j in range(ns)])
    x=(ij+.5+motif[ij[:,0]%4,ij[:,1]%4])/ns
    mass=areas[ij[:,0]%4,ij[:,1]%4]/ns**2
    volume=mass.copy() # rho=1 at this static test state; fixed cell masses
    c=1+.2*np.cos(2*np.pi*x[:,0])
    tree=cKDTree(x,boxsize=1.);rows=[];cols=[];values=[];moment=0.;condition=0.
    for i in range(len(x)):
        ids,r,meta=build_row(x,i,tree,box=1.)
        rows.extend([i]*len(ids));cols.extend(ids)
        values.extend(np.sqrt(mass[i])*c[i]/np.sqrt(volume[i])*r)
        moment=max(moment,meta['moment_max']);condition=max(condition,meta['condition'])
    B=csr_matrix((values,(rows,cols)),shape=(len(x),len(x)))
    velocity={
        'smooth_shear':np.column_stack((np.sin(2*np.pi*x[:,1]),np.zeros(len(x)))),
        'smooth_compression':np.column_stack((np.sin(2*np.pi*x[:,0]),np.sin(2*np.pi*x[:,1]))),
        'grid_alternating':np.column_stack(((-1.)**ij[:,0],np.zeros(len(x))))}
    cases={}
    for name,v in velocity.items():
        dz=B@v;k=(B.T@dz)/mass[:,None]
        norm2=float(np.sum(mass[:,None]*v*v))
        omega2=float(np.sum(dz*dz)/norm2)
        cases[name]={'auxiliary_stiffness_rayleigh':omega2,
                     'stiffness_times_h_squared':omega2/ns**2,
                     'acceleration_derivative_mass_RMS':float(np.sqrt(np.sum(mass[:,None]*k*k)/mass.sum())),
                     'Bv_L2':float(np.linalg.norm(dz)),
                     'quadratic_form_identity_error':float(np.sum(mass[:,None]*v*k)-np.sum(dz*dz))}
    row={'nside':ns,'cells':len(x),'h':1/ns,'max_moment':moment,
         'max_condition':condition,'cases':cases}
    results.append(row)
    print('completed nside='+str(ns),flush=True,file=__import__('sys').stderr)
orders={}
for name in velocity:
    orders[name]={}
    for key in ('auxiliary_stiffness_rayleigh','acceleration_derivative_mass_RMS','Bv_L2'):
        vals=[row['cases'][name][key] for row in results]
        orders[name][key]=[float(np.log2(vals[i]/vals[i+1])) for i in range(2)]
checks={
    'quadratic_moments':all(row['max_moment']<1e-12 for row in results),
    'smooth_stiffness_decreases':all(min(orders[name]['auxiliary_stiffness_rayleigh'])>2
        for name in ('smooth_shear','smooth_compression')),
    'smooth_pointwise_force_norm_decreases':all(min(orders[name]['acceleration_derivative_mass_RMS'])>0
        for name in ('smooth_shear','smooth_compression')),
    'grid_restoring_scale_retained':all(row['cases']['grid_alternating']['stiffness_times_h_squared']>.01 for row in results)}
print(json.dumps({'scope':__doc__,'rows':results,'observed_orders':orders,'checks':checks,
    'limitations':'Static B^T B response with z initially zero; no guarantee for accumulated memory, moving geometry, shocks, arbitrary disorder or full hydro. This irregular motif is a specific refinement family, not a universal convergence proof.'},indent=2))
assert all(checks.values()),'Spatial consistency gate failed; inspect saved evidence.'
