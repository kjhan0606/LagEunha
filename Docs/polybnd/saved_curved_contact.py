"""Preserved curved physical IC and independent exact-solution scoring.

The reference polygon is used only for validation/scoring, never by the density
provider, weight solve, or numerical RHS. The stored initial mass is authoritative.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import install_affine_shear_domain
from moment_aligned_density_audit import polygon_moments


def checked_json(path,digest):
    data=Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest()!=digest:
        raise ValueError('preserved curved input hash mismatch: '+str(path))
    return json.loads(data)


def initialize_saved(cfg):
    if cfg.get('kind')!='saved_curved_contact':raise ValueError('wrong saved IC kind')
    source=checked_json(cfg['initial_state'],cfg['initial_state_sha256'])
    reference=checked_json(cfg['reference'],cfg['reference_sha256'])
    n=cfg['n'];mass=np.asarray(source['mass'],dtype=float);state=np.asarray(source['state'],dtype=float)
    raw=np.asarray(reference['raw'],dtype=float)
    if source['t']!=0. or source['n']!=n or mass.shape!=(n*n,) or state.shape!=(20*n*n+4,):
        raise ValueError('invalid saved curved initial state identity')
    if not np.all(np.isfinite(state)) or not np.all(np.isfinite(mass)) or np.any(mass<=0):
        raise ValueError('invalid saved curved state values')
    if raw.shape!=(n*n,7) or not np.array_equal(mass,raw[:,0]):
        raise ValueError('saved curved reference mass mismatch')
    if not np.array_equal(state[:-4].reshape(-1,20)[:,3:9],raw[:,1:7]):
        raise ValueError('saved curved position moments mismatch')
    boundary=np.asarray(reference['reference_boundary'])
    if boundary.shape!=(128,2) or not Polygon(boundary).is_valid:
        raise ValueError('invalid physical reference polygon')
    install_affine_shear_domain()
    return mass.copy(),state.copy()


def advected_boundary(boundary,t):
    result=np.asarray(boundary,dtype=float).copy()
    result[:,0]+=t*(.2+1.1*result[:,1])
    return result


def reference_raw(cells,boundary,t):
    """Integrate the advected 128-gon, not a circular approximation."""
    high=Polygon(advected_boundary(boundary,t));exact=[];geometry=[]
    for cell in cells:
        geo=polygon_moments(cell);intersection=Polygon(cell).intersection(high)
        inside=np.zeros(7) if intersection.area<1e-28 else polygon_moments(np.asarray(intersection.exterior.coords[:-1]))
        geometry.append(geo);exact.append(geo+99*inside)
    return np.asarray(exact),np.asarray(geometry),high


def score(cp,cfg,label='final'):
    import moment_ale_evolve as model
    from predictive_evolution_batch import density_options
    from predictive_boundary_density import BoundaryDensity
    from face_density_selection_audit import one_sided_interval
    from curved_contact_rhs import affine_state_from_position_moments
    original,_=initialize_saved(cfg);mass=np.asarray(cp['mass']);np.testing.assert_array_equal(mass,original)
    reference=checked_json(cfg['reference'],cfg['reference_sha256']);t=float(cp['t'])
    state=np.asarray(cp['state']);a=state[:-4].reshape(-1,20)
    cells,edges=model.geometry(a[:,:2],a[:,2],t)
    exact,geo,high=reference_raw(cells,reference['reference_boundary'],t)
    _,exact_state,_=affine_state_from_position_moments(state,exact,geo)
    eraw=exact_state[:-4].reshape(-1,20)
    provider=BoundaryDensity(**density_options(cfg,True))
    _,closure=provider(cells,edges,np.column_stack([mass,a[:,3:]]))
    error=0.;face_error=0.;face_norm=0.
    for parts,rhos in provider.last_cells:
        for vertices,rho in zip(parts,rhos):
            poly=Polygon(vertices);inside=poly.intersection(high).area
            error+=abs(rho-100.)*inside+abs(rho-1.)*max(0.,poly.area-inside)
    boundary=advected_boundary(reference['reference_boundary'],t)
    for edge_id,(i,j,start,end,normal,length) in enumerate(edges):
        for cell,outward in ((i,normal),(j,-normal)):
            if cell<0:continue
            interval=one_sided_interval(boundary,start,end,-outward)
            profile=provider.last_profiles[(cell,edge_id)];S=len(profile)
            for k,rho in enumerate(profile):
                inside=0. if interval is None else max(0.,min((k+1)/S,interval[1])-max(k/S,interval[0]))
                outside=max(0.,1/S-inside)
                face_error+=length*(abs(rho-100.)*inside+abs(rho-1.)*outside)
                face_norm+=length*(100.*inside+outside)
    total=mass.sum();domain=model.domain(t)
    expected=polygon_moments(domain)+99*polygon_moments(boundary)
    partition=float(np.max(abs(exact.sum(0)-expected)/np.maximum(1.,abs(expected))))
    global_mass=float(abs(exact[:,0].sum()-total)/total)
    cell=float(abs(exact[:,0]-mass).sum()/total)
    valid=partition<1e-11 and global_mass<1e-11 and error/total+1e-8>=cell
    return dict(label=label,t=t,reference_valid=bool(valid),density_L1_relative=float(error/total),
        cell_mass_L1_relative=cell,max_cell_mass_relative=float(np.max(abs(exact[:,0]-mass)/mass)),
        q_L1_over_total_mass=float(abs(exact[:,1:3]-a[:,3:5]).sum()/total),
        Q_L1_over_total_mass=float(abs(exact[:,3:7]-a[:,5:9]).sum()/total),
        momentum_L1_over_total_mass=float(abs(eraw[:,9:11]-a[:,9:11]).sum()/total),
        mixed_L1_over_total_mass=float(abs(eraw[:,11:15]-a[:,11:15]).sum()/total),
        kinetic_L1_over_total_mass=float(abs(eraw[:,15:19]-a[:,15:19]).sum()/total),
        energy_L1_relative=float(abs(eraw[:,19]-a[:,19]).sum()/abs(eraw[:,19]).sum()),
        face_pointwise_L1_relative=float(face_error/face_norm),reference_partition_error=partition,
        exact_global_mass_relative_error=global_mass,closure=closure)
