"""Sensitivity of the prior map at its roundoff branch, not an evolution.
Only inferred mean inputs are perturbed. No stored state/mass is modified.
"""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
from adaptive_coupled_mass import atomic
from neighbor_density_audit import interface_prior
from face_density_selection_audit import line_prior_mean
from affine_kinetic_mass import make_bins,solve,arrays,evaluate
from kinetic_factor_rhs import ALPHA
import moment_ale_evolve as model


def run():
    import sys
    smooth='--smooth' in sys.argv
    source=json.loads(Path('checkpoint.initial.json').read_text());cfg=source['config'];mass=np.asarray(source['mass'])
    setup('shear_contact',cfg['n'],cfg['amplitude']);record=json.loads(Path('step_80.json').read_text())['record']
    y=np.asarray(record['y']);a=y[:-4].reshape(-1,20);L=[np.asarray(f) for f in record['L']]
    polys,edges=model.geometry(a[:,:2],a[:,2],record['t']);raw=np.column_stack([mass,a[:,3:]])
    states=[model.decode(v,r) for v,r in zip(polys,raw)];gen=np.array([model.velocity(s,p) for s,p in zip(states,a[:,:2])])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,polys)]);order=np.argsort(means);cell=int(order[1])
    low=float(means.min());high=float(means.max());threshold=32*np.finfo(float).eps*max(1.,abs(low),abs(high))
    centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in polys]);neighbors=[set() for _ in mass]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].add(j);neighbors[j].add(i)
    S=64;rows=[];baseline_profiles=None;baseline_z=None
    for multiplier in (.5,.99,1.01,2.,4.,16.,100.):
        modified=means.copy();modified[cell]=low+threshold*multiplier
        priors=[interface_prior(i,v,centers,modified,neighbors,smooth_phase_floor=smooth)[1] for i,v in enumerate(polys)];profiles={}
        for edge_id,(i,j,va,vb,*_) in enumerate(edges):
            for owner in (i,j):
                if owner<0:continue
                profiles[(owner,edge_id)]=np.array([line_prior_mean(va+(vb-va)*k/S,va+(vb-va)*(k+1)/S,priors[owner],modified[owner],low,high) for k in range(S)])
        def trace(owner,edge_id,parameter):return profiles[(owner,edge_id)][min(S-1,int(parameter*S))]
        bins=make_bins(edges,a[:,:2],states,L,mass,gen,trace,S,model.boundary_velocity,model.velocity,ALPHA)
        z,_,_=solve(bins,len(mass))
        if baseline_profiles is None:baseline_profiles=profiles;baseline_z=z
        max_profile=max(float(np.max(abs(v-baseline_profiles[k]))) for k,v in profiles.items())
        residual,_=evaluate(arrays(bins),len(mass),baseline_z)
        row=dict(multiplier=multiplier,actual_gap=float(modified[cell]-low),flat=priors[cell]['flat'],blend=priors[cell].get('blend',0. if priors[cell]['flat'] else 1.),
            max_profile_difference=max_profile,max_weight_rate_difference=float(np.max(abs(z-baseline_z))),
            normalized_frozen_weight_mass_rate=1e-4*float(np.max(abs(residual)/mass)))
        rows.append(row);print(row,flush=True)
    atomic('threshold_probe_report.json',dict(completed=True,smooth_phase_floor=smooth,cell=cell,t=record['t'],threshold=threshold,rows=rows,
        note='Prior-map sensitivity only: mean inputs perturbed while geometry/velocity remain fixed. Not an admissibility test of perturbed full states, not an evolution, and not a change to saved masses. Comparison of frozen old weights reveals map discontinuity sensitivity.'))


if __name__=='__main__':run()
