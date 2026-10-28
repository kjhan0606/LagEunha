"""Uniform-density necessary-condition audit on actual canonical face geometry.

This tests a proposed pure weight-remap mass law, not full ALE evolution.
Compute on lageunha. Never edits the production state or source data.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    source = Path('/gpfs/kjhan/LagForce/Codex_gradient_operator_20261017/operator')
    prior = Path('/gpfs/kjhan/LagForce/Codex_shared_flux_weight_response_20261017/response.npz')
    a = np.loadtxt(source/'rhs.dat')
    a = a[np.argsort(a[:,0])]
    n = len(a)
    z = np.load(prior)
    i,j = z['i'],z['j']
    f = np.loadtxt(source/'raw_faces.dat.common.rank0000')
    keys = f[:,0].astype(int)*n+f[:,1].astype(int)
    order = np.argsort(keys)
    keys,f = keys[order],f[order]
    idx = np.searchsorted(keys,i*n+j)
    assert np.array_equal(keys[idx],i*n+j)
    f = f[idx]
    rvec = f[:,2:4]
    distance = np.linalg.norm(rvec,axis=1)
    normal = rvec/distance[:,None]
    area = f[:,4]
    t = f[:,5:7]-.5*rvec
    t -= np.sum(t*normal,axis=1)[:,None]*normal
    velocity = a[:,3:5]
    geom_flux = area*(.5*np.sum((velocity[i]+velocity[j])*normal,axis=1)
                      +np.sum((velocity[i]-velocity[j])*t,axis=1)/distance)
    div = lambda flux:np.bincount(i,weights=flux,minlength=n)-np.bincount(j,weights=flux,minlength=n)
    vx = div(geom_flux)
    vw = div(area/(2*distance)*(z['weight_rate'][i]-z['weight_rate'][j]))
    volume = a[:,6]
    radius = np.linalg.norm(a[:,1:3]-3,axis=1)
    masks = {'inner_background':radius<.5,'disk_inner_edge':(radius>=.5)&(radius<.6),
             'disk_core':(radius>=.6)&(radius<1.9),'all_active':radius<2.5}
    rho = .01
    mass = rho*volume
    target = -vx/volume
    cases = []
    for scale in (0.,.25,1.,2.):
        vdot = vx+scale*vw
        mdot = rho*scale*vw
        fractional_density_rate = mdot/mass-vdot/volume
        error = fractional_density_rate-target
        assert np.max(abs(error[masks['all_active']]))<1e-12
        regions = {}
        for name,s in masks.items():
            rms = lambda v:float(np.sqrt(np.average(v[s]**2,weights=mass[s])))
            regions[name] = {'cells':int(s.sum()),'fractional_volume_rate_rms':rms(vdot/volume),
                             'fractional_mass_rate_rms':rms(mdot/mass),
                             'fractional_density_rate_rms':rms(fractional_density_rate),
                             'cancellation_identity_max_error':float(abs(error[s]).max())}
        cases.append({'weight_rate_scale':scale,'paired_global_mass_rate':float(mdot.sum()),'regions':regions})
    result = {'scope':'Necessary-condition diagnostic at actual t2 circular-velocity geometry, with manufactured uniform density. Full density jump and thermodynamics not tested.',
              'source':str(source),'weight_rate_source':str(prior),
              'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source/'rhs.dat',source/'raw_faces.dat.common.rank0000',prior)},
              'proposed_restricted_mass_law':'mdot_i = rho * Vdot_i_weight; fixed mass if weight rates vanish',
              'geometry':'Vdot_i = Vdot_i_generator + Vdot_i_weight',
              'identity':'rhodot_i/rho_i = mdot_i/m_i - Vdot_i/V_i = -Vdot_i_generator/V_i',
              'conclusion':'At a fixed instantaneous geometry and uniform density, mass exchange that is only conservative weight remapping cancels the weight contribution to density change. It cannot directly cancel generator-induced geometric compression by tuning weight rates.',
              'limits':['Not a proof that every weight-coupled mass scheme is impossible. Weights can alter future geometry and generator-induced volume response.',
                        'A more general mass law requires a specified physical density/transport interpretation, not arbitrary reassignment.',
                        'Does not establish pressure response, stability, or a complete conservative time integrator.',
                        'Paired internal faces used; unpaired distant frozen cap terms are outside active region.'],
              'cases':cases}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f: json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'output':str(args.output),'inner':[{**{'scale':c['weight_rate_scale']},**c['regions']['inner_background']} for c in cases]},indent=2))


if __name__=='__main__': main()
