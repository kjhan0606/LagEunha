"""Small static diagnostics for the subcell reference; no dynamics."""
import json
import numpy as np
from test_subcell_moments import setup_case
from subcell_moments import transfer, moments
from subcell_quadrature import quadrature, features, compress, subdivide

for kind in ('rotation', 'shear', 'kepler'):
    for ratio in (1., .01):
        fragments, _, receivers = setup_case(kind, ratio)
        fragments = subdivide(fragments)
        before = moments(fragments); groups, closure = transfer(fragments, receivers)
        after = moments([f for g in groups for f in g])
        spin = 0.; unresolved = 0.; constant_velocity_loss = 0.; max_moment_error = 0.
        max_points_before = 0; max_points_after = 0; phi_change = 0.; phi_change17 = 0.; max_points17 = 0
        for g in groups:
            m = moments(g); spin += m['spin']; unresolved += m['residual_energy']
            constant_velocity_loss += m['kinetic']-np.dot(m['momentum'], m['momentum'])/(2*m['mass'])
            p, w = quadrature(g); before_m = features(p)@w
            potential = lambda x: -1/np.linalg.norm(x, axis=1)
            p17, w17 = compress(p, w, potential)
            phi_before = np.sum(-w/np.linalg.norm(p[:, :2], axis=1))
            phi_change17 += np.dot(w17, potential(p17[:, :2]))-phi_before
            max_points17 = max(max_points17, len(w17))
            max_points_before = max(max_points_before, len(w))
            p, w = compress(p, w); max_points_after = max(max_points_after, len(w))
            max_moment_error = max(max_moment_error, float(np.max(abs(features(p)@w-before_m))))
            phi_change += np.sum(-w/np.linalg.norm(p[:, :2], axis=1))-phi_before
        print(json.dumps(dict(kind=kind, density_ratio=ratio, exact_coverage_error=closure,
              mass_relative=after['mass']/before['mass']-1,
              momentum_error=float(np.linalg.norm(after['momentum']-before['momentum'])),
              angular_error=after['angular']-before['angular'],
              kinetic_relative=after['kinetic']/before['kinetic']-1,
              internal_relative=after['internal']/before['internal']-1,
              retained_spin=spin, residual_energy=unresolved,
              constant_velocity_loss_over_internal=constant_velocity_loss/before['internal'],
              quadrature_points_before=max_points_before, quadrature_points_after=max_points_after,
              max_absolute_moment_error=max_moment_error,
              unpreserved_kepler_potential_change=phi_change,
              preserved_fixed_potential_change=phi_change17, quadrature_points17=max_points17)), flush=True)
